"""Command-line interface for yelmo-config."""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

from . import __version__
from . import namelist as nl
from .checking import check as run_check
from .completion import script as completion_script
from .diffing import diff_all, diff_groups, empty_namelist
from .merge import complete_namelist
from .resolve import parse_selector, resolve_groups
from .schema import build_schema


def _jsonable(raw: str):
    """Typed, JSON-serialisable view of a raw namelist value (tuples -> lists)."""
    v = nl.normalize(raw)
    return list(v) if isinstance(v, tuple) else v


def _emit_json(obj) -> None:
    print(json.dumps(obj, indent=2))


# --------------------------------------------------------------------------- #
# Tiny ANSI colour helper (no dependency)
# --------------------------------------------------------------------------- #
class C:
    enabled = True
    @classmethod
    def _w(cls, code, s):
        return f"\033[{code}m{s}\033[0m" if cls.enabled else s
    @classmethod
    def bold(cls, s): return cls._w("1", s)
    @classmethod
    def dim(cls, s): return cls._w("2", s)
    @classmethod
    def red(cls, s): return cls._w("31", s)
    @classmethod
    def green(cls, s): return cls._w("32", s)
    @classmethod
    def yellow(cls, s): return cls._w("33", s)
    @classmethod
    def cyan(cls, s): return cls._w("36", s)


def _load(path: Path) -> nl.Namelist:
    if not path.is_file():
        sys.exit(f"error: file not found: {path}")
    return nl.parse_file(path)


# --------------------------------------------------------------------------- #
# Commands
# --------------------------------------------------------------------------- #
def cmd_groups(args, schema):
    cmap = schema.component_map()
    if args.format == "json":
        _emit_json({
            "groups": [{"name": g, "n_params": len(schema.group_params(g))}
                       for g in schema.groups],
            "component_map": {k: str(v) for k, v in cmap.items()},
        })
        return
    if args.format == "compact":
        for g in schema.groups:
            print(g)
        return
    print(C.bold("Canonical library groups:"))
    for g in schema.groups:
        n = len(schema.group_params(g))
        print(f"  {C.cyan(g):<24} {C.dim(f'{n} parameters')}")
    print()
    print(C.bold("Component group mapping (set in &yelmo):"))
    for key, canonical in cmap.items():
        print(f"  {key:<16} -> &{canonical}")


def _param_dict(ps) -> dict:
    d = {
        "group": ps.group, "name": ps.name, "qualified": ps.qualified,
        "default": _jsonable(ps.default_raw), "default_raw": ps.default_raw,
        "type": ps.type, "doc": ps.doc, "section": ps.section, "enum": ps.enum,
    }
    if ps.enum_constraints and any(ec.conditional for ec in ps.enum_constraints):
        d["enum_branches"] = [
            {"allowed": ec.allowed,
             "when": [{"param": p, "value": v} for p, v in ec.conditions]}
            for ec in ps.enum_constraints
        ]
    return d


def _emit_params(params, fmt: str) -> None:
    """JSON / compact rendering shared by `list` and `search`."""
    if fmt == "json":
        _emit_json([_param_dict(ps) for ps in params])
    elif fmt == "compact":
        for ps in params:
            print(f"{ps.qualified}={ps.default_raw}")


def cmd_list(args, schema):
    groups = [args.group] if args.group else schema.groups
    if args.group and args.group not in schema.nml.groups:
        sys.exit(f"error: unknown group '{args.group}'. Try `yelmo-config groups`.")

    if args.format in ("json", "compact"):
        params = [ps for g in groups for ps in schema.group_params(g)]
        _emit_params(params, args.format)
        return

    for g in groups:
        if g not in schema.nml.groups:
            sys.exit(f"error: unknown group '{g}'. Try `yelmo-config groups`.")
        params = schema.group_params(g)
        print(C.bold(f"&{g}") + C.dim(f"  ({len(params)} parameters)"))
        last_section = object()
        for ps in params:
            if ps.section != last_section:
                last_section = ps.section
                if ps.section:
                    print(C.dim(f"  ── {ps.section} ──"))
            doc = f"  {C.dim(ps.doc)}" if ps.doc else ""
            enum = C.yellow(f"  [{('|').join(ps.enum)}]") if ps.enum else ""
            print(f"  {C.cyan(ps.name):<22} = {ps.default_raw:<18}{doc}{enum}")
        print()


def cmd_show(args, schema):
    token = args.param
    if "." not in token:
        sys.exit("error: specify a parameter as group.name, e.g. ydyn.solver")
    group, name = token.split(".", 1)
    ps = schema.get(group, name)
    if ps is None:
        sys.exit(f"error: no such parameter '{token}'. Try `yelmo-config search {name}`.")

    notes = [rc.note for rc in schema.constraints.ranges
             if rc.group == group and rc.name == name]
    notes += [oc.note for oc in schema.constraints.orders
              if oc.group == group and name in (oc.left, oc.right)]

    if args.format == "json":
        d = _param_dict(ps)
        d["constraints"] = notes
        _emit_json(d)
        return
    if args.format == "compact":
        print(f"{ps.qualified}={ps.default_raw}")
        return

    print(C.bold(ps.qualified))
    print(f"  default : {ps.default_raw}")
    print(f"  type    : {ps.type}")
    if ps.section:
        print(f"  section : {ps.section}")
    if ps.doc:
        print(f"  doc     : {ps.doc}")
    if ps.enum_constraints and any(ec.conditional for ec in ps.enum_constraints):
        print("  allowed :")
        for ec in ps.enum_constraints:
            if ec.conditions:
                cond = ", ".join(
                    f"{p}={'.true.' if v is True else '.false.' if v is False else repr(v)}"
                    for p, v in ec.conditions)
                print(f"    when {cond}: {', '.join(ec.allowed)}")
            else:
                print(f"    always: {', '.join(ec.allowed)}")
    elif ps.enum:
        print(f"  allowed : {', '.join(ps.enum)}")
    for rc in schema.constraints.ranges:
        if rc.group == group and rc.name == name:
            print(f"  range   : {rc.note}")
    for oc in schema.constraints.orders:
        if oc.group == group and oc.left == name:
            print(f"  order   : {oc.note}")
        if oc.group == group and oc.right == name:
            print(f"  order   : {oc.note}")


def cmd_search(args, schema):
    pat = re.compile(args.pattern, re.I)
    matches = [ps for ps in schema.params.values()
               if pat.search(f"{ps.name} {ps.doc or ''}") or pat.search(ps.qualified)]

    if args.format in ("json", "compact"):
        _emit_params(matches, args.format)
        return

    for ps in matches:
        doc = f"  {C.dim(ps.doc)}" if ps.doc else ""
        print(f"{C.cyan(ps.qualified):<34} = {ps.default_raw:<16}{doc}")
    if not matches:
        print(C.dim("no matches"))


def cmd_write(args, schema):
    user = _load(Path(args.user)) if args.user else None
    out_nml = complete_namelist(user, schema)

    if args.format == "json":
        obj = {g.name: {p.name: _jsonable(p.raw)
                        for p in g.items if isinstance(p, nl.Param)}
               for g in out_nml.groups.values()}
        text = json.dumps(obj, indent=2) + "\n"
    elif args.format == "compact":
        lines = [f"{g.name}.{p.name}={p.raw}"
                 for g in out_nml.groups.values()
                 for p in g.items if isinstance(p, nl.Param)]
        text = "\n".join(lines) + "\n"
    else:
        text = nl.dump(out_nml, align=not args.no_align)

    if args.output:
        Path(args.output).write_text(text)
        print(C.green(f"wrote {args.output}"))
    else:
        sys.stdout.write(text)


def _fmt_val(v):
    return v if v is not None else C.dim("(absent)")


def cmd_diff(args, schema):
    sel_a = parse_selector(args.a)
    nml_a = _load(sel_a.path)
    label_a = sel_a.path.name

    if args.b:
        sel_b = parse_selector(args.b)
        nml_b = _load(sel_b.path)
        label_b = sel_b.path.name
    else:
        sel_b = None
        nml_b = empty_namelist()
        label_b = "defaults"

    # group-scoped diff: at least one side names a group
    group_scoped = sel_a.group is not None or (sel_b and sel_b.group is not None)

    if group_scoped:
        ga = sel_a.group
        gb = (sel_b.group if sel_b else None) or ga  # default: same group name
        if ga is None:
            ga = gb
        if sel_b is None:
            # diffing one group against the defaults: use its canonical group
            gb = resolve_groups(nml_a, schema).get(ga, ga)
        la = f"{label_a}:{ga}"
        lb = f"{label_b}:{gb}"
        gds = [diff_groups(nml_a, ga, nml_b, gb, schema, la, lb, raw=args.raw)]
    else:
        gds = diff_all(nml_a, nml_b, schema, label_a, label_b, raw=args.raw)

    n = sum(gd.n_changed for gd in gds)

    if args.format == "json":
        _emit_json([
            {"left": gd.left_label, "right": gd.right_label,
             "diffs": [{"name": d.name, "left": d.left, "right": d.right,
                        "status": d.status} for d in gd.diffs]}
            for gd in gds if gd.diffs
        ])
    elif args.format == "compact":
        for gd in gds:
            grp = gd.left_label.rpartition(":")[2]
            for d in gd.diffs:
                left = d.left if d.left is not None else "(absent)"
                right = d.right if d.right is not None else "(absent)"
                print(f"{grp}.{d.name}: {left} -> {right}")
    else:
        for gd in gds:
            _print_group_diff(gd)
        print(C.dim(f"\n{n} difference(s)."))

    if args.exit_code and n:
        sys.exit(1)


def _print_group_diff(gd):
    if not gd.diffs:
        return
    print(C.bold(f"{gd.left_label}  vs  {gd.right_label}"))
    width = max(len(d.name) for d in gd.diffs)
    for d in gd.diffs:
        if d.status == "changed":
            print(f"  ~ {d.name:<{width}}  {C.red(d.left)}  ->  {C.green(d.right)}")
        elif d.status == "only_left":
            print(f"  - {d.name:<{width}}  {C.red(_fmt_val(d.left))}  "
                  f"{C.dim('(only in ' + gd.left_label + ')')}")
        else:
            print(f"  + {d.name:<{width}}  {C.green(_fmt_val(d.right))}  "
                  f"{C.dim('(only in ' + gd.right_label + ')')}")
    print()


def cmd_check(args, schema):
    user = _load(Path(args.file))
    root = Path(args.root) if args.root else schema.repo_root
    rep = run_check(user, schema, check_files=args.files, root=root)

    if args.format in ("json", "compact"):
        if args.format == "json":
            _emit_json({
                "file": str(args.file), "ok": rep.ok,
                "errors": len(rep.errors), "warnings": len(rep.warnings),
                "issues": [{"severity": i.severity, "group": i.group,
                            "name": i.name, "message": i.message}
                           for i in rep.issues],
                "passthrough_groups": rep.passthrough_groups,
            })
        else:
            for i in rep.issues:
                loc = f"{i.group}.{i.name}" if i.name else i.group
                print(f"{i.severity} {loc}: {i.message}")
        if not rep.ok:
            sys.exit(1)
        return

    for issue in rep.issues:
        loc = f"{issue.group}.{issue.name}" if issue.name else issue.group
        tag = C.red("error") if issue.severity == "error" else C.yellow("warning")
        print(f"  {tag}  {C.cyan(loc)}: {issue.message}")

    if rep.passthrough_groups:
        print(C.dim(f"\n  pass-through (non-library) groups: "
                    f"{', '.join(rep.passthrough_groups)}"))

    n_err, n_warn = len(rep.errors), len(rep.warnings)
    print()
    if rep.ok:
        msg = C.green("OK") if not n_warn else C.yellow(f"OK with {n_warn} warning(s)")
        print(f"  {msg} — {Path(args.file).name}")
    else:
        print(f"  {C.red('FAILED')} — {n_err} error(s), {n_warn} warning(s)")
        sys.exit(1)


# yelmo-config lives in the tools/yelmo-config subdirectory of the yelmo repo,
# so pip must be told to build from that subdirectory of the git checkout.
UPDATE_REPO = "git+https://github.com/fesmc/yelmo"
UPDATE_SUBDIR = "tools/yelmo-config"


def cmd_update(args) -> None:
    """`yelmo-config update [ref]` — self-update by reinstalling from the git
    repo with ``pip install -U``. With no ``ref``, pip pulls the default branch;
    pass a branch, tag, or commit SHA (e.g. ``yelmo-config update dev``) to
    install that ref instead."""
    url = UPDATE_REPO
    if args.ref:
        url = f"{url}@{args.ref}"
    url = f"{url}#subdirectory={UPDATE_SUBDIR}"
    cmd = [sys.executable, "-m", "pip", "install", "-U", url]

    if args.dry_run:
        print("DRY RUN — no changes will be made.")
        print(f"  would run: {' '.join(cmd)}")
        return

    print(f"yelmo-config update (currently {__version__})")
    print(f"  running: {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True)
    except (subprocess.CalledProcessError, OSError) as e:
        sys.exit(f"yelmo-config: self-update failed: {e}")


def cmd_snapshot(args) -> None:
    """`yelmo-config snapshot` — refresh the defaults + enum snapshots from a
    local Yelmo checkout into that checkout's yelmo_config/data/ directory. For
    maintainers: run from a checkout, then commit the updated data files.

    The destination is derived from the located checkout, not from this
    package's own data/ dir: under a non-editable install those differ, and
    writing to the latter would update site-packages while leaving the
    maintainer's tree untouched."""
    import json
    from . import namelist as nl
    from .constraints import enums_to_json, extract_enums
    from .locate import SNAPSHOT_RELPATH, find_local_defaults, find_src

    try:
        dp = find_local_defaults(args.defaults)
    except FileNotFoundError as e:
        sys.exit(f"yelmo-config snapshot: {e}")
    sp = Path(args.src) if args.src else find_src(dp)
    if not (sp and sp.is_dir()):
        sys.exit("yelmo-config snapshot: could not locate src/ for enum extraction "
                 "(pass --src). A full checkout is required.")

    out_dir = Path(dp).resolve().parent.parent / SNAPSHOT_RELPATH
    if not out_dir.is_dir():
        sys.exit(f"yelmo-config snapshot: no snapshot directory in the checkout: "
                 f"{out_dir}. A full checkout is required.")
    out_defaults = out_dir / "yelmo_defaults.nml"
    out_enums = out_dir / "enums.json"

    out_defaults.write_text(Path(dp).read_text())
    enums = extract_enums(Path(sp))
    out_enums.write_text(json.dumps(enums_to_json(enums), indent=2) + "\n")

    n_enum = sum(len(v) for v in enums.values())
    print(f"snapshot written to {out_dir}")
    print(f"  defaults : {dp}")
    print(f"  enums    : {sp}  ({n_enum} constraint(s) across {len(enums)} parameter(s))")


# --------------------------------------------------------------------------- #
# Argument parsing
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="yelmo-config",
        description="Discover, manage, compare and validate Yelmo parameter files.")
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--defaults", help="path to yelmo_defaults.nml (overrides auto-detection)")
    p.add_argument("--src", help="path to Yelmo src/ dir (for enum constraints)")
    p.add_argument("--no-color", action="store_true", help="disable coloured output")

    # --format lives on a shared parent so it may follow the subcommand
    # (e.g. `yelmo-config list ydyn --format json`).
    fmt = argparse.ArgumentParser(add_help=False)
    fmt.add_argument("--format", choices=["text", "json", "compact"], default="text",
                     help="output format (default: text)")

    sub = p.add_subparsers(dest="command", required=True)

    sp = sub.add_parser("groups", parents=[fmt],
                        help="list canonical groups and component mapping")
    sp.set_defaults(func=cmd_groups)

    sp = sub.add_parser("list", parents=[fmt],
                        help="list parameters (optionally for one group)")
    sp.add_argument("group", nargs="?", help="restrict to one group, e.g. ydyn")
    sp.set_defaults(func=cmd_list)

    sp = sub.add_parser("show", parents=[fmt],
                        help="show details for one parameter (group.name)")
    sp.add_argument("param", help="e.g. ydyn.solver")
    sp.set_defaults(func=cmd_show)

    sp = sub.add_parser("search", parents=[fmt],
                        help="search names and docs (regex, case-insensitive)")
    sp.add_argument("pattern")
    sp.set_defaults(func=cmd_search)

    sp = sub.add_parser("write", parents=[fmt],
                        help="write a complete par file from defaults (+overrides)")
    sp.add_argument("user", nargs="?", help="optional user par file whose overrides to apply")
    sp.add_argument("-o", "--output", help="output path (default: stdout)")
    sp.add_argument("--no-align", action="store_true", help="do not column-align output")
    sp.set_defaults(func=cmd_write)

    sp = sub.add_parser("diff", parents=[fmt],
                        help="compare two par files, or one against the defaults")
    sp.add_argument("a", help="file or file:group")
    sp.add_argument("b", nargs="?", help="file or file:group (default: the schema defaults)")
    sp.add_argument("--raw", action="store_true",
                    help="compare only params literally present (not the effective config)")
    sp.add_argument("--exit-code", action="store_true",
                    help="exit non-zero if any differences are found")
    sp.set_defaults(func=cmd_diff)

    sp = sub.add_parser("check", parents=[fmt],
                        help="validate a user par file against schema + constraints")
    sp.add_argument("file")
    sp.add_argument("--files", action="store_true",
                    help="also check that *_path inputs exist (needs domain/grid resolved)")
    sp.add_argument("--root", help="root dir to resolve relative input paths against")
    sp.set_defaults(func=cmd_check)

    sp = sub.add_parser("update",
                        help="self-update yelmo-config (pip install -U from its git repo)")
    sp.add_argument("ref", nargs="?", default=None,
                    help="git branch, tag, or commit SHA to install (e.g. 'dev'). "
                         "Default: the repo's default branch.")
    sp.add_argument("--dry-run", action="store_true",
                    help="print the pip command without running it")
    sp.set_defaults(func=cmd_update)

    sp = sub.add_parser("snapshot",
                        help="refresh the bundled defaults/enums from a checkout (maintainers)")
    sp.set_defaults(func=cmd_snapshot)

    sp = sub.add_parser("completion", help="emit a shell completion script (bash|zsh)")
    sp.add_argument("shell", choices=["bash", "zsh"])

    sp = sub.add_parser("_complete")   # hidden helper used by the completion scripts
    sp.add_argument("kind", choices=["groups", "params"])

    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    args = parser.parse_args(argv)

    if (args.no_color or os.environ.get("NO_COLOR") or not sys.stdout.isatty()
            or getattr(args, "format", "text") != "text"):
        C.enabled = False

    # these do not need a built schema (and work without a Yelmo checkout)
    if args.command == "completion":
        sys.stdout.write(completion_script(args.shell))
        return
    if args.command == "update":
        cmd_update(args)
        return
    if args.command == "snapshot":
        cmd_snapshot(args)
        return

    try:
        schema = build_schema(args.defaults, args.src)
    except FileNotFoundError as e:
        if args.command == "_complete":
            return   # silent: completion just offers nothing outside a checkout
        sys.exit(f"error: {e}")

    for w in schema.warnings:
        print(f"warning: {w}", file=sys.stderr)

    if args.command == "_complete":
        if args.kind == "groups":
            print("\n".join(schema.groups))
        else:
            print("\n".join(f"{g}.{n}" for (g, n) in schema.params))
        return

    args.func(args, schema)


if __name__ == "__main__":
    main()
