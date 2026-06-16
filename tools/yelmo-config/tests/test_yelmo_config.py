"""Tests for yelmo-config.

Run with:  python -m pytest tools/yelmo-config/tests
These exercise parsing/normalization/diff/check logic against a small inline
schema, plus a smoke test against the repo's real yelmo_defaults.nml if found.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from yelmo_config import namelist as nl
from yelmo_config.checking import check
from yelmo_config.constraints import (Constraints, OrderConstraint,
                                      RangeConstraint, EnumConstraint)
from yelmo_config.diffing import diff_all, diff_groups, empty_namelist
from yelmo_config.merge import complete_namelist, effective_group
from yelmo_config.resolve import parse_selector, resolve_groups
from yelmo_config.schema import Schema

DEFAULTS_TEXT = """\
! header comment
&yelmo
    domain    = "None"
    grid_name = "None"
    nml_ydyn  = "ydyn"
    nml_ytill = "ytill"
    cfl_max   = 0.1    ! Maximum value
/

&ydyn
    solver = "diva"    ! "sia", "diva"
    beta_q = 1.0
/

&ytill
    z0 = -300.0
    z1 = 200.0
/
"""


def make_schema() -> Schema:
    nml = nl.parse_string(DEFAULTS_TEXT)
    constraints = Constraints(
        enums={("ydyn", "solver"): [EnumConstraint("ydyn", "solver", ["sia", "diva"])]},
        ranges=[RangeConstraint("yelmo", "cfl_max", min=0.0, max=1.0,
                                min_inclusive=False, note="must be in (0,1]")],
        orders=[OrderConstraint("ytill", "z0", "lt", "z1", note="z0 < z1")],
    )
    from yelmo_config.schema import ParamSchema
    params = {}
    for g in nml.groups.values():
        for it in g.items:
            if isinstance(it, nl.Param):
                ecs = constraints.enum_for(g.name, it.name)
                allowed = ecs[0].allowed if ecs else None
                params[(g.name, it.name)] = ParamSchema(
                    g.name, it.name, it.raw, nl.infer_type(it.raw),
                    it.comment, None, allowed, ecs or None)
    return Schema(nml=nml, params=params, constraints=constraints)


# --------------------------------------------------------------------------- #
# namelist parsing / normalization
# --------------------------------------------------------------------------- #
def test_parse_roundtrip_preserves_comments():
    nml = nl.parse_string(DEFAULTS_TEXT)
    out = nl.dump(nml)
    assert "! header comment" in out
    assert "Maximum value" in out
    assert "&ydyn" in out and "&ytill" in out


def test_normalize_types():
    assert nl.normalize("1e3") == 1000.0
    assert nl.normalize('"x"') == "x"
    assert nl.normalize("True") is True
    assert nl.normalize(".false.") is False
    assert nl.normalize("11.7, 29.0, 57.0") == (11.7, 29.0, 57.0)
    assert nl.normalize('"a" "b"') == ("a", "b")


def test_inline_comment_not_split_inside_quotes():
    code, comment = nl._split_comment('opt = "-i minres ! foo"  ! real comment')
    assert "minres ! foo" in code
    assert comment.strip() == "real comment"


def test_infer_type():
    assert nl.infer_type("0.1") == "real"
    assert nl.infer_type("5") == "int"
    assert nl.infer_type('"x"') == "string"
    assert nl.infer_type("True") == "bool"
    assert nl.infer_type("1.0, 2.0") == "real[]"


# --------------------------------------------------------------------------- #
# group resolution / selectors
# --------------------------------------------------------------------------- #
def test_resolve_renamed_group():
    schema = make_schema()
    user = nl.parse_string('&yelmo\n nml_ydyn = "ydyn_north"\n/\n&ydyn_north\n solver="sia"\n/\n')
    mapping = resolve_groups(user, schema)
    assert mapping["ydyn_north"] == "ydyn"
    assert mapping["ytill"] == "ytill"   # default, untouched


def test_parse_selector():
    assert parse_selector("a.nml").group is None
    s = parse_selector("a.nml:ydyn_north")
    assert s.path == Path("a.nml") and s.group == "ydyn_north"
    assert parse_selector("dir/a.nml").group is None


# --------------------------------------------------------------------------- #
# merge / effective config
# --------------------------------------------------------------------------- #
def test_effective_group_overlays_defaults():
    schema = make_schema()
    user = nl.parse_string('&ydyn\n beta_q = 0.5\n/\n')
    eff = effective_group(user, "ydyn", schema)
    assert eff["beta_q"] == "0.5"          # overridden
    assert eff["solver"] == '"diva"'        # from defaults


def test_complete_namelist_applies_override():
    schema = make_schema()
    user = nl.parse_string('&ydyn\n solver = "sia"\n/\n')
    full = complete_namelist(user, schema)
    assert full["ydyn"].get("solver").raw == '"sia"'
    assert full["ytill"].get("z0").raw == "-300.0"   # default preserved


# --------------------------------------------------------------------------- #
# diff
# --------------------------------------------------------------------------- #
def test_diff_against_defaults():
    schema = make_schema()
    user = nl.parse_string('&ydyn\n solver = "sia"\n/\n')
    gds = diff_all(user, empty_namelist(), schema, "user", "defaults")
    changed = {d.name for gd in gds for d in gd.diffs}
    assert "solver" in changed


def test_group_scoped_cross_group_diff():
    schema = make_schema()
    a = nl.parse_string('&yelmo\n nml_ydyn="ydyn_north"\n/\n&ydyn_north\n solver="sia"\n/\n')
    b = nl.parse_string('&ydyn\n solver="diva"\n/\n')
    gd = diff_groups(a, "ydyn_north", b, "ydyn", schema, "A:ydyn_north", "B:ydyn")
    names = {d.name: d for d in gd.diffs}
    assert names["solver"].left == '"sia"' and names["solver"].right == '"diva"'


# --------------------------------------------------------------------------- #
# check
# --------------------------------------------------------------------------- #
def test_check_flags_violations():
    schema = make_schema()
    user = nl.parse_string(
        '&yelmo\n domain="Ant"\n grid_name="g"\n cfl_max=1.5\n/\n'
        '&ydyn\n solver="bogus"\n junk=1\n/\n'
        '&ytill\n z0=100.0\n z1=-50.0\n/\n')
    rep = check(user, schema)
    msgs = " ".join(i.message for i in rep.errors)
    assert "cfl_max" in msgs            # range
    assert "bogus" in msgs              # enum
    assert "junk" in msgs              # unknown
    assert "z0" in msgs                 # ordering
    assert not rep.ok


def test_check_passes_clean_file():
    schema = make_schema()
    user = nl.parse_string('&yelmo\n domain="Ant"\n grid_name="g"\n/\n&ydyn\n solver="sia"\n/\n')
    rep = check(user, schema)
    assert rep.ok


def test_check_reports_dummy_domain():
    schema = make_schema()
    user = nl.parse_string('&ydyn\n solver="sia"\n/\n')
    rep = check(user, schema)
    assert any(i.name == "domain" for i in rep.errors)


def test_passthrough_groups():
    schema = make_schema()
    user = nl.parse_string('&yelmo\n domain="A"\n grid_name="g"\n/\n&ctrl\n time_end=1.0\n/\n')
    rep = check(user, schema)
    assert "ctrl" in rep.passthrough_groups


# --------------------------------------------------------------------------- #
# smoke test against the real defaults file, if present
# --------------------------------------------------------------------------- #
def test_conditional_enum_extraction(tmp_path):
    from yelmo_config.constraints import extract_enums
    (tmp_path / "x.f90").write_text(
        "    subroutine ycalv_par_load(group_ycalv)\n"
        "        if (par%use_lsf) then\n"
        "            call yelmo_check_enum(group_ycalv,\"m\",par%m,\"a|b\")\n"
        "        else\n"
        "            call yelmo_check_enum(group_ycalv,\"m\",par%m,\"c|d\")\n"
        "        end if\n"
        "        call yelmo_check_enum(group_ycalv,\"g\",par%g,\"zero|x\")\n"
        "    end subroutine\n")
    enums = extract_enums(tmp_path)
    branches = enums[("ycalv", "m")]
    assert len(branches) == 2
    by_cond = {b.conditions[0]: b.allowed for b in branches}
    assert by_cond[("use_lsf", True)] == ["a", "b"]
    assert by_cond[("use_lsf", False)] == ["c", "d"]
    assert enums[("ycalv", "g")][0].conditions == []   # unconditional, after endif


def test_order_guard_holds():
    from yelmo_config.constraints import OrderConstraint
    eq = OrderConstraint("ymat", "a", "lt", "b", when_name="m", when_value="x")
    assert eq.guard_holds({"m": '"x"'})
    assert not eq.guard_holds({"m": '"y"'})
    sub = OrderConstraint("ymat", "a", "lt", "b", when_name="m", when_contains="-tracer")
    assert sub.guard_holds({"m": '"shear3D-tracer"'})
    assert not sub.guard_holds({"m": '"shear3D"'})
    assert not sub.guard_holds({})                  # guard param absent -> not enforced
    assert OrderConstraint("ymat", "a", "lt", "b").guard_holds({})  # unconditional


def test_enum_maybe_applies():
    from yelmo_config.constraints import EnumConstraint
    ec = EnumConstraint("ycalv", "m", ["a"], conditions=[("use_lsf", True)])
    assert ec.maybe_applies({"use_lsf": "True"})
    assert not ec.maybe_applies({"use_lsf": "False"})
    assert ec.maybe_applies({})                     # param absent -> lenient


def test_jsonable_and_completion():
    from yelmo_config.cli import _jsonable
    from yelmo_config.completion import script
    assert _jsonable("1e3") == 1000.0
    assert _jsonable('"x"') == "x"
    assert _jsonable("1.0, 2.0") == [1.0, 2.0]
    assert "complete -F" in script("bash")
    assert "#compdef" in script("zsh")
    assert "update" in script("bash") and "update" in script("zsh")
    with pytest.raises(ValueError):
        script("fish")


def test_update_dry_run(capsys):
    from yelmo_config import cli
    cli.main(["update", "--dry-run"])
    out = capsys.readouterr().out
    assert "DRY RUN" in out
    assert "pip install -U" in out
    assert "git+https://github.com/fesmc/yelmo#subdirectory=tools/yelmo-config" in out

    cli.main(["update", "dev", "--dry-run"])
    out = capsys.readouterr().out
    assert "git+https://github.com/fesmc/yelmo@dev#subdirectory=tools/yelmo-config" in out


def _run_cli(tmp_path, argv):
    """Run the CLI with a temp defaults file; return (exit_code, stdout)."""
    import io
    import contextlib
    from yelmo_config import cli
    dp = tmp_path / "yelmo_defaults.nml"
    dp.write_text(DEFAULTS_TEXT)
    buf = io.StringIO()
    code = 0
    with contextlib.redirect_stdout(buf):
        try:
            cli.main(["--defaults", str(dp), *argv])
        except SystemExit as e:
            code = e.code or 0
    return code, buf.getvalue()


def test_cli_list_compact(tmp_path):
    _, out = _run_cli(tmp_path, ["list", "ydyn", "--format", "compact"])
    assert 'ydyn.solver="diva"' in out
    assert "ydyn.beta_q=1.0" in out


def test_cli_check_json(tmp_path):
    # range (cfl_max) + unknown param: both work without the Fortran src/ present
    import json as _json
    user = tmp_path / "u.nml"
    user.write_text('&yelmo\n domain="A"\n grid_name="g"\n cfl_max=1.5\n/\n'
                    '&ydyn\n junk=1\n/\n')
    code, out = _run_cli(tmp_path, ["check", str(user), "--format", "json"])
    data = _json.loads(out)
    assert data["ok"] is False
    names = {i["name"] for i in data["issues"]}
    assert "cfl_max" in names and "junk" in names
    assert code == 1


def test_enums_json_roundtrip():
    from yelmo_config.constraints import (EnumConstraint, enums_from_json,
                                          enums_to_json)
    enums = {
        ("ydyn", "solver"): [EnumConstraint("ydyn", "solver", ["sia", "diva"])],
        ("ycalv", "m"): [
            EnumConstraint("ycalv", "m", ["a", "b"], conditions=[("use_lsf", True)]),
            EnumConstraint("ycalv", "m", ["c"], conditions=[("use_lsf", False)]),
        ],
    }
    back = enums_from_json(enums_to_json(enums))
    assert enums_to_json(back) == enums_to_json(enums)
    assert back[("ycalv", "m")][0].conditions == [("use_lsf", True)]


def test_bundled_files_present_and_loadable():
    # the committed snapshots must exist and parse
    from yelmo_config.locate import BUNDLED_DEFAULTS, BUNDLED_ENUMS
    from yelmo_config.constraints import load_bundled_enums
    from yelmo_config import namelist as nl
    assert BUNDLED_DEFAULTS.is_file() and BUNDLED_ENUMS.is_file()
    assert "ydyn" in nl.parse_file(BUNDLED_DEFAULTS).groups
    enums = load_bundled_enums(BUNDLED_ENUMS)
    assert enums[("ydyn", "solver")][0].allowed  # non-empty


def test_local_differs_warning(tmp_path, monkeypatch):
    # a local (auto-discovered) defaults file that differs from the bundle warns
    (tmp_path / "input").mkdir()
    (tmp_path / "input" / "yelmo_defaults.nml").write_text(
        '&yelmo\n domain="None"\n cfl_max=0.999\n/\n')
    monkeypatch.chdir(tmp_path)
    from yelmo_config.schema import build_schema
    schema = build_schema()                 # no explicit path -> source "local"
    assert schema.defaults_source == "local"
    assert any("differ" in w for w in schema.warnings)


def test_snapshot_writes_bundle(tmp_path, monkeypatch):
    import argparse
    from yelmo_config import cli
    from yelmo_config import locate
    from yelmo_config.constraints import load_bundled_enums
    # fake checkout
    (tmp_path / "input").mkdir()
    dp = tmp_path / "input" / "yelmo_defaults.nml"
    dp.write_text('&ycalv\n use_lsf=False\n calv_flt_method="vm-l19"\n/\n')
    (tmp_path / "src").mkdir()
    (tmp_path / "src" / "x.f90").write_text(
        '    subroutine ycalv_par_load(group_ycalv)\n'
        '        call yelmo_check_enum(group_ycalv,"calv_flt_method",par%x,"vm-l19|kill")\n'
        '    end subroutine\n')
    # redirect the bundle target into tmp
    out_def = tmp_path / "out" / "yelmo_defaults.nml"
    out_enm = tmp_path / "out" / "enums.json"
    monkeypatch.setattr(locate, "BUNDLED_DEFAULTS", out_def)
    monkeypatch.setattr(locate, "BUNDLED_ENUMS", out_enm)

    cli.cmd_snapshot(argparse.Namespace(defaults=str(dp), src=str(tmp_path / "src")))
    assert out_def.is_file() and out_enm.is_file()
    enums = load_bundled_enums(out_enm)
    assert enums[("ycalv", "calv_flt_method")][0].allowed == ["vm-l19", "kill"]


def test_real_defaults_smoke():
    from yelmo_config.locate import DEFAULTS_RELPATH
    here = Path(__file__).resolve()
    root = None
    for d in here.parents:
        if (d / DEFAULTS_RELPATH).is_file():
            root = d
            break
    if root is None:
        pytest.skip("real yelmo_defaults.nml not found")
    from yelmo_config.schema import build_schema
    schema = build_schema(root / DEFAULTS_RELPATH)
    assert "ydyn" in schema.groups
    # enum extracted from Fortran source, if src/ is present
    if schema.src_path:
        assert schema.get("ydyn", "solver").enum
