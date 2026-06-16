"""Validate a user par file against the schema and constraints.

Checks performed (only over the canonical library groups; driver groups such as
``&ctrl`` / ``&set_*`` are reported once as informational pass-throughs):

  * unknown parameters   — names not declared in the defaults for that group
  * enum violations      — value outside the ``yelmo_check_enum`` allowed set
  * range violations     — scalar outside the curated numeric range
  * ordering violations  — inter-parameter ordering (e.g. sd_min < sd_max)
  * type mismatches      — e.g. a string where the default is numeric
  * missing file inputs  — ``*_path`` whose ``*_load`` is true but file absent
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from . import namelist as nl
from .merge import effective_all
from .resolve import resolve_groups
from .schema import Schema

# severity: "error" mirrors a Fortran `stop`; "warning" is advisory.
@dataclass
class Issue:
    severity: str
    group: str
    name: str | None
    message: str


@dataclass
class Report:
    issues: list = field(default_factory=list)
    passthrough_groups: list = field(default_factory=list)

    @property
    def errors(self) -> list:
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> list:
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def ok(self) -> bool:
        return not self.errors


def _type_compatible(default_type: str, value_raw: str) -> bool:
    vt = nl.infer_type(value_raw)
    if default_type == vt:
        return True
    numeric = {"int", "real"}
    if default_type in numeric and vt in numeric:
        return True   # int<->real is fine
    # arrays: compare base
    if default_type.endswith("[]") or vt.endswith("[]"):
        return default_type.rstrip("[]") == vt.rstrip("[]") or "mixed" in (default_type, vt)
    return False


def _subst_path(raw: str, domain: str | None, grid: str | None) -> str | None:
    s = nl.normalize(raw)
    if not isinstance(s, str):
        return None
    if domain:
        s = s.replace("{domain}", domain)
    if grid:
        s = s.replace("{grid_name}", grid)
    if "{" in s:
        return None  # unresolved template -> skip existence check
    return s


def check(user: nl.Namelist, schema: Schema, *,
          check_files: bool = False, root: Path | None = None) -> Report:
    rep = Report()
    mapping = resolve_groups(user, schema)             # actual -> canonical
    known_actual = set(mapping)

    # domain / grid for path substitution
    yelmo = user.get("yelmo")
    domain = grid = None
    if yelmo is not None:
        dp, gp = yelmo.get("domain"), yelmo.get("grid_name")
        domain = nl.normalize(dp.raw) if dp else None
        grid = nl.normalize(gp.raw) if gp else None

    # effective config (defaults + overrides) — needed to evaluate conditional
    # enum guards (e.g. calv_flt_method depends on use_lsf) and value checks.
    eff = effective_all(user, schema)   # canonical -> {name: raw}

    # --- per-group structural checks: unknown params, enums, types ---------- #
    for gname, group in user.groups.items():
        if gname not in known_actual:
            rep.passthrough_groups.append(gname)
            continue
        canonical = mapping[gname]

        for p in group.params().values():
            ps = schema.get(canonical, p.name)
            if ps is None:
                rep.issues.append(Issue(
                    "error", gname, p.name,
                    f"unknown parameter '{p.name}' (not declared in &{canonical} defaults)"))
                continue
            # type
            if not _type_compatible(ps.type, p.raw):
                rep.issues.append(Issue(
                    "warning", gname, p.name,
                    f"type mismatch: default is {ps.type}, got {nl.infer_type(p.raw)} ({p.raw})"))
            # enum (honour conditional guards against the effective config)
            ecs = schema.constraints.enum_for(canonical, p.name)
            if ecs:
                gmap = eff.get(canonical, {})
                acceptable: list = []
                for ec in ecs:
                    if ec.maybe_applies(gmap):
                        for a in ec.allowed:
                            if a not in acceptable:
                                acceptable.append(a)
                val = nl.normalize(p.raw)
                if str(val) not in acceptable:
                    rep.issues.append(Issue(
                        "error", gname, p.name,
                        f"invalid value '{val}' — allowed: {', '.join(acceptable)}"))

    actual_by_canon = {c: a for a, c in mapping.items()}

    def report_group(canonical: str) -> str:
        return actual_by_canon.get(canonical, canonical)

    for rc in schema.constraints.ranges:
        gmap = eff.get(rc.group)
        if not gmap or rc.name not in gmap:
            continue
        val = nl.normalize(gmap[rc.name])
        if rc.violates(val):
            rep.issues.append(Issue(
                "error", report_group(rc.group), rc.name,
                f"{rc.name} = {val} out of range — {rc.note}"))

    for oc in schema.constraints.orders:
        gmap = eff.get(oc.group)
        if not gmap or oc.left not in gmap or oc.right not in gmap:
            continue
        if not oc.guard_holds(gmap):
            continue
        lval, rval = nl.normalize(gmap[oc.left]), nl.normalize(gmap[oc.right])
        if not oc.satisfied(lval, rval):
            rep.issues.append(Issue(
                "error", report_group(oc.group), oc.left,
                f"{oc.left}={lval}, {oc.right}={rval} — {oc.note}"))

    # --- missing 'None' overrides (domain/grid still dummy) ---------------- #
    if domain == "None" or domain is None:
        rep.issues.append(Issue("error", "yelmo", "domain",
                                "domain is still 'None' — user par file must set it"))
    if grid == "None" or grid is None:
        rep.issues.append(Issue("error", "yelmo", "grid_name",
                                "grid_name is still 'None' — user par file must set it"))

    # --- file-existence checks (*_load true => *_path must exist) ----------- #
    if check_files:
        for canonical, gmap in eff.items():
            for name, raw in gmap.items():
                m = re.fullmatch(r"(.+)_load", name)
                if not m:
                    continue
                if nl.normalize(raw) is not True:
                    continue
                path_key = f"{m.group(1)}_path"
                if path_key not in gmap:
                    continue
                resolved = _subst_path(gmap[path_key], domain, grid)
                if resolved is None:
                    continue
                full = (root / resolved) if (root and not Path(resolved).is_absolute()) else Path(resolved)
                if not full.exists():
                    rep.issues.append(Issue(
                        "error", report_group(canonical), path_key,
                        f"{name}=True but file not found: {resolved}"))
    return rep
