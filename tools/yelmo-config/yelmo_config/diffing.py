"""Compare two parameter sets (or one set against the defaults)."""

from __future__ import annotations

from dataclasses import dataclass

from . import namelist as nl
from .merge import effective_group
from .resolve import Selector, resolve_groups
from .schema import Schema


@dataclass
class ParamDiff:
    name: str
    left: str | None       # raw value on the left, or None if absent
    right: str | None
    status: str            # "changed" | "only_left" | "only_right"


@dataclass
class GroupDiff:
    left_label: str
    right_label: str
    diffs: list            # list[ParamDiff]

    @property
    def n_changed(self) -> int:
        return len(self.diffs)


def _diff_maps(left: dict, right: dict) -> list:
    out: list = []
    for name in sorted(set(left) | set(right)):
        lraw = left.get(name)
        rraw = right.get(name)
        if name in left and name in right:
            if nl.normalize(lraw) != nl.normalize(rraw):
                out.append(ParamDiff(name, lraw, rraw, "changed"))
        elif name in left:
            out.append(ParamDiff(name, lraw, None, "only_left"))
        else:
            out.append(ParamDiff(name, None, rraw, "only_right"))
    return out


def diff_groups(left_nml: nl.Namelist, left_group: str,
                right_nml: nl.Namelist, right_group: str,
                schema: Schema,
                left_label: str, right_label: str,
                raw: bool = False) -> GroupDiff:
    """Diff a single group on each side, aligned by parameter name."""
    lmap = effective_group(left_nml, left_group, schema, include_defaults=not raw)
    rmap = effective_group(right_nml, right_group, schema, include_defaults=not raw)
    return GroupDiff(left_label, right_label, _diff_maps(lmap, rmap))


def diff_all(left_nml: nl.Namelist, right_nml: nl.Namelist,
             schema: Schema,
             left_label: str, right_label: str,
             raw: bool = False) -> list:
    """Diff every canonical group's values. Returns list[GroupDiff]."""
    out: list = []
    lmap_actual = {c: a for a, c in resolve_groups(left_nml, schema).items()}
    rmap_actual = {c: a for a, c in resolve_groups(right_nml, schema).items()}

    for canonical in schema.groups:
        la = lmap_actual.get(canonical, canonical)
        ra = rmap_actual.get(canonical, canonical)
        lmap = effective_group(left_nml, la, schema, include_defaults=not raw)
        rmap = effective_group(right_nml, ra, schema, include_defaults=not raw)
        diffs = _diff_maps(lmap, rmap)
        if diffs:
            out.append(GroupDiff(f"{left_label}:{canonical}",
                                 f"{right_label}:{canonical}", diffs))
    return out


def empty_namelist() -> nl.Namelist:
    """A namelist with no overrides -> effective config equals pure defaults."""
    return nl.Namelist()
