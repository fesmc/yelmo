"""The Yelmo parameter *schema*: the canonical defaults plus their metadata.

Built from ``input/yelmo_defaults.nml`` (names, default values, inline-comment
docs, section headers) enriched with the enum/range/order constraints.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from . import namelist as nl
from .constraints import Constraints, load_constraints
from .locate import find_defaults, find_src

# &yelmo keys that name the group a component is read from (default value =
# canonical group name). Used to resolve user-renamed groups, e.g.
# nml_ydyn = "ydyn_north".
COMPONENT_KEYS = [
    "nml_ytopo", "nml_ycalv", "nml_ydyn", "nml_ytill", "nml_ymat",
    "nml_ytherm", "nml_yhyd", "nml_masks", "nml_init_topo", "nml_data",
]

_SECTION_RE = re.compile(r"^\s*!\s*(?:={2,}|-{2,})\s*(.*?)\s*(?:={2,}|-{2,})\s*$")


@dataclass
class ParamSchema:
    group: str
    name: str
    default_raw: str
    type: str
    doc: str | None = None
    section: str | None = None
    enum: list | None = None         # union of allowed values, if any
    enum_constraints: list | None = None  # per-branch EnumConstraint list (conditional)

    @property
    def qualified(self) -> str:
        return f"{self.group}.{self.name}"


@dataclass
class Schema:
    nml: nl.Namelist
    params: dict = field(default_factory=dict)   # (group, name) -> ParamSchema
    constraints: Constraints = field(default_factory=Constraints)
    defaults_path: Path | None = None
    src_path: Path | None = None

    # -- group/component helpers ------------------------------------------- #
    @property
    def groups(self) -> list:
        return list(self.nml.groups.keys())

    def component_map(self) -> dict:
        """Map ``nml_*`` key -> default (canonical) group name."""
        yelmo = self.nml.get("yelmo")
        out = {}
        if yelmo is None:
            return out
        for key in COMPONENT_KEYS:
            p = yelmo.get(key)
            if p is not None:
                out[key] = nl.normalize(p.raw)
        return out

    # -- lookups ----------------------------------------------------------- #
    def get(self, group: str, name: str) -> ParamSchema | None:
        return self.params.get((group, name))

    def group_params(self, group: str) -> list:
        return [s for (g, _), s in self.params.items() if g == group]


def _section_label(text: str) -> str | None:
    m = _SECTION_RE.match(text)
    if m and m.group(1):
        return m.group(1)
    return None


def build_schema(defaults_path: str | Path | None = None,
                 src_path: str | Path | None = None) -> Schema:
    dp = find_defaults(defaults_path)
    nml = nl.parse_file(dp)
    sp = Path(src_path) if src_path else find_src(dp)

    constraints = load_constraints(sp)

    params: dict = {}
    for group in nml.groups.values():
        section: str | None = None
        for item in group.items:
            if isinstance(item, nl.Comment):
                lbl = _section_label(item.text)
                if lbl is not None:
                    section = lbl
                continue
            ecs = constraints.enum_for(group.name, item.name)
            enum_union = None
            if ecs:
                enum_union = []
                for ec in ecs:
                    for a in ec.allowed:
                        if a not in enum_union:
                            enum_union.append(a)
            params[(group.name, item.name)] = ParamSchema(
                group=group.name,
                name=item.name,
                default_raw=item.raw,
                type=nl.infer_type(item.raw),
                doc=item.comment,
                section=section,
                enum=enum_union,
                enum_constraints=ecs or None,
            )

    return Schema(nml=nml, params=params, constraints=constraints,
                  defaults_path=dp, src_path=sp)
