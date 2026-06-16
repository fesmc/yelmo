"""Parameter constraints, assembled from two sources.

* **Enum** constraints are extracted from the ``yelmo_check_enum`` calls in the
  Fortran ``src/`` (single source of truth — the model itself enforces them).
* **Range** and **ordering** constraints are loaded from the curated
  ``data/constraints.toml`` (see that file for why they are not scraped).

Everything is keyed by canonical group name (``ytopo``, ``ydyn``, ...).
"""

from __future__ import annotations

import json
import re
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent / "data"


# --------------------------------------------------------------------------- #
# Constraint records
# --------------------------------------------------------------------------- #
@dataclass
class EnumConstraint:
    group: str
    name: str
    allowed: list           # list[str]
    # Guard conditions under which this allowed-set applies, as (param, value)
    # pairs that must ALL hold (value is True/False for a flag, or a string).
    # A param of None marks an unparseable guard. Empty => unconditional.
    conditions: list = field(default_factory=list)

    def maybe_applies(self, gmap: dict) -> bool:
        """True unless a known guard provably fails for the given config.

        Lenient by design: an unparseable guard, or one referencing a parameter
        absent from ``gmap``, does not exclude the branch — avoids false errors.
        """
        from . import namelist as _nl
        for param, value in self.conditions:
            if param is None or param not in gmap:
                return True
            if _nl.normalize(gmap[param]) != value:
                return False
        return True

    @property
    def conditional(self) -> bool:
        return bool(self.conditions)


@dataclass
class RangeConstraint:
    group: str
    name: str
    min: float | None = None
    max: float | None = None
    min_inclusive: bool = True
    max_inclusive: bool = True
    note: str = ""

    def violates(self, value) -> bool:
        try:
            v = float(value)
        except (TypeError, ValueError):
            return False
        if self.min is not None:
            if v < self.min or (not self.min_inclusive and v == self.min):
                return True
        if self.max is not None:
            if v > self.max or (not self.max_inclusive and v == self.max):
                return True
        return False


@dataclass
class OrderConstraint:
    group: str
    left: str
    op: str          # lt | le | gt | ge
    right: str
    note: str = ""
    when_name: str | None = None
    when_value: str | None = None       # guard holds when when_name == when_value
    when_contains: str | None = None    # guard holds when when_name contains this substring

    def guard_holds(self, gmap: dict) -> bool:
        """Whether this ordering check applies for the given config."""
        if self.when_name is None:
            return True
        from . import namelist as _nl
        raw = gmap.get(self.when_name)
        if raw is None:
            return False
        val = _nl.normalize(raw)
        if self.when_value is not None:
            return str(val) == self.when_value
        if self.when_contains is not None:
            return self.when_contains in str(val)
        return True

    def satisfied(self, lval, rval) -> bool:
        try:
            a, b = float(lval), float(rval)
        except (TypeError, ValueError):
            return True  # can't compare -> don't flag
        return {
            "lt": a < b, "le": a <= b, "gt": a > b, "ge": a >= b,
        }[self.op]


@dataclass
class Constraints:
    enums: dict = field(default_factory=dict)   # (group, name) -> list[EnumConstraint]
    ranges: list = field(default_factory=list)  # list[RangeConstraint]
    orders: list = field(default_factory=list)  # list[OrderConstraint]

    def enum_for(self, group: str, name: str) -> list:
        return self.enums.get((group, name), [])


# --------------------------------------------------------------------------- #
# Enum extraction from Fortran source
# --------------------------------------------------------------------------- #
# group argument -> canonical group name, for the cases where the call passes a
# bare `group` variable rather than a `group_<kind>` one. Resolved instead via
# the enclosing `subroutine <kind>_par_load` name; this map is the fallback.
_SUBROUTINE_RE = re.compile(r"^\s*subroutine\s+(\w+)", re.I)
_SUBROUTINE_GROUP_RE = re.compile(r"^\s*subroutine\s+(\w+?)_par_load\b", re.I)

# Block-structure / condition recognition for conditional enum guards.
_IF_THEN_RE = re.compile(r"^\s*if\s*\((.*)\)\s*then\s*$", re.I)
_ELSEIF_RE = re.compile(r"^\s*else\s*if\b", re.I)
_ELSE_RE = re.compile(r"^\s*else\s*$", re.I)
_ENDIF_RE = re.compile(r"^\s*end\s*if\s*$", re.I)


def _parse_cond(expr: str):
    """Parse a simple Fortran condition into a (param, value) guard.

    Recognises ``par%flag``, ``.not. par%flag`` and
    ``[trim(]par%name[)] == "value"``. Anything else -> (None, None).
    """
    e = expr.strip()
    m = re.fullmatch(r"\.not\.\s*par%(\w+)", e, re.I)
    if m:
        return (m.group(1), False)
    m = re.fullmatch(r"par%(\w+)", e, re.I)
    if m:
        return (m.group(1), True)
    m = re.fullmatch(r'(?:trim\()?\s*par%(\w+)\)?\s*==\s*["\']([^"\']*)["\']', e, re.I)
    if m:
        return (m.group(1), m.group(2))
    return (None, None)


def _negate(cond):
    param, value = cond
    if param is not None and isinstance(value, bool):
        return (param, not value)
    return (None, None)   # negation of a string-equality / unknown guard
_ENUM_CALL_RE = re.compile(
    r"""yelmo_check_enum\s*\(\s*
        (?P<group>[\w"']+)\s*,\s*       # group token: group_x | group | "literal"
        ["'](?P<name>\w+)["']\s*,\s*    # parameter name
        [^,]+,\s*                        # par%<member>  (ignored)
        ["'](?P<allowed>[^"']*)["']\s*   # pipe-separated allowed set
        \)""",
    re.VERBOSE,
)


def _join_continuations(text: str) -> list[tuple[int, str]]:
    """Join Fortran ``&`` continuation lines; return (lineno, joined) pairs."""
    out: list[tuple[int, str]] = []
    buf, start = "", 0
    for i, line in enumerate(text.splitlines(), start=1):
        code = line.split("!", 1)[0].rstrip()
        if buf:
            code = code.lstrip()
            if code.startswith("&"):
                code = code[1:]
        if code.endswith("&"):
            if not buf:
                start = i
            buf += code[:-1]
            continue
        if buf:
            out.append((start, buf + code))
            buf = ""
        else:
            out.append((i, code))
    if buf:
        out.append((start, buf))
    return out


def _resolve_group(token: str, current_sub_group: str | None) -> str | None:
    token = token.strip()
    if token.startswith(("'", '"')):
        return token.strip("'\"")
    if token.startswith("group_"):
        return token[len("group_"):]
    if token == "group":
        return current_sub_group
    return None


def extract_enums(src_dir: Path) -> dict:
    """Return ``{(group, name): [EnumConstraint, ...]}`` parsed from ``src_dir``.

    A parameter may have more than one entry when its check is guarded by an
    ``if (...) / else`` block (e.g. ``calv_flt_method`` branches on ``use_lsf``);
    each branch carries the guard conditions under which it applies.
    """
    enums: dict = {}
    for f in sorted(src_dir.glob("*.f90")):
        text = f.read_text(errors="replace")
        current_sub_group: str | None = None
        cond_stack: list = []   # active (param, value) guards, innermost last
        for _, line in _join_continuations(text):
            if _SUBROUTINE_RE.match(line):
                cond_stack = []   # conditions never cross subroutine boundaries
                mg = _SUBROUTINE_GROUP_RE.match(line)
                current_sub_group = mg.group(1) if mg else None
                continue

            if _ELSEIF_RE.match(line):
                if cond_stack:
                    cond_stack[-1] = (None, None)   # disjunction -> treat leniently
                continue
            if _ELSE_RE.match(line):
                if cond_stack:
                    cond_stack[-1] = _negate(cond_stack[-1])
                continue
            if _ENDIF_RE.match(line):
                if cond_stack:
                    cond_stack.pop()
                continue
            mif = _IF_THEN_RE.match(line)
            if mif:
                cond_stack.append(_parse_cond(mif.group(1)))
                continue

            for m in _ENUM_CALL_RE.finditer(line):
                group = _resolve_group(m.group("group"), current_sub_group)
                if not group:
                    continue
                allowed = [a.strip() for a in m.group("allowed").split("|") if a.strip()]
                conditions = [c for c in cond_stack if c != (None, None)]
                enums.setdefault((group, m.group("name")), []).append(
                    EnumConstraint(group=group, name=m.group("name"),
                                   allowed=allowed, conditions=conditions)
                )
    return enums


# --------------------------------------------------------------------------- #
# Enum snapshot (bundled fallback when no src/ is available)
# --------------------------------------------------------------------------- #
def enums_to_json(enums: dict) -> list:
    """Flatten ``{(group, name): [EnumConstraint]}`` to a JSON-serialisable list."""
    out = []
    for (g, n), ecs in sorted(enums.items()):
        for ec in ecs:
            out.append({"group": ec.group, "name": ec.name, "allowed": ec.allowed,
                        "conditions": [[p, v] for p, v in ec.conditions]})
    return out


def enums_from_json(data: list) -> dict:
    enums: dict = {}
    for d in data:
        conds = [(p, v) for p, v in d.get("conditions", [])]
        enums.setdefault((d["group"], d["name"]), []).append(
            EnumConstraint(group=d["group"], name=d["name"],
                           allowed=list(d["allowed"]), conditions=conds))
    return enums


def load_bundled_enums(path: Path) -> dict:
    if not path.is_file():
        return {}
    return enums_from_json(json.loads(path.read_text()))


# --------------------------------------------------------------------------- #
# Range / ordering from curated TOML
# --------------------------------------------------------------------------- #
def _load_toml(path: Path) -> tuple[list, list]:
    data = tomllib.loads(path.read_text())
    ranges = [
        RangeConstraint(
            group=r["group"], name=r["name"],
            min=r.get("min"), max=r.get("max"),
            min_inclusive=r.get("min_inclusive", True),
            max_inclusive=r.get("max_inclusive", True),
            note=r.get("note", ""),
        )
        for r in data.get("range", [])
    ]
    orders = [
        OrderConstraint(
            group=o["group"], left=o["left"], op=o["op"], right=o["right"],
            note=o.get("note", ""),
            when_name=o.get("when_name"), when_value=o.get("when_value"),
            when_contains=o.get("when_contains"),
        )
        for o in data.get("order", [])
    ]
    return ranges, orders


def load_constraints(enums: dict, toml_path: Path | None = None) -> Constraints:
    """Build Constraints from an already-resolved ``enums`` dict plus the curated
    range/ordering TOML. Enum resolution (live ``src/`` vs bundled snapshot) is
    decided by the caller (see ``schema.build_schema``)."""
    toml_path = toml_path or (DATA_DIR / "constraints.toml")
    ranges, orders = _load_toml(toml_path)
    return Constraints(enums=enums, ranges=ranges, orders=orders)
