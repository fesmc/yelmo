"""Resolve user-file group names to canonical schema groups.

A user par file may rename a component group, e.g. ``yelmo%nml_ydyn =
"ydyn_north"`` means the ``&ydyn_north`` group should be read/validated against
the canonical ``ydyn`` schema. This module builds that mapping from a user
namelist (falling back to the schema defaults for any key the user omits).

It also parses ``file[:group]`` selector tokens used by the ``diff`` command.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from . import namelist as nl
from .schema import COMPONENT_KEYS, Schema


def resolve_groups(user: nl.Namelist, schema: Schema) -> dict:
    """Return ``{actual_group_name_in_file: canonical_group}``.

    Always includes ``yelmo`` -> ``yelmo``. For every component key, the active
    group name is the user's value if present, else the schema default.
    """
    mapping = {"yelmo": "yelmo"}
    defaults_cmap = schema.component_map()           # nml_* -> canonical
    user_yelmo = user.get("yelmo")

    for key in COMPONENT_KEYS:
        canonical = defaults_cmap.get(key)
        if canonical is None:
            continue
        actual = canonical
        if user_yelmo is not None:
            p = user_yelmo.get(key)
            if p is not None:
                actual = nl.normalize(p.raw)
        mapping[str(actual)] = str(canonical)
    return mapping


def canonical_of(group_name: str, user: nl.Namelist, schema: Schema) -> str | None:
    """Canonical schema group for an actual group name in ``user`` (or None)."""
    return resolve_groups(user, schema).get(group_name)


# --------------------------------------------------------------------------- #
# Selector parsing:  PATH  or  PATH:GROUP
# --------------------------------------------------------------------------- #
@dataclass
class Selector:
    path: Path
    group: str | None = None  # actual group name within the file, if specified


def parse_selector(token: str) -> Selector:
    """Parse ``file`` or ``file:group``.

    A ``:`` is only treated as a group separator when what follows looks like a
    group name (no path separator) and the part before it is an existing path or
    has a namelist-ish extension — so Windows-style ``C:\\...`` and paths with
    colons are not mis-split. In practice Yelmo files have no colons.
    """
    if ":" in token:
        head, _, tail = token.rpartition(":")
        if head and tail and "/" not in tail and "\\" not in tail:
            return Selector(path=Path(head), group=tail)
    return Selector(path=Path(token))
