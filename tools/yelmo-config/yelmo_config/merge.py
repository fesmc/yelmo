"""Merge user overrides onto the schema defaults to get the *effective* config.

This mirrors what the Yelmo library does at load time: read every parameter
from ``yelmo_defaults.nml``, then overwrite any the user par file specifies.
Group renames (``nml_ydyn = "ydyn_north"``) are respected.
"""

from __future__ import annotations

import copy

from . import namelist as nl
from .resolve import resolve_groups
from .schema import Schema


def effective_group(user: nl.Namelist, actual_group: str, schema: Schema,
                    include_defaults: bool = True) -> dict:
    """Effective ``{name: raw}`` for one *actual* group in ``user``.

    The defaults for the group's canonical schema are overlaid with whatever the
    user provides in ``&actual_group``. Unknown user params are included too.
    With ``include_defaults=False`` only the params literally present in the
    user file's group are returned (used by ``diff --raw``).
    """
    mapping = resolve_groups(user, schema)
    canonical = mapping.get(actual_group, actual_group)

    out: dict = {}
    default_grp = schema.nml.get(canonical)
    if include_defaults and default_grp is not None:
        for p in default_grp.items:
            if isinstance(p, nl.Param):
                out[p.name] = p.raw

    user_grp = user.get(actual_group)
    if user_grp is not None:
        for p in user_grp.items:
            if isinstance(p, nl.Param):
                out[p.name] = p.raw
    return out


def effective_all(user: nl.Namelist, schema: Schema) -> dict:
    """Effective ``{canonical_group: {name: raw}}`` across all library groups."""
    mapping = resolve_groups(user, schema)
    actual_by_canonical = {canon: actual for actual, canon in mapping.items()}

    out: dict = {}
    for canonical in schema.groups:
        actual = actual_by_canonical.get(canonical, canonical)
        out[canonical] = effective_group(user, actual, schema)
    return out


def complete_namelist(user: nl.Namelist | None, schema: Schema) -> nl.Namelist:
    """A full namelist = defaults (comments/sections intact) + user overrides.

    Group names stay canonical (matching the documented defaults structure).
    """
    out = copy.deepcopy(schema.nml)
    if user is None:
        return out

    mapping = resolve_groups(user, schema)
    actual_by_canonical = {canon: actual for actual, canon in mapping.items()}

    for canonical, group in out.groups.items():
        actual = actual_by_canonical.get(canonical, canonical)
        user_grp = user.get(actual)
        if user_grp is None:
            continue
        for p in user_grp.items:
            if isinstance(p, nl.Param):
                group.set_value(p.name, p.raw)
    return out
