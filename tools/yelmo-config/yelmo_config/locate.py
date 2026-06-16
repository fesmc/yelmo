"""Locate the Yelmo repository artifacts the tool reads.

Resolution order for the defaults namelist:
  1. explicit ``--defaults`` / argument
  2. ``$YELMO_DEFAULTS`` env var
  3. ``$YELMO_ROOT/input/yelmo_defaults.nml`` env var
  4. walk up from the current working directory
  5. walk up from this installed package (covers editable installs)

The Fortran ``src/`` directory (used to extract enum constraints) is resolved
as the sibling ``src/`` of the located ``input/`` directory.
"""

from __future__ import annotations

import os
from pathlib import Path

DEFAULTS_RELPATH = Path("input") / "yelmo_defaults.nml"


def _walk_up(start: Path) -> Path | None:
    for d in [start, *start.parents]:
        cand = d / DEFAULTS_RELPATH
        if cand.is_file():
            return cand
    return None


def find_defaults(explicit: str | os.PathLike | None = None) -> Path:
    if explicit:
        p = Path(explicit).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"defaults file not found: {p}")
        return p.resolve()

    env_defaults = os.environ.get("YELMO_DEFAULTS")
    if env_defaults:
        p = Path(env_defaults).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"$YELMO_DEFAULTS set but not found: {p}")
        return p.resolve()

    env_root = os.environ.get("YELMO_ROOT")
    if env_root:
        p = Path(env_root).expanduser() / DEFAULTS_RELPATH
        if not p.is_file():
            raise FileNotFoundError(f"$YELMO_ROOT/{DEFAULTS_RELPATH} not found: {p}")
        return p.resolve()

    found = _walk_up(Path.cwd()) or _walk_up(Path(__file__).resolve().parent)
    if found is None:
        raise FileNotFoundError(
            "Could not locate input/yelmo_defaults.nml. Run from inside a Yelmo "
            "checkout, or set $YELMO_ROOT / $YELMO_DEFAULTS, or pass --defaults."
        )
    return found.resolve()


def repo_root(defaults: Path) -> Path:
    # defaults = <root>/input/yelmo_defaults.nml
    return defaults.parent.parent


def find_src(defaults: Path) -> Path | None:
    src = repo_root(defaults) / "src"
    return src if src.is_dir() else None
