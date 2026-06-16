"""Locate the Yelmo artifacts the tool reads, with a bundled fallback.

The tool needs ``input/yelmo_defaults.nml`` (params/defaults/docs) and the enum
constraints extracted from ``src/*.f90``. When run inside a Yelmo checkout it
uses those live files; otherwise it falls back to snapshots bundled in the
package (``data/yelmo_defaults.nml`` / ``data/enums.json``), so ``yelmo-config``
works without a local copy of Yelmo.

Resolution order for the defaults namelist:
  1. explicit ``--defaults`` / argument            -> source "explicit"
  2. ``$YELMO_DEFAULTS`` / ``$YELMO_ROOT`` env var  -> source "env"
  3. walk up from the current working directory     -> source "local"
  4. walk up from this installed package            -> source "local"
  5. the bundled snapshot                           -> source "bundled"
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

DEFAULTS_RELPATH = Path("input") / "yelmo_defaults.nml"

DATA_DIR = Path(__file__).resolve().parent / "data"
BUNDLED_DEFAULTS = DATA_DIR / "yelmo_defaults.nml"
BUNDLED_ENUMS = DATA_DIR / "enums.json"


@dataclass
class DefaultsResolution:
    path: Path           # the defaults file actually used
    source: str          # "explicit" | "env" | "local" | "bundled"
    repo_root: Path | None   # checkout root (for src/ and input paths), or None
    src: Path | None         # src/ dir, if present


def _walk_up(start: Path) -> Path | None:
    for d in [start, *start.parents]:
        cand = d / DEFAULTS_RELPATH
        if cand.is_file():
            return cand
    return None


def _find_local_defaults(explicit: str | os.PathLike | None = None):
    """Return ``(path, source)`` for a real checkout, or ``None`` if not found.

    Never returns the bundled snapshot; raises only when an explicit/env path is
    given but does not exist.
    """
    if explicit:
        p = Path(explicit).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"defaults file not found: {p}")
        return p.resolve(), "explicit"

    env_defaults = os.environ.get("YELMO_DEFAULTS")
    if env_defaults:
        p = Path(env_defaults).expanduser()
        if not p.is_file():
            raise FileNotFoundError(f"$YELMO_DEFAULTS set but not found: {p}")
        return p.resolve(), "env"

    env_root = os.environ.get("YELMO_ROOT")
    if env_root:
        p = Path(env_root).expanduser() / DEFAULTS_RELPATH
        if not p.is_file():
            raise FileNotFoundError(f"$YELMO_ROOT/{DEFAULTS_RELPATH} not found: {p}")
        return p.resolve(), "env"

    found = _walk_up(Path.cwd()) or _walk_up(Path(__file__).resolve().parent)
    if found is not None:
        return found.resolve(), "local"
    return None


def find_local_defaults(explicit=None) -> Path:
    """Like ``_find_local_defaults`` but require a checkout (no bundled fallback)."""
    res = _find_local_defaults(explicit)
    if res is None:
        raise FileNotFoundError(
            "Could not locate input/yelmo_defaults.nml. Run from inside a Yelmo "
            "checkout, or set $YELMO_ROOT / $YELMO_DEFAULTS, or pass --defaults."
        )
    return res[0]


def resolve_defaults(explicit: str | os.PathLike | None = None) -> DefaultsResolution:
    """Resolve the defaults file, falling back to the bundled snapshot."""
    res = _find_local_defaults(explicit)
    if res is None:
        return DefaultsResolution(BUNDLED_DEFAULTS, "bundled", repo_root=None, src=None)
    path, source = res
    root = path.parent.parent
    src = root / "src"
    return DefaultsResolution(path, source, repo_root=root,
                              src=src if src.is_dir() else None)


def find_src(defaults: Path) -> Path | None:
    src = defaults.parent.parent / "src"
    return src if src.is_dir() else None
