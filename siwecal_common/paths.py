"""
Single source of truth for filesystem locations across the monorepo.

Every path the three packages need -- where the run/event ROOT files live, the
calibration files, the geometry maps, the ``data_reference`` configs, the output
directory, and the optional valcache directory -- is resolved here from a single
``settings.yml`` at the repo root. No analysis module hard-codes an absolute
path; they all ask this module instead.

Resolution rules
----------------
* The repo root is the directory that contains ``settings.yml`` (by default the
  parent of this package). Override the settings file with the ``SIWECAL_SETTINGS``
  environment variable.
* Relative values in ``settings.yml`` (e.g. ``./geometry``) resolve against the
  repo root. Absolute values are used verbatim.
* ``data_dir`` may be a single path or a list of search roots; ``data_dirs()``
  always returns a list, and :func:`resolve_input` walks them in order.
* A missing or empty ``settings.yml`` is not fatal: built-in defaults relative to
  the repo root are used, so the repo works out-of-the-box.
"""

import os
from functools import lru_cache
from typing import List, Optional

import yaml

# ------------------------------------------------------------------ roots ---
PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(PACKAGE_DIR)

# Built-in defaults (used when a key is absent from settings.yml). Relative to
# the repo root; mirror the keys documented in settings.example.yml.
_DEFAULTS = {
    "data_dir": ["."],
    "calib_dir": "./calibration",
    "geometry_dir": "./geometry",
    "configs_dir": "./configs/data",
    "output_dir": "./validation_output",
    "cache_dir": None,
    "pid_dir": None,
}


def settings_file() -> str:
    """Path to the active ``settings.yml`` (``$SIWECAL_SETTINGS`` or repo root)."""
    return os.environ.get("SIWECAL_SETTINGS",
                          os.path.join(REPO_ROOT, "settings.yml"))


@lru_cache(maxsize=1)
def _load() -> dict:
    """Parse ``settings.yml`` once; absent/empty file -> built-in defaults only."""
    path = settings_file()
    document = {}
    if os.path.exists(path):
        with open(path) as handle:
            document = yaml.safe_load(handle) or {}
    merged = dict(_DEFAULTS)
    merged.update({k: v for k, v in document.items() if v is not None
                   or k in ("cache_dir", "pid_dir")})
    return merged


def reload() -> None:
    """Drop the cached settings (call after editing ``settings.yml`` at runtime)."""
    _load.cache_clear()


# ------------------------------------------------------------- resolution ---
def resolve(path: str) -> str:
    """Absolute path: relative values resolve against the repo root."""
    return path if os.path.isabs(path) else os.path.normpath(
        os.path.join(REPO_ROOT, path))


def _get(key: str):
    return _load().get(key, _DEFAULTS.get(key))


# ------------------------------------------------------------- public API ---
def data_dirs() -> List[str]:
    """All run/event search roots, as a list of absolute paths (>=1)."""
    value = _get("data_dir")
    items = value if isinstance(value, (list, tuple)) else [value]
    return [resolve(str(item)) for item in items if item]


def data_dir() -> str:
    """The primary (first) data search root."""
    roots = data_dirs()
    return roots[0] if roots else REPO_ROOT


def resolve_input(rel_or_abs: str) -> str:
    """Locate an input file/run given a path that may be relative to a data root.

    Absolute paths are returned unchanged. Relative paths are looked up under
    each :func:`data_dirs` root, returning the first that exists; if none exists
    the path under the first root is returned (so callers get a sensible,
    reportable "not found" path).
    """
    if os.path.isabs(rel_or_abs):
        return rel_or_abs
    roots = data_dirs()
    for root in roots:
        candidate = os.path.join(root, rel_or_abs)
        if os.path.exists(candidate):
            return candidate
    return os.path.join(roots[0], rel_or_abs) if roots else rel_or_abs


def calib_dir() -> str:
    """Directory holding the calibration (pedestal / MIP) files."""
    return resolve(str(_get("calib_dir")))


def geometry_dir() -> str:
    """Directory holding the pad maps and slab-geometry YAMLs."""
    return resolve(str(_get("geometry_dir")))


def geometry_file(name: str) -> str:
    """Absolute path of a file inside :func:`geometry_dir`."""
    return os.path.join(geometry_dir(), name)


def configs_dir() -> str:
    """Directory holding the ``data_reference*.yml`` run-list configs."""
    return resolve(str(_get("configs_dir")))


def config_file(name: str) -> str:
    """Absolute path of a file inside :func:`configs_dir`."""
    return os.path.join(configs_dir(), name)


def output_dir() -> str:
    """Default output base for validation plots/results."""
    return resolve(str(_get("output_dir")))


def cache_dir() -> Optional[str]:
    """Directory for ``*.valcache.root`` files, or ``None`` to keep them next to
    the input (the default behaviour)."""
    value = _get("cache_dir")
    return resolve(str(value)) if value else None


def pid_dir() -> Optional[str]:
    """Directory for the Gaudi ``ecal_pid_*.root`` (EDM4hep) outputs, or ``None``
    to keep them next to each input ``ecal_<run>.root`` (the default)."""
    value = _get("pid_dir")
    return resolve(str(value)) if value else None
