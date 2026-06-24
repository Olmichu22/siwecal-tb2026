"""
Optional YAML configuration loader for the SiW-ECAL event builder.

By design the source of truth for every tunable is the dataclass default in
:mod:`siwecal_eventbuilder.config` / :mod:`siwecal_eventbuilder.geometry`. This
module adds an *optional* layer on top: a ``config.yml`` file whose presence is
not required and which only needs to list the keys the user wants to change.

Precedence (lowest to highest), with the CLI layer applied later in
:mod:`siwecal_eventbuilder.cli`::

    dataclass defaults  ->  config.yml  ->  command-line flags

The recognised top-level sections of ``config.yml`` are:

* ``builder``     -- overrides for :class:`~siwecal_eventbuilder.config.BuilderConfig`.
* ``geometry``    -- overrides for :class:`~siwecal_eventbuilder.geometry.DetectorGeometry`.
* ``paths``       -- optional ``main_path`` / ``calibration_dir`` overrides.
* ``calibration`` -- optional calibration source (mode + files + mip_run).
* ``mapping``     -- optional pad → (x, y) position map (default + per-slab files).

Anything else at the top level is reported as a warning so typos in section
names do not pass silently.
"""

import os
import sys
from dataclasses import dataclass, field
from typing import Optional

import yaml

from .config import BuilderConfig
from .geometry import DetectorGeometry

KNOWN_SECTIONS = frozenset({"builder", "geometry", "paths", "calibration",
                            "mapping"})


def load_config_file(path: Optional[str]) -> Optional[dict]:
    """Parse a YAML config file, returning ``None`` if it is absent.

    A missing file is not an error: the builder simply runs on dataclass
    defaults. A present-but-empty file parses to ``None`` and is treated the
    same way.
    """
    if not path or not os.path.exists(path):
        return None
    with open(path) as handle:
        return yaml.safe_load(handle)


@dataclass(frozen=True)
class AppSettings:
    """Resolved configuration bundle (defaults merged with ``config.yml``).

    This is the ``defaults -> config.yml`` layer only. The CLI applies its own
    overrides on top of these values afterwards.
    """

    config: BuilderConfig = field(default_factory=BuilderConfig)
    geometry: DetectorGeometry = field(default_factory=DetectorGeometry)
    paths: dict = field(default_factory=dict)
    calibration: dict = field(default_factory=dict)
    mapping: dict = field(default_factory=dict)

    @classmethod
    def from_yaml(cls, path: Optional[str]) -> "AppSettings":
        """Build settings from an optional ``config.yml``.

        Parameters
        ----------
        path:
            Path to the YAML file, or ``None``. If the file does not exist the
            returned settings are exactly the dataclass defaults.

        Raises
        ------
        ValueError
            Propagated from :meth:`BuilderConfig.from_mapping` /
            :meth:`DetectorGeometry.from_mapping` when an option name is invalid.
        """
        document = load_config_file(path)
        if not document:
            return cls()

        unknown_sections = [key for key in document if key not in KNOWN_SECTIONS]
        if unknown_sections:
            print(
                "WARNING: ignoring unknown section(s) in config file: "
                f"{', '.join(sorted(unknown_sections))}. "
                f"Known sections: {', '.join(sorted(KNOWN_SECTIONS))}.",
                file=sys.stderr,
            )

        return cls(
            config=BuilderConfig.from_mapping(document.get("builder")),
            geometry=DetectorGeometry.from_mapping(document.get("geometry")),
            paths=dict(document.get("paths") or {}),
            calibration=dict(document.get("calibration") or {}),
            mapping=dict(document.get("mapping") or {}),
        )
