"""
Configuration for the event viewer: where the data and geometry files live.

All filesystem locations default to the shared ``settings.yml`` resolved by
:mod:`siwecal_common.paths`, so the viewer follows the same data/geometry roots
as the event builder and the validation. Each can still be overridden from the
command line (see ``event_viewer/__main__.py``) without touching the rest.
"""

import os
from dataclasses import dataclass, field
from typing import Dict, List

from . import _paths as paths

from . import PROJECT_ROOT


def _default_pad_map_files() -> Dict:
    """``(chip,channel)->(x,y)`` maps: a mandatory default + per-slab overrides.

    Mirrors the mapping the event builder uses: FEV10 wafers everywhere, with the
    FEV11 chip-on-board (rotated) wafer overriding slab 12.
    """
    return {
        "default": "fev10_rotate_chip_channel_x_y_mapping.txt",
        12: "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt",
    }


def _default_data_dirs() -> List[str]:
    """Run/event search roots, plus the optional valcache dir, from settings.yml."""
    roots = list(paths.data_dirs())
    cache = paths.cache_dir()
    if cache and cache not in roots:
        roots.append(cache)
    return roots


@dataclass
class ViewerConfig:
    """Paths and display constants for one viewer session."""

    project_root: str = PROJECT_ROOT
    # Directories scanned (recursively) for the file dropdown. Defaults to the
    # settings.yml data roots (+ cache_dir); a single --data-dir overrides them.
    data_dirs: List[str] = field(default_factory=_default_data_dirs)
    geometry_dir: str = field(default_factory=paths.geometry_dir)
    tree_name: str = "ecal"

    # Geometry inputs.
    slab_z_yaml: str = field(
        default_factory=lambda: paths.geometry_file("slab_z_positions.yml"))
    pad_map_files: Dict = field(default_factory=_default_pad_map_files)
    n_layers: int = 15
    x0_mm: float = 3.5            # tungsten radiation length (for weighted energy)

    # Interactive MIP cut: when a file has no pre-computed MIP-cut branches the
    # slider falls back to an in-memory per-event recompute, which is too slow
    # for large runs. The slider is disabled above this many events (generate a
    # valcache with ``siwecal_validation`` to re-enable the cut on big files).
    max_recompute_events: int = 10000

    # Display.
    colorscale: str = "Viridis"
    host: str = "127.0.0.1"
    port: int = 8050

    @property
    def data_dir(self) -> str:
        """Primary data directory (kept for backwards compatibility)."""
        return self.data_dirs[0] if self.data_dirs else self.project_root

    @data_dir.setter
    def data_dir(self, value: str) -> None:
        """Setting a single data dir replaces the scanned roots with just it."""
        self.data_dirs = [value]

    def resolve(self, path: str) -> str:
        """Absolute path, interpreting relative paths against ``project_root``."""
        return path if os.path.isabs(path) else os.path.join(self.project_root, path)

    @property
    def slab_z_yaml_path(self) -> str:
        return self.resolve(self.slab_z_yaml)
