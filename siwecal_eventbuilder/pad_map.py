"""
Pad → (x, y) position mapping for the SiW-ECAL prototype.

The converted ``ecal`` tree identifies each hit by ``(slab, chip, channel)`` but
stores no physical position. This module turns the ``(chip, channel)`` part into
a transverse position ``(x, y)`` in millimetres, reading a text map of the form::

    chip x0 y0 channel x y
    0 -67.05 -67.05 16 -86.3 -86.3
    ...

(one row per ``(chip, channel)``; 16 chips × 64 channels = 1024 rows). The map is
the same for every slab built from the same FEV board, but some slabs use a
different board/orientation (e.g. an FEV11 COB rotated wafer), so :class:`PadMap`
supports a per-slab override on top of a mandatory ``default`` map -- exactly the
dictionary mechanism validated in the standalone tests.

``chip`` here is the *geometric* chip index (0..15), which is how the map file is
keyed; :mod:`siwecal_eventbuilder.hit_collector` looks positions up with that same
geometric index (not the hardware ``chip_id`` written to ``hit_chip``).
"""

import math
import os


def load_mapping_file(path: str) -> dict:
    """Parse ``chip x0 y0 channel x y`` into ``{(chip, channel): (x, y)}``.

    Skips the header line, blank lines and ``#`` comments (mirrors the calibration
    file parsing in :mod:`siwecal_eventbuilder.calibration`).
    """
    mapping = {}
    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            items = stripped.split()
            if len(items) < 6 or not items[0].lstrip("-").isdigit():
                continue  # header row ("chip ...") or short/garbled line
            chip, channel = int(items[0]), int(items[3])
            x, y = float(items[4]), float(items[5])
            mapping[(chip, channel)] = (x, y)
    print(f"[PadMap] Loaded {len(mapping)} (chip,channel) entries from {path}")
    return mapping


class PadMap:
    """``(slab, chip, channel) → (x, y)`` with an optional per-slab override map.

    Construct via :meth:`from_files`. ``position`` returns ``(nan, nan)`` for a
    pad absent from the relevant map and keeps a running count of such misses so
    a summary warning can be emitted after a run.
    """

    def __init__(self, default_map: dict, slab_maps: dict = None):
        self._default = default_map
        self._slab_maps = slab_maps or {}
        self.n_unmapped = 0

    @classmethod
    def from_files(cls, files: dict, base_dir: str = "") -> "PadMap":
        """Build from a ``files`` dict mixing ``"default"`` with int slab keys.

        ``files`` must contain a ``"default"`` entry; any integer key is a per-slab
        override. Relative paths are resolved against ``base_dir``. A missing file
        raises :class:`FileNotFoundError`.
        """
        if "default" not in files:
            raise ValueError("PadMap files dict must contain a 'default' entry")

        def resolve(path):
            full = path if os.path.isabs(path) else os.path.join(base_dir, path)
            if not os.path.exists(full):
                raise FileNotFoundError(f"pad mapping file not found: {full}")
            return full

        default_map = load_mapping_file(resolve(files["default"]))
        slab_maps = {}
        for key, path in files.items():
            if key == "default":
                continue
            slab = int(key)
            slab_maps[slab] = load_mapping_file(resolve(path))
            print(f"[PadMap] slab {slab} uses override map {path}")
        return cls(default_map, slab_maps)

    def position(self, slab: int, chip: int, channel: int):
        """Return ``(x, y)`` in mm for this pad, or ``(nan, nan)`` if unmapped."""
        mapping = self._slab_maps.get(slab, self._default)
        xy = mapping.get((chip, channel))
        if xy is None:
            self.n_unmapped += 1
            return (math.nan, math.nan)
        return xy
