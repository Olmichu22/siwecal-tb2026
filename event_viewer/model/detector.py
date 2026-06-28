"""
Detector geometry for drawing: silicon pad planes and tungsten absorber blocks.

``DetectorModel`` turns the existing geometry pieces -- :class:`DetectorGeometry`
(layer z positions), :class:`PadMap` (chip/channel -> x,y) and the
``slab_z_positions.yml`` (z and W thickness per slab) -- into the numbers the
visualisers need:

* ``slab_z_array`` : per-hit z from the hit's slab index.
* ``pads_for_slab``: the (x, y) centres of every pad in a layer (grey grid).
* ``silicon_quads``/``tungsten_boxes``: corner coordinates for the 3-D Mesh3d
  representation of each Si plane and W plate.

The class owns no plotting code; it only produces geometry arrays.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

import numpy as np
import yaml

from .._geometry import DetectorGeometry
from .._pad_map import PadMap

# The W plate's downstream face sits this far upstream of its Si plane (mm);
# mirrors the convention documented in slab_z_positions.yml.
W_GAP_MM = 5.0


class DetectorModel:
    """Geometry provider for the silicon planes and tungsten absorbers."""

    def __init__(self, geometry: DetectorGeometry, pad_map: PadMap,
                 slab_z_mm: Tuple[float, ...], w_thickness_mm: Tuple[float, ...]):
        self.geometry = geometry
        self.pad_map = pad_map
        self.slab_z_mm = np.asarray(slab_z_mm, dtype=float)
        self.w_thickness_mm = np.asarray(w_thickness_mm, dtype=float)
        self._pads_cache: Dict[int, np.ndarray] = {}
        self._all_pads = self._collect_all_pads()
        self.pad_pitch = self._estimate_pitch(self._all_pads)
        self.x_extent, self.y_extent = self._plane_extent(self._all_pads)

    # --------------------------------------------------------- construction --
    @classmethod
    def from_config(cls, config) -> "DetectorModel":
        """Build from a :class:`event_viewer.config.ViewerConfig`."""
        with open(config.slab_z_yaml_path) as handle:
            doc = yaml.safe_load(handle) or {}
        slab_z = tuple(float(z) for z in doc.get("slab_z_mm", ()))
        w_thick = tuple(float(t) for t in doc.get("w_thickness_mm", ()))

        geometry = DetectorGeometry.from_mapping(
            {"slab_z_mm": slab_z} if slab_z else None)
        pad_map = PadMap.from_files(config.pad_map_files,
                                    base_dir=config.geometry_dir)
        if not slab_z:
            slab_z = geometry.slab_z_mm
        if not w_thick:
            w_thick = tuple(0.0 for _ in slab_z)
        return cls(geometry, pad_map, slab_z, w_thick)

    # ------------------------------------------------------------- z lookup --
    def slab_z(self, slab: int) -> float:
        return self.geometry.slab_z(slab)

    def slab_z_array(self, slab: np.ndarray) -> np.ndarray:
        """Vectorised per-hit z [mm] from an array of slab indices."""
        slab = np.asarray(slab, dtype=np.int64)
        out = np.full(slab.shape, np.nan, dtype=float)
        valid = (slab >= 0) & (slab < self.slab_z_mm.size)
        out[valid] = self.slab_z_mm[slab[valid]]
        return out

    # --------------------------------------------------------------- pads --
    def pads_for_slab(self, slab: int) -> np.ndarray:
        """``(N, 2)`` array of pad centres ``(x, y)`` [mm] for one layer."""
        if slab in self._pads_cache:
            return self._pads_cache[slab]
        points = []
        for chip in range(self.geometry.n_chips_per_slab):
            for channel in range(self.geometry.n_channels_per_chip):
                x, y = self.pad_map.position(slab, chip, channel)
                if np.isfinite(x) and np.isfinite(y):
                    points.append((x, y))
        pads = np.array(points, dtype=float) if points else np.empty((0, 2))
        self._pads_cache[slab] = pads
        return pads

    def _collect_all_pads(self) -> np.ndarray:
        """Union of pad centres across all layers (for pitch / extent)."""
        chunks = [self.pads_for_slab(s)
                  for s in range(self.geometry.n_slab_positions)]
        chunks = [c for c in chunks if c.size]
        return np.vstack(chunks) if chunks else np.empty((0, 2))

    @staticmethod
    def _estimate_pitch(pads: np.ndarray) -> float:
        """Pad pitch [mm] = smallest positive gap between unique x coordinates."""
        if pads.size == 0:
            return 5.5
        xs = np.unique(np.round(pads[:, 0], 3))
        diffs = np.diff(xs)
        diffs = diffs[diffs > 1e-6]
        return float(diffs.min()) if diffs.size else 5.5

    @staticmethod
    def _plane_extent(pads: np.ndarray) -> Tuple[Tuple[float, float],
                                                 Tuple[float, float]]:
        if pads.size == 0:
            return (-90.0, 90.0), (-90.0, 90.0)
        margin = 3.0
        x0, x1 = pads[:, 0].min() - margin, pads[:, 0].max() + margin
        y0, y1 = pads[:, 1].min() - margin, pads[:, 1].max() + margin
        return (float(x0), float(x1)), (float(y0), float(y1))

    # ------------------------------------------------------- 3-D primitives --
    def silicon_quads(self) -> List[Tuple[int, float]]:
        """``(slab, z)`` for each silicon plane to draw as a flat quad."""
        return [(s, float(self.slab_z_mm[s])) for s in range(self.slab_z_mm.size)]

    def tungsten_boxes(self) -> List[Tuple[int, float, float]]:
        """``(slab, z_downstream, z_upstream)`` of each W plate [mm].

        The plate's downstream face is ``W_GAP_MM`` upstream of its Si plane;
        the upstream face is a further ``w_thickness`` away. Plates of zero
        thickness are skipped.
        """
        boxes = []
        for s in range(self.slab_z_mm.size):
            thickness = float(self.w_thickness_mm[s]) if s < self.w_thickness_mm.size else 0.0
            if thickness <= 0:
                continue
            z_down = float(self.slab_z_mm[s]) - W_GAP_MM
            z_up = z_down - thickness
            boxes.append((s, z_down, z_up))
        return boxes
