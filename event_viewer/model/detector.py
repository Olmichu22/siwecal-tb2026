"""
Detector geometry for drawing: silicon pad planes and tungsten absorber blocks.

``DetectorModel`` turns the existing geometry pieces -- :class:`DetectorGeometry`
(layer z positions), :class:`PadMap` (chip/channel -> x,y) and the
``slab_z_positions.yml`` (z and W thickness per slab) -- into the numbers the
visualisers need:

* ``slab_z_array`` : per-hit z from the hit's slab index.
* ``pads_for_slab``: the (x, y) centres of every pad in a layer (grey grid).
* ``wafer_rects``  : the active area of each silicon sensor in a layer.
* ``silicon_quads``/``tungsten_boxes``: corner coordinates for the 3-D Mesh3d
  representation of each Si plane and W plate.

A layer is not a continuous silicon plane: it is an array of silicon sensors,
each carrying 16x16 pads with an inactive 0.61 mm rim around them, so the
boundary between two butted sensors is a 1.22 mm dead cross. ``wafer_rects``
recovers that segmentation from the pad map itself -- no wafer count or gap
width is hardcoded here -- so the drawn geometry follows whatever mapping file
the run uses. The rectangle it returns is the sensor's active area, and what
falls between two of them is the dead region, where no hit can appear.

This is the test-beam viewer, reading the real pad maps
(mappings/fev10*/fev11*). Those tables are written on the measured spec-sheet
grid (see mappings/regenerate_pad_maps.py): pads 5.53 mm apart inside a sensor
and 6.75 mm across a boundary, the boundary being the two 0.61 mm sensor rims.
The inference below uses only the ratio of those two steps, so it is
independent of the exact numbers and also handled the older nominal map.

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

# A step between neighbouring pad coordinates wider than this many pitches is
# read as the boundary between two sensors rather than the pad pitch. The map
# has 5.53 mm inside a sensor and 6.75 mm across a boundary -- 1.22 pitches --
# because the dead region is only the two 0.61 mm rims. Anything between ~1.05
# and ~1.2 separates the two; the older nominal map, with its 7.6 mm boundary,
# also passes comfortably.
WAFER_GAP_PITCHES = 1.1

# Active area of one pad [mm], inside the pixel pitch: 5.52 of diode plus a
# 0.01 inter-pad margin gives the 5.53 pitch (spec-sheet numbers). Only the 3-D
# view uses the active size, to draw the square of a hit; the 2-D view fills
# whole pitch cells so that they stay contiguous. The pitch itself is measured
# off the pad map and is only a fallback here.
PAD_ACTIVE_MM = 5.52
PAD_PITCH_MM = 5.53


class DetectorModel:
    """Geometry provider for the silicon planes and tungsten absorbers."""

    def __init__(self, geometry: DetectorGeometry, pad_map: PadMap,
                 slab_z_mm: Tuple[float, ...], w_thickness_mm: Tuple[float, ...]):
        self.geometry = geometry
        self.pad_map = pad_map
        self.slab_z_mm = np.asarray(slab_z_mm, dtype=float)
        self.w_thickness_mm = np.asarray(w_thickness_mm, dtype=float)
        self._pads_cache: Dict[int, np.ndarray] = {}
        self._grid_cache: Dict[int, Tuple] = {}
        self._all_pads = self._collect_all_pads()
        self.pad_pitch = self._estimate_pitch(self._all_pads)
        # A pad can never be drawn wider than the pitch of the map it sits on,
        # or neighbouring pads would overlap.
        self.pad_active = min(PAD_ACTIVE_MM, self.pad_pitch)
        self.x_extent, self.y_extent = self._plane_extent(self._all_pads)
        self.wafer_rects = self._compute_wafer_rects(self._all_pads,
                                                     self.pad_pitch)

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
            return PAD_PITCH_MM
        xs = np.unique(np.round(pads[:, 0], 3))
        diffs = np.diff(xs)
        diffs = diffs[diffs > 1e-6]
        return float(diffs.min()) if diffs.size else PAD_PITCH_MM

    # ------------------------------------------------------------- wafers --
    @staticmethod
    def _bands(values: np.ndarray, pitch: float) -> List[Tuple[float, float]]:
        """Group pad coordinates into wafers along one axis.

        Pads are evenly spaced inside a sensor, so every step wider than
        ``WAFER_GAP_PITCHES * pitch`` marks the dead region where one sensor
        ends and the next begins. Returns the ``(first, last)`` pad coordinate
        of each sensor, in increasing order.
        """
        uniq = np.unique(np.round(values, 3))
        if uniq.size == 0:
            return []
        cuts = np.flatnonzero(np.diff(uniq) > WAFER_GAP_PITCHES * pitch) + 1
        return [(float(band[0]), float(band[-1]))
                for band in np.split(uniq, cuts) if band.size]

    def _compute_wafer_rects(self, pads: np.ndarray, pitch: float
                             ) -> List[Tuple[float, float, float, float]]:
        """``(x0, x1, y0, y1)`` of every silicon sensor, in mm.

        The pads butt against each other inside a sensor, so the active area
        reaches exactly half a pitch beyond the outermost pad *centre*. What is
        left between two such rectangles is the pair of inactive rims, i.e. the
        real dead cross; drawing the sensors' outer edges instead would show
        nothing, since butted sensors touch. A map with no gap at all (a single
        sensor, or a hypothetical continuous plane) degrades to one rectangle
        covering the whole active area, which is the correct picture.
        """
        if pads.size == 0:
            return [(self.x_extent[0], self.x_extent[1],
                     self.y_extent[0], self.y_extent[1])]
        half = pitch / 2.0
        x_bands = self._bands(pads[:, 0], pitch)
        y_bands = self._bands(pads[:, 1], pitch)
        return [(x0 - half, x1 + half, y0 - half, y1 + half)
                for (y0, y1) in y_bands for (x0, x1) in x_bands]

    # --------------------------------------------------------- cell grid --
    @staticmethod
    def _axis_cells(centres: np.ndarray, pitch: float
                    ) -> Tuple[np.ndarray, Dict[float, int]]:
        """Cell edges along one axis, plus ``rounded centre -> column``.

        Every pad becomes a cell exactly one pitch wide, and each boundary
        between two sensors becomes one extra cell of its own, never filled.
        Without that
        spacer the edges would not be contiguous, and a renderer asked for
        contiguous edges would silently stretch the pads either side of the gap
        to cover it, painting over the dead region.
        """
        half = pitch / 2.0
        edges = [float(centres[0]) - half]
        index: Dict[float, int] = {}
        for n, centre in enumerate(centres):
            centre = float(centre)
            if n and centre - float(centres[n - 1]) > WAFER_GAP_PITCHES * pitch:
                edges.append(centre - half)   # closes the guard-ring cell
            index[round(centre, 3)] = len(edges) - 1
            edges.append(centre + half)
        return np.asarray(edges, dtype=float), index

    def pad_cell_grid(self, slab: int) -> Tuple:
        """``(x_edges, y_edges, x_index, y_index)`` for one layer.

        Drives the 2-D view: the edges give every pad its true pitch-wide
        footprint at any zoom level (a marker's size is fixed in pixels and
        cannot), and the indices place a hit at ``(x, y)`` on the right cell.
        """
        cached = self._grid_cache.get(slab)
        if cached is not None:
            return cached
        pads = self.pads_for_slab(slab)
        if pads.size == 0:
            grid = (np.empty(0), np.empty(0), {}, {})
        else:
            x_edges, x_index = self._axis_cells(
                np.unique(np.round(pads[:, 0], 3)), self.pad_pitch)
            y_edges, y_index = self._axis_cells(
                np.unique(np.round(pads[:, 1], 3)), self.pad_pitch)
            grid = (x_edges, y_edges, x_index, y_index)
        self._grid_cache[slab] = grid
        return grid

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
