"""
A single reconstructed event: its per-hit arrays plus its per-event metrics.

This is a plain data container produced by :class:`EventDataset`; the visualisers
consume it directly. ``z`` is filled from the hit's slab via the detector model,
so the event already carries physical 3-D positions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict

import numpy as np


@dataclass
class Event:
    """Per-hit arrays (all aligned, length = number of hit channels) + metrics."""

    index: int                       # tree entry index
    x: np.ndarray                    # hit x [mm]
    y: np.ndarray                    # hit y [mm]
    z: np.ndarray                    # hit z [mm] (from slab)
    slab: np.ndarray                 # hit layer index
    chip: np.ndarray
    chan: np.ndarray
    energy: np.ndarray               # hit energy [MIP]
    hg: np.ndarray                   # high-gain ADC
    lg: np.ndarray                   # low-gain ADC
    metrics: Dict[str, float] = field(default_factory=dict)

    @property
    def n_hits(self) -> int:
        return int(self.x.size)
