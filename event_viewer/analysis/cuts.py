"""
Dynamic per-event selection cuts.

A :class:`CutModel` is an ordered set of ``variable in [lo, hi]`` ranges. It maps
a per-event DataFrame to a boolean mask / index list, and round-trips to a plain
list of dicts so it can live inside a Dash ``dcc.Store``.

Cutting on a variable naturally drops events whose value is NaN (e.g. shower
variables on non-showers), because ``NaN`` fails both comparisons -- which is the
desired behaviour when the user explicitly selects a range on that variable.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Cut:
    """A single inclusive range cut on one event-level variable."""

    variable: str
    lo: float
    hi: float

    def mask(self, df: pd.DataFrame) -> np.ndarray:
        if self.variable not in df.columns:
            return np.ones(len(df), dtype=bool)
        col = df[self.variable].to_numpy(dtype=float)
        return (col >= self.lo) & (col <= self.hi)


class CutModel:
    """An AND-combination of :class:`Cut` ranges."""

    def __init__(self, cuts: List[Cut] = None):
        self.cuts: List[Cut] = list(cuts) if cuts else []

    @property
    def is_empty(self) -> bool:
        return not self.cuts

    # ----------------------------------------------------------- selection --
    def mask(self, df: pd.DataFrame) -> np.ndarray:
        """Boolean mask of events passing every cut."""
        keep = np.ones(len(df), dtype=bool)
        for cut in self.cuts:
            keep &= cut.mask(df)
        return keep

    def passing_indices(self, df: pd.DataFrame) -> np.ndarray:
        """Integer positions (tree entries) of events passing every cut."""
        return np.flatnonzero(self.mask(df))

    # -------------------------------------------------------- (de)serialise --
    def to_store(self) -> List[dict]:
        return [{"variable": c.variable, "lo": c.lo, "hi": c.hi} for c in self.cuts]

    @classmethod
    def from_store(cls, data) -> "CutModel":
        if not data:
            return cls()
        return cls([Cut(d["variable"], float(d["lo"]), float(d["hi"]))
                    for d in data])
