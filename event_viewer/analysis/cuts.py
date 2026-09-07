"""
Dynamic per-event selection cuts.

A :class:`CutModel` is an ordered set of ``variable in [lo, hi]`` ranges. It maps
a per-event DataFrame to a boolean mask / index list, and round-trips to a plain
list of dicts so it can live inside a Dash ``dcc.Store``.

Cutting on a variable naturally drops events whose value is NaN (e.g. shower
variables on non-showers), because ``NaN`` fails both comparisons -- which is the
desired behaviour when the user explicitly selects a range on that variable.

Complementary selections
------------------------
A cut can be flagged ``invert``, and the flagged ones are complemented **as a
group**: the selection is every event that passes all the plain cuts *and fails
the combination of the inverted ones*.

    keep = AND(plain cuts) AND NOT( AND(inverted cuts) )

Complementing as a group rather than one cut at a time is the distinction that
makes the feature useful, because ``NOT(A AND B)`` is not ``NOT A AND NOT B``.
With a window on ``n_hits`` left plain and a cut on ``dl_score`` inverted, the
selection is "the events inside the occupancy window that the score threw away"
-- the question you actually ask of a discriminant. Inverting both instead gives
everything either cut rejected.

NaN is excluded from the complement too. An event with no value for an inverted
variable did not *fail* that cut, it has nothing to compare, and letting it in
would fill the complement with events that were never candidates.

Quantile cuts
-------------
A cut can instead be expressed in ``mode="quantile"``, where ``lo``/``hi`` are
fractions in [0, 1] and the actual edges are the corresponding percentiles of the
variable. Absolute thresholds are not comparable across runs -- the same
``dl_score < 0.52`` keeps 0.5% of one beam energy and 40% of another -- while
"the lowest 10%" means the same thing everywhere.

The percentiles are resolved at ``mask()`` time rather than being frozen into a
value when the cut is created, because the reference sample depends on the MIP
threshold: moving that slider would otherwise leave a stale number silently
claiming to be a percentile it no longer is.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Cut:
    """A single inclusive range cut on one event-level variable.

    ``invert`` does not negate this cut on its own; it marks it as belonging to
    the group that :class:`CutModel` complements together.

    With ``mode="quantile"`` the range is given as fractions in [0, 1] of the
    variable's own distribution instead of absolute values: ``lo=0.1, hi=0.9``
    keeps the events between the 10th and the 90th percentile.
    """

    variable: str
    lo: float
    hi: float
    invert: bool = False
    mode: str = "value"

    def resolve_range(self, df: pd.DataFrame):
        """The ``(lo, hi)`` edges of this cut in the units of the variable.

        For a quantile cut the edges are the ``lo``/``hi`` percentiles of the
        finite values of the column, so the caller (e.g. the histogram shading)
        sees the same numbers the mask compares against. Returns ``None`` when
        the column has no finite value to take percentiles of.

        ``df`` is expected to be the *whole* table the cuts are applied to, not
        a subset already filtered by the other cuts: percentiles taken over a
        pre-filtered subset would depend on the order the cuts were applied and
        become circular once combined with the complement. Every caller today
        passes the full (MIP-threshold filtered) table.
        """
        if self.mode != "quantile":
            return float(self.lo), float(self.hi)
        if self.variable not in df.columns:
            return None
        col = df[self.variable].to_numpy(dtype=float)
        col = col[np.isfinite(col)]
        if col.size == 0:
            return None
        qlo, qhi = sorted((min(max(float(self.lo), 0.0), 1.0),
                           min(max(float(self.hi), 0.0), 1.0)))
        return float(np.quantile(col, qlo)), float(np.quantile(col, qhi))

    def mask(self, df: pd.DataFrame) -> np.ndarray:
        if self.variable not in df.columns:
            return np.ones(len(df), dtype=bool)
        edges = self.resolve_range(df)
        if edges is None:
            return np.zeros(len(df), dtype=bool)
        lo, hi = edges
        col = df[self.variable].to_numpy(dtype=float)
        return (col >= lo) & (col <= hi)

    def defined(self, df: pd.DataFrame) -> np.ndarray:
        """Events that have a value for this variable at all."""
        if self.variable not in df.columns:
            return np.ones(len(df), dtype=bool)
        return np.isfinite(df[self.variable].to_numpy(dtype=float))


class CutModel:
    """An AND-combination of :class:`Cut` ranges, with an inverted subgroup."""

    def __init__(self, cuts: List[Cut] = None):
        self.cuts: List[Cut] = list(cuts) if cuts else []

    @property
    def is_empty(self) -> bool:
        return not self.cuts

    @property
    def inverted(self) -> List[Cut]:
        return [c for c in self.cuts if c.invert]

    # ----------------------------------------------------------- selection --
    def mask(self, df: pd.DataFrame) -> np.ndarray:
        """Boolean mask: passes every plain cut, fails the inverted group."""
        keep = np.ones(len(df), dtype=bool)
        for cut in self.cuts:
            if not cut.invert:
                keep &= cut.mask(df)

        inverted = self.inverted
        if inverted:
            group = np.ones(len(df), dtype=bool)
            defined = np.ones(len(df), dtype=bool)
            for cut in inverted:
                group &= cut.mask(df)
                defined &= cut.defined(df)
            keep &= defined & ~group
        return keep

    def passing_indices(self, df: pd.DataFrame) -> np.ndarray:
        """Integer positions (tree entries) of events passing the selection."""
        return np.flatnonzero(self.mask(df))

    # -------------------------------------------------------- (de)serialise --
    def to_store(self) -> List[dict]:
        return [{"variable": c.variable, "lo": c.lo, "hi": c.hi,
                 "invert": bool(c.invert), "mode": c.mode} for c in self.cuts]

    @classmethod
    def from_store(cls, data) -> "CutModel":
        if not data:
            return cls()
        # ``invert`` and ``mode`` are read with defaults so a store written by an
        # older session (or a bookmarked layout) still loads.
        return cls([Cut(d["variable"], float(d["lo"]), float(d["hi"]),
                        bool(d.get("invert", False)),
                        str(d.get("mode", "value")))
                    for d in data])
