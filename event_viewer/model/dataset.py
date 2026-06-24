"""
A per-file dataset: the event table plus on-demand single-event access.

``EventDataset`` ties together an :class:`EventFileReader` (raw I/O) and a
:class:`DetectorModel` (geometry), exposing exactly what the controller and
visualisers need: the per-event DataFrame, the indices passing a cut, and a fully
built :class:`Event` (per-hit arrays with physical z + the row of metrics).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .event import Event


class EventDataset:
    """Per-file view combining the reader's table with single-event access."""

    def __init__(self, reader, detector):
        self.reader = reader
        self.detector = detector
        self._table = reader.event_table(0.0)

    # -------------------------------------------------------------- table --
    @property
    def table(self) -> pd.DataFrame:
        return self._table

    @property
    def n_events(self) -> int:
        return self.reader.n_events

    @property
    def has_metrics(self) -> bool:
        return self.reader.has_metrics

    def feature_columns(self):
        return self.reader.feature_columns()

    def passing(self, cut_model) -> np.ndarray:
        """Tree-entry indices passing ``cut_model`` (all events if empty)."""
        if cut_model is None or cut_model.is_empty:
            return np.arange(self.n_events)
        return cut_model.passing_indices(self._table)

    # ------------------------------------------------------------- event --
    def get_event(self, index: int) -> Event:
        """Build the :class:`Event` at tree entry ``index``."""
        hits = self.reader.read_hits(index)
        slab = np.asarray(hits.get("hit_slab", []), dtype=np.int64)
        row = self._table.iloc[index].to_dict() if index < len(self._table) else {}
        return Event(
            index=index,
            x=np.asarray(hits.get("hit_x", []), dtype=float),
            y=np.asarray(hits.get("hit_y", []), dtype=float),
            z=self.detector.slab_z_array(slab),
            slab=slab,
            chip=np.asarray(hits.get("hit_chip", []), dtype=np.int64),
            chan=np.asarray(hits.get("hit_chan", []), dtype=np.int64),
            energy=np.asarray(hits.get("hit_energy", []), dtype=float),
            hg=np.asarray(hits.get("hit_hg", []), dtype=float),
            lg=np.asarray(hits.get("hit_lg", []), dtype=float),
            metrics=row,
        )

    def accumulate(self, indices) -> Event:
        """Aggregate the hits of many events into one pseudo-:class:`Event`.

        Hits from all ``indices`` are pooled and summed per pad (same ``slab`` and
        rounded ``x``/``y``), so the resulting "event" shows the *accumulated*
        energy deposited in each pad across the whole set -- the average shower
        footprint of, e.g., a cluster. ``chip``/``chan`` are not meaningful for an
        aggregate and are set to ``-1``.
        """
        hits = self.reader.all_hits()
        indices = list(indices)
        if not indices or "hit_x" not in hits or "hit_y" not in hits:
            empty = np.empty(0)
            return Event(index=-1, x=empty, y=empty, z=empty,
                         slab=empty.astype(np.int64), chip=empty.astype(np.int64),
                         chan=empty.astype(np.int64), energy=empty,
                         hg=empty, lg=empty, metrics={})
        slab = np.concatenate([np.asarray(hits["hit_slab"][i]) for i in indices])
        x = np.concatenate([np.asarray(hits["hit_x"][i]) for i in indices])
        y = np.concatenate([np.asarray(hits["hit_y"][i]) for i in indices])
        energy = np.concatenate([np.asarray(hits["hit_energy"][i]) for i in indices])

        pooled = pd.DataFrame({"slab": slab.astype(int),
                               "xr": np.round(x, 2), "yr": np.round(y, 2),
                               "e": energy.astype(float)})
        agg = pooled.groupby(["slab", "xr", "yr"], as_index=False)["e"].sum()
        slab_g = agg["slab"].to_numpy(dtype=np.int64)
        return Event(
            index=-1,
            x=agg["xr"].to_numpy(dtype=float),
            y=agg["yr"].to_numpy(dtype=float),
            z=self.detector.slab_z_array(slab_g),
            slab=slab_g,
            chip=np.full(slab_g.size, -1, dtype=np.int64),
            chan=np.full(slab_g.size, -1, dtype=np.int64),
            energy=agg["e"].to_numpy(dtype=float),
            hg=np.empty(0), lg=np.empty(0), metrics={},
        )
