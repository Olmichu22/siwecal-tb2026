"""
Reader over a ``k4SiWEcalReco`` EDM4hep file, exposing the same API as
:class:`~event_viewer.io.valcache_reader.EventFileReader` so the rest of the
viewer is unchanged.

Per-event metrics are read straight from each ``Cluster``'s ``shapeParameters``
(no recomputation); per-hit arrays come from the ``CalorimeterHit`` collection
plus the parallel ``UserDataCollection``s (chip/chan/sca/hg/lg). All of this is
delegated to :class:`siwecal_common.edm4hep_pid.PidFileReader`; this class only
shapes it into the per-event ``pandas`` table the viewer expects.
"""

from __future__ import annotations

from typing import List

import pandas as pd

from siwecal_common.edm4hep_pid import PidFileReader

# Identifier-like columns: never offered as cut/cluster variables.
IDENTIFIER_COLUMNS = ("run", "event", "spill", "bcid", "nhit_chan")
# MIP hit-energy threshold levels -> shapeParameter name prefix.
_THRESHOLD_PREFIXES = {0.5: "mip05", 1.0: "mip1"}
# Stored as a boolean column (the rest of the scalars are floats).
_BOOL_COLUMNS = ("is_shower",)


class Edm4hepEventReader:
    """Lazy reader over one EDM4hep PID file (viewer-facing API)."""

    def __init__(self, path: str, tree_name: str = "events", n_layers: int = 15):
        self.path = path
        self.tree_name = tree_name
        self.n_layers = n_layers
        self._pid = PidFileReader(path, n_layers=n_layers)

        # Base per-event scalars = shape parameters that are neither per-layer
        # blocks nor MIP-cut variants. Note the MIP-cut prefixes are ``mip05_``/
        # ``mip1_`` -- do NOT exclude the base scalar ``mip_likeness``.
        _mip_cut = tuple(f"{p}_" for p in _THRESHOLD_PREFIXES.values())
        self._base_scalars = [
            n for n in self._pid.shape_names
            if "_per_layer_" not in n and not n.startswith(_mip_cut)]
        self.has_metrics = True
        self.has_mip_thresholds = any(n.startswith("mip05_")
                                      for n in self._pid.shape_names)
        self._tables: dict = {}

    # ------------------------------------------------------------- metadata --
    @property
    def n_events(self) -> int:
        return self._pid.n_events

    @property
    def perhit_branches(self) -> List[str]:
        from siwecal_common.edm4hep_pid import PERHIT_FIELDS
        return list(PERHIT_FIELDS)

    @property
    def scalar_branches(self) -> List[str]:
        return list(self._base_scalars) + list(IDENTIFIER_COLUMNS)

    # ---------------------------------------------------------- event table --
    def event_table(self, threshold: float = 0.0) -> pd.DataFrame:
        """One row per event with scalar metrics for the given hit_energy threshold.

        ``threshold = 0`` returns the base variables; ``0.5``/``1.0`` return the
        pre-computed ``mip05_*``/``mip1_*`` variants renamed to the base names so
        the rest of the viewer is unaware of the threshold.
        """
        key = round(float(threshold), 4)
        if key in self._tables:
            return self._tables[key]

        if key <= 0.0 or not self.has_mip_thresholds or key not in _THRESHOLD_PREFIXES:
            cols = self._pid.scalar_columns(self._base_scalars)
            df = self._frame(cols)
        else:
            prefix = _THRESHOLD_PREFIXES[key]
            mip_names = [f"{prefix}_{n}" for n in self._base_scalars]
            cols = self._pid.scalar_columns(mip_names)
            renamed = {n: cols[f"{prefix}_{n}"] for n in self._base_scalars}
            renamed.update(self._pid.identifiers())
            df = self._frame(renamed)

        self._tables[key] = df
        return df

    @staticmethod
    def _frame(cols: dict) -> pd.DataFrame:
        df = pd.DataFrame(cols)
        for name in _BOOL_COLUMNS:
            if name in df.columns:
                df[name] = df[name].astype(bool)
        return df

    def feature_columns(self) -> List[str]:
        """Numeric scalar columns usable as cut / clustering variables."""
        df = self.event_table(0.0)
        return [name for name in df.columns
                if name not in IDENTIFIER_COLUMNS
                and (pd.api.types.is_numeric_dtype(df[name])
                     or pd.api.types.is_bool_dtype(df[name]))]

    # ------------------------------------------------------------ per event --
    def read_hits(self, index: int) -> dict:
        return self._pid.read_hits(index)

    def all_hits(self) -> dict:
        return self._pid.all_hits()

    def close(self) -> None:
        self._pid.close()
