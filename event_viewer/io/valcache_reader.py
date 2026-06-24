"""
Read a reconstructed ``ecal`` ROOT file (optionally a ``.valcache.root``) with
uproot, exposing a clean separation between *per-event* and *per-hit* branches.

``EventFileReader`` is deliberately thin: it opens the file, classifies the
branches, builds the per-event table once (lazily), and reads the per-hit arrays
of a single event on demand. It does **not** know about geometry, plotting or
Dash; those layers consume the numpy/pandas it returns.

Branch classification (by uproot interpretation)
------------------------------------------------
* per-hit  : jagged arrays of length ``nhit_chan`` (``hit_x``, ``hit_energy``, …).
* per-layer: fixed-length ``double[15]`` profiles (``energy_per_layer``, …).
* scalar   : one value per event (identifiers + per-event metrics).

A file "has metrics" when the marker branch ``mip_likeness`` is present (written
by the validation cache). Plain ``ecal_*.root`` files lack it; the reader then
augments its scalar table with cheap derivable quantities computed from the
per-hit arrays (see :mod:`event_viewer.analysis.derived`).
"""

from __future__ import annotations

from typing import List, Optional

import pandas as pd
import uproot

# Branch present only in the validation cache; used as the "has metrics" marker.
METRICS_MARKER = "mip_likeness"

# Identifier-like scalar columns that are never useful as cut/cluster variables.
IDENTIFIER_COLUMNS = ("run", "event", "spill", "bcid")

# Pre-computed MIP-cut threshold levels: float → branch-name prefix.
_THRESHOLD_PREFIXES = {0.5: "mip05", 1.0: "mip1"}


class EventFileReader:
    """Lazy reader over one ``ecal`` tree.

    Parameters
    ----------
    path : str
        Path to the ROOT file.
    tree_name : str
        Name of the TTree to read (default ``"ecal"``).
    n_layers : int
        Number of detector layers, used when computing cheap metrics.
    """

    def __init__(self, path: str, tree_name: str = "ecal", n_layers: int = 15):
        self.path = path
        self.tree_name = tree_name
        self.n_layers = n_layers

        self._file = uproot.open(path)
        if tree_name not in self._file:
            keys = [k.split(";")[0] for k in self._file.keys()]
            raise KeyError(
                f"Tree '{tree_name}' not found in {path}. Available: {keys}")
        self.tree = self._file[tree_name]

        self._perhit, self._perlayer, self._scalar = self._classify()
        self.has_metrics = METRICS_MARKER in self._scalar
        # True when the valcache contains pre-computed mip05_* branches.
        self.has_mip_thresholds = "mip05_nhit" in self._scalar
        self._tables: dict = {}          # threshold (float) -> pd.DataFrame
        self._all_hits: Optional[dict] = None        # built lazily

    # ------------------------------------------------------------- metadata --
    @property
    def n_events(self) -> int:
        return int(self.tree.num_entries)

    @property
    def perhit_branches(self) -> List[str]:
        return list(self._perhit)

    @property
    def scalar_branches(self) -> List[str]:
        return list(self._scalar)

    def _classify(self):
        """Split branch names into (per-hit, per-layer, scalar) lists."""
        perhit, perlayer, scalar = [], [], []
        for name in self.tree.keys():
            branch = self.tree[name]
            interp_kind = type(branch.interpretation).__name__
            if interp_kind == "AsJagged":
                perhit.append(name)
            elif "[" in branch.typename:        # fixed-length array, e.g. double[15]
                perlayer.append(name)
            else:
                scalar.append(name)
        return perhit, perlayer, scalar

    # ---------------------------------------------------------- event table --
    def _threshold_scalar_names(self, prefix: str) -> List[str]:
        """Branch names in the tree that carry the given mip prefix."""
        return [n for n in self._scalar if n.startswith(f"{prefix}_")]

    def event_table(self, threshold: float = 0.0) -> pd.DataFrame:
        """One row per event with scalar metrics for the given hit_energy threshold.

        * threshold = 0.0 (default): all scalar branches from the tree; for
          files without the validation cache, cheap derivable quantities are
          appended on-the-fly.
        * threshold > 0 and ``has_mip_thresholds``: reads the pre-computed
          prefixed branches (e.g. ``mip05_*``) and returns them renamed to
          the standard names so the rest of the viewer is unaware of the
          threshold.  Identifier columns are taken from the threshold-0 table.

        Results are cached per threshold.
        """
        key = round(float(threshold), 4)
        if key in self._tables:
            return self._tables[key]

        if threshold <= 0.0 or not self.has_mip_thresholds:
            # Standard path (threshold = 0 or no pre-computed branches).
            # Exclude mip-prefix branches so they don't pollute the feature list.
            base_scalar = [n for n in self._scalar
                           if not any(n.startswith(f"{p}_")
                                      for p in _THRESHOLD_PREFIXES.values())]
            df = self.tree.arrays(base_scalar, library="pd")
            if not self.has_metrics:
                from ..analysis.derived import compute_cheap_metrics
                cheap = compute_cheap_metrics(self)
                new_cols = [c for c in cheap.columns if c not in df.columns]
                df = pd.concat([df, cheap[new_cols]], axis=1)
            self._tables[key] = df
            return df

        # threshold > 0: read prefixed branches and rename to standard names.
        prefix = _THRESHOLD_PREFIXES.get(key)
        if prefix is None:
            # Unknown threshold level — fall back to threshold=0.
            return self.event_table(0.0)

        prefixed = self._threshold_scalar_names(prefix)
        if not prefixed:
            return self.event_table(0.0)

        thr_df = self.tree.arrays(prefixed, library="pd")
        # Rename mipXX_name → name.
        rename = {b: b[len(prefix) + 1:] for b in prefixed}
        thr_df = thr_df.rename(columns=rename)

        # Merge identifier columns from the threshold-0 table.
        base_df = self.event_table(0.0)
        id_cols = [c for c in IDENTIFIER_COLUMNS if c in base_df.columns
                   and c not in thr_df.columns]
        if id_cols:
            thr_df = pd.concat([base_df[id_cols].reset_index(drop=True),
                                 thr_df.reset_index(drop=True)], axis=1)

        self._tables[key] = thr_df
        return thr_df

    def feature_columns(self) -> List[str]:
        """Numeric scalar columns usable as cut / clustering variables.

        Drops pure identifiers and any non-numeric column.  Always derived
        from the threshold-0 table so the list is stable regardless of which
        threshold is currently active.
        """
        df = self.event_table(0.0)
        cols = []
        for name in df.columns:
            if name in IDENTIFIER_COLUMNS:
                continue
            if pd.api.types.is_numeric_dtype(df[name]) or \
                    pd.api.types.is_bool_dtype(df[name]):
                cols.append(name)
        return cols

    # ------------------------------------------------------------ per event --
    def read_hits(self, index: int) -> dict:
        """Per-hit numpy arrays for a single event ``index``.

        Returns a dict keyed by branch name (``hit_x``, ``hit_energy``, …). Empty
        events yield empty arrays.
        """
        arrays = self.tree.arrays(
            self._perhit, entry_start=index, entry_stop=index + 1, library="np")
        return {name: arrays[name][0] for name in self._perhit}

    def all_hits(self) -> dict:
        """All events' position/energy per-hit arrays at once (cached).

        Returns a dict of jagged numpy object-arrays (one entry per event) for
        the branches needed to accumulate hits across many events. Reading the
        whole columns once is far cheaper than per-event reads when aggregating
        thousands of events (e.g. for the per-cluster accumulated views).
        """
        if self._all_hits is None:
            names = [n for n in ("hit_slab", "hit_x", "hit_y", "hit_energy")
                     if n in self._perhit]
            self._all_hits = self.tree.arrays(names, library="np")
        return self._all_hits

    def close(self) -> None:
        self._file.close()
