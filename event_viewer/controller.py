"""
Session controller: the single object the Dash callbacks talk to.

It owns the heavy, non-serialisable state -- the detector geometry, the figure
builders, and a cache of :class:`EventDataset` keyed by file path -- so that the
``dcc.Store`` values can stay light (just the current path, event position, cuts
and cluster labels). Switching files at runtime is just asking the controller for
a different dataset; the previous one stays cached.
"""

from __future__ import annotations

import glob
import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from .analysis.clustering import UNCLUSTERED, ClusteringService
from .analysis.cuts import CutModel
from .config import ViewerConfig
from .io import make_reader
from .model import DetectorModel, Event, EventDataset
from .viz import DetectorScene3D, DistributionPlots, LayerGrid2D


def _filter_event(event: Event, mask: np.ndarray) -> Event:
    """Return a new Event with per-hit arrays filtered by boolean mask."""
    def _m(arr):
        return arr[mask] if arr.size == mask.size else arr
    return Event(
        index=event.index,
        x=_m(event.x), y=_m(event.y), z=_m(event.z),
        slab=_m(event.slab), chip=_m(event.chip), chan=_m(event.chan),
        energy=_m(event.energy), hg=_m(event.hg), lg=_m(event.lg),
        metrics=event.metrics,
    )


class ViewerController:
    """Holds geometry, builders and a per-path dataset cache."""

    def __init__(self, config: ViewerConfig):
        self.config = config
        self.detector = DetectorModel.from_config(config)
        self.scene3d = DetectorScene3D(self.detector, config.colorscale)
        self.layers2d = LayerGrid2D(self.detector, config.colorscale)
        self.distributions = DistributionPlots()
        self.clustering = ClusteringService()
        self._datasets: Dict[str, EventDataset] = {}
        self._accum_cache: Dict = {}   # (path, token, label) -> accumulated Event
        self._table_cache: Dict = {}   # (path, threshold) -> pd.DataFrame

    # ----------------------------------------------------------- file mgmt --
    def list_files(self) -> List[str]:
        """All ``*.valcache.root`` then ``ecal_*.root`` files under every data dir.

        The ``ecal_*.root`` pattern also matches the ``k4SiWEcalReco`` outputs
        (``ecal_<run>.edm4hep.root`` / ``ecal_<run>.valtree.root``), which
        ``make_reader`` then reads with the right reader. Scans each configured
        root (``settings.yml`` data roots + the optional cache dir).
        """
        seen, files = set(), []
        for data_dir in self.config.data_dirs:
            patterns = [
                os.path.join(data_dir, "**", "*.valcache.root"),
                os.path.join(data_dir, "**", "ecal_*.root"),
            ]
            for pattern in patterns:
                for path in sorted(glob.glob(pattern, recursive=True)):
                    if path.endswith(".valcache.root") or ".valcache." not in path:
                        if path not in seen:
                            seen.add(path)
                            files.append(path)
        return files

    def dataset(self, path: str) -> EventDataset:
        """Return (and cache) the :class:`EventDataset` for ``path``."""
        if path not in self._datasets:
            reader = make_reader(path, self.config.tree_name, self.config.n_layers)
            self._datasets[path] = EventDataset(reader, self.detector)
        return self._datasets[path]

    # --------------------------------------------------- hit-energy threshold --
    def _filtered_table(self, path: str, threshold: float):
        """Per-event DataFrame with metrics for hit_energy >= threshold.

        At threshold = 0.0 returns the original cached table.
        For non-zero thresholds uses pre-computed branches in the valcache
        (instant); falls back to in-memory recompute only for plain ecal files
        that have no pre-computed branches.
        """
        if threshold <= 0.0:
            return self.dataset(path).table
        key = (path, round(float(threshold), 4))
        if key in self._table_cache:
            return self._table_cache[key]
        reader = self.dataset(path).reader
        if reader.has_mip_thresholds:
            df = reader.event_table(float(threshold))
        else:
            # Fallback: in-memory recompute (slow; only for ecal files with no
            # pre-computed branches).
            from .analysis.recompute import recompute_all_metrics
            df = recompute_all_metrics(
                reader, self.detector.w_thickness_mm,
                float(threshold), self.config.n_layers)
        # Only a subset of columns is recomputed per threshold; the variable
        # menus, however, are built from the full threshold-0 table. Carry over
        # any base columns the threshold table lacks (e.g. raw branches like
        # ``nhit_slab``) so every selectable variable stays plottable. Rows are
        # positionally aligned (one per event, same order/length).
        base = self.dataset(path).table
        missing = [c for c in base.columns if c not in df.columns]
        if missing:
            df = pd.concat([df.reset_index(drop=True),
                            base[missing].reset_index(drop=True)], axis=1)
        self._table_cache[key] = df
        return df

    # ------------------------------------------------------- event figures --
    def event_figures(self, path: str, index: int, color_clip: bool,
                      hit_threshold: float = 0.0, show_moliere: bool = False,
                      show_axis: bool = False):
        """``(scene3d_fig, layers2d_fig, metrics_rows)`` for one event."""
        event = self.dataset(path).get_event(index)
        if hit_threshold > 0.0:
            mask = event.energy >= hit_threshold
            event = _filter_event(event, mask)
            table = self._filtered_table(path, hit_threshold)
            event.metrics = (table.iloc[index].to_dict()
                             if index < len(table) else {})
        scene = self.scene3d.event_figure(
            event, color_clip, show_moliere=show_moliere, show_axis=show_axis)
        layers = self.layers2d.build(event, color_clip)
        rows = self._metric_rows(event)
        return scene, layers, rows

    def _metric_rows(self, event) -> List[dict]:
        """Per-event scalar metrics as ``[{variable, value}]`` for the table."""
        rows = []
        for key, value in event.metrics.items():
            if isinstance(value, (list, tuple, np.ndarray)):
                continue  # skip per-layer vectors
            if isinstance(value, float):
                shown = f"{value:.4g}"
            else:
                shown = str(value)
            rows.append({"variable": key, "value": shown})
        return rows

    # ------------------------------------------------------------- passing --
    def passing_indices(self, path: str, cut_model: CutModel,
                        hit_threshold: float = 0.0) -> np.ndarray:
        if hit_threshold <= 0.0:
            return self.dataset(path).passing(cut_model)
        df = self._filtered_table(path, hit_threshold)
        if cut_model is None or cut_model.is_empty:
            return np.arange(len(df))
        return cut_model.passing_indices(df)

    # -------------------------------------------------------- distributions --
    def histogram(self, path: str, variable: str, cut_model: CutModel,
                  nbins: int = 60, cluster=None, hit_threshold: float = 0.0):
        """Histogram of ``variable``; stacked by cluster when ``cluster`` given.

        Without ``cluster`` the histogram covers the events passing the current
        cuts. With ``cluster`` (a ``{passing, labels}`` snapshot) the histogram is
        built over exactly those clustered events and split per label, so the
        cluster separation matches the scatter.
        """
        df = self._filtered_table(path, hit_threshold)
        if variable not in df.columns:
            return self.distributions.histogram(np.empty(0), variable)

        cut_range = None
        for cut in (cut_model.cuts if cut_model else []):
            if cut.variable == variable:
                cut_range = (cut.lo, cut.hi)

        if cluster is not None:
            sub = df.iloc[cluster["passing"]]
            values = sub[variable].to_numpy(dtype=float)
            return self.distributions.histogram(
                values, variable, cut_range, nbins, labels=cluster["labels"])

        keep = cut_model.mask(df) if cut_model and not cut_model.is_empty \
            else np.ones(len(df), bool)
        values = df.loc[keep, variable].to_numpy(dtype=float)
        return self.distributions.histogram(values, variable, cut_range, nbins)

    def histogram_split(self, path: str, variable: str, cut_model: CutModel,
                        nbins: int = 60, hit_threshold: float = 0.0):
        """Full distribution of ``variable`` coloured by the cut decision.

        Drives the Event tab: the whole distribution is shown, with the events
        passing *every* cut highlighted and the rest muted, so the effect of the
        cuts (which also limit the one-by-one navigation) is visible without
        hiding the removed events.
        """
        df = self._filtered_table(path, hit_threshold)
        if variable not in df.columns:
            return self.distributions.histogram_split(
                np.empty(0), np.empty(0, bool), variable, nbins)
        keep = cut_model.mask(df) if cut_model and not cut_model.is_empty \
            else np.ones(len(df), bool)
        values = df[variable].to_numpy(dtype=float)
        return self.distributions.histogram_split(values, keep, variable, nbins)

    def variable_range(self, path: str, variable: str,
                       hit_threshold: float = 0.0):
        """``(min, max)`` of a finite variable, for slider bounds."""
        df = self._filtered_table(path, hit_threshold)
        if variable not in df.columns:
            return 0.0, 1.0
        col = df[variable].to_numpy(dtype=float)
        col = col[np.isfinite(col)]
        if col.size == 0:
            return 0.0, 1.0
        return float(col.min()), float(col.max())

    # ----------------------------------------------------------- clustering --
    def run_clustering(self, path: str, cut_model: CutModel, features: List[str],
                       algo: str, n_clusters: int, eps: float, min_samples: int,
                       hit_threshold: float = 0.0):
        """Cluster the passing events; return ``(passing_indices, labels)``."""
        df = self._filtered_table(path, hit_threshold)
        if cut_model is None or cut_model.is_empty:
            keep = np.arange(len(df))
        else:
            keep = cut_model.passing_indices(df)
        sub = df.iloc[keep]
        labels = self.clustering.fit(sub, features, algo, n_clusters,
                                     eps, min_samples)
        return keep.tolist(), labels.tolist()

    def _accumulated_event(self, path: str, cluster, label: int):
        """Cached accumulated pseudo-event for one cluster of a run.

        Keyed by the cluster run ``token`` so repeated threshold tweaks reuse the
        (expensive) accumulation instead of recomputing it.
        """
        key = (path, cluster.get("token"), int(label))
        if key not in self._accum_cache:
            passing = np.asarray(cluster["passing"])
            labels = np.asarray(cluster["labels"])
            members = passing[labels == int(label)]
            self._accum_cache[key] = self.dataset(path).accumulate(members.tolist())
        return self._accum_cache[key]

    def cluster_panels(self, path: str, cluster):
        """``[(label, n_events, e_max), ...]`` for every group, incl. unclustered.

        The unclustered group (label ``-1``: DBSCAN noise, or events dropped from
        the fit because a feature was NaN) is included too -- it still has hits, so
        its accumulated profile is meaningful -- and is ordered last. ``e_max``
        (the largest accumulated pad energy) sizes each panel's threshold slider.
        Accumulations are computed (and cached) here.
        """
        labels = np.asarray(cluster["labels"])
        # Normal clusters first (ascending), the unclustered group last.
        order = sorted(set(labels.tolist()), key=lambda lab: (lab == UNCLUSTERED, lab))
        panels = []
        for lab in order:
            event = self._accumulated_event(path, cluster, lab)
            e_max = float(np.nanmax(event.energy)) if event.energy.size else 1.0
            panels.append((int(lab), int((labels == lab).sum()), e_max))
        return panels

    def cluster_scene(self, path: str, cluster, label: int, threshold: float = 0.0):
        """Accumulated 3-D scene for one cluster, fading hits below ``threshold``."""
        event = self._accumulated_event(path, cluster, int(label))
        thr = float(threshold) if threshold else None
        return self.scene3d.event_figure(event, color_clip=True, threshold=thr)

    def cluster_scatter(self, path: str, xvar: str, yvar: str,
                        passing: Optional[List[int]], labels: Optional[List[int]],
                        cut_model: CutModel, hit_threshold: float = 0.0):
        """Scatter of two variables; colour by cluster label when available."""
        df = self._filtered_table(path, hit_threshold)
        if passing is not None and labels is not None:
            sub = df.iloc[passing]
            return self.distributions.scatter(sub, xvar, yvar, labels)
        keep = self.passing_indices(path, cut_model, hit_threshold)
        return self.distributions.scatter(df.iloc[keep], xvar, yvar, None)
