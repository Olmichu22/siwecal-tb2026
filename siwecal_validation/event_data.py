"""
Per-event variables read from the reconstructed ``ecal`` tree.

:class:`EventData` is the in-memory representation of one sample (a run or an
energy point): a set of aligned numpy arrays, one entry per accepted event. It
knows how to load itself from a ROOT file (computing every metric via
:mod:`siwecal_validation.metrics`) and how to return a filtered copy given a
:class:`~siwecal_validation.selection.CutSet`.

All per-event scalars and per-layer profiles are computed in a single pass over
the tree; the physics lives in :mod:`metrics`, so this module is just I/O and
bookkeeping.
"""

from dataclasses import dataclass, fields, replace

import numpy as np
import ROOT

from . import metrics
from . import vars_cache


def _append_scalar_metrics(ecols, slab, energy, x, y,
                           w_over_x0, n_layers, profile_kind,
                           shower_thr, peak_min, start_frac, containment,
                           hit_threshold):
    """Compute per-event scalar metrics after filtering hits below hit_threshold
    and append one value to each list in ecols.  Called inside the main event
    loop for each extra MIP-cut threshold level.
    """
    NAN = float("nan")
    mask = energy >= hit_threshold
    s, e = slab[mask], energy[mask]
    fx, fy = x[mask], y[mask]

    if s.size == 0:
        ecols["nhit"].append(0)
        ecols["energy"].append(0.0)
        ecols["zbary"].append(NAN)
        ecols["mip_likeness"].append(0.0)
        ecols["weighte"].append(0.0)
        ecols["bar_x"].append(NAN)
        ecols["bar_y"].append(NAN)
        ecols["bar_r"].append(NAN)
        ecols["moliere"].append(NAN)
        ecols["transverse_rms"].append(NAN)
        ecols["is_shower"].append(False)
        ecols["shower_start"].append(NAN)
        ecols["shower_max"].append(NAN)
        ecols["shower_end"].append(NAN)
        ecols["shower_start_10"].append(NAN)
        ecols["shower_end_10"].append(NAN)
        ecols["shower_length"].append(NAN)
        ecols["first_layer"].append(NAN)
        ecols["last_layer"].append(NAN)
        ecols["n_layers_hit"].append(0)
        ecols["e_over_nhit"].append(NAN)
        return

    hits_layer = metrics.hits_per_layer(s, n_layers)
    e_layer    = metrics.energy_per_layer(s, e, n_layers)
    w_layer    = metrics.weighte_per_layer(s, e, w_over_x0, n_layers)
    first, last, n_layers_hit = metrics.layer_extent(hits_layer)

    pos = e > 0
    pw = metrics.hit_weights(s[pos], e[pos], w_over_x0)
    bx, by, br = metrics.barycenter_xy(fx[pos], fy[pos], pw)

    profile = {"nhit": hits_layer, "sume": e_layer, "weighte": w_layer}[profile_kind]
    sh = metrics.shower_features(profile, shower_thr, peak_min, start_frac)

    moliere = (metrics.moliere_radius(fx[pos], fy[pos], pw, bx, by, containment)
               if sh.is_shower and np.isfinite(bx) else 0.0)

    ecols["nhit"].append(s.size)
    ecols["energy"].append(float(e.sum()))
    ecols["zbary"].append(metrics.zbary(s, e))
    ecols["mip_likeness"].append(metrics.mip_likeness(hits_layer, n_layers))
    ecols["weighte"].append(metrics.weighte_total(s, e, w_over_x0))
    ecols["bar_x"].append(bx)
    ecols["bar_y"].append(by)
    ecols["bar_r"].append(br)
    ecols["moliere"].append(moliere)
    ecols["transverse_rms"].append(
        metrics.transverse_rms(fx[pos], fy[pos], pw, bx, by)
        if np.isfinite(bx) else NAN)
    ecols["is_shower"].append(sh.is_shower)
    ecols["shower_start"].append(sh.start_layer)
    ecols["shower_max"].append(sh.max_layer)
    ecols["shower_end"].append(sh.end_layer)
    ecols["shower_start_10"].append(sh.start_layer_10)
    ecols["shower_end_10"].append(sh.end_layer_10)
    ecols["shower_length"].append(sh.length)
    ecols["first_layer"].append(first)
    ecols["last_layer"].append(last)
    ecols["n_layers_hit"].append(n_layers_hit)
    ecols["e_over_nhit"].append(float(e.sum()) / s.size)


@dataclass
class EventData:
    """Aligned per-event arrays for one sample."""

    label: str

    # --- original variables (unchanged) ------------------------------------
    nhit: np.ndarray            # total number of hit channels per event
    zbary: np.ndarray           # energy-weighted mean slab index (layer units)
    energy: np.ndarray          # summed hit energy per event (MIP units)
    hits_per_layer: np.ndarray  # (n_events, n_layers): hits per layer
    mip_likeness: np.ndarray    # inverse-hit-rate MIP-likeness score per event

    # --- particle-discrimination variables ---------------------------------
    energy_per_layer: np.ndarray    # (n_events, n_layers): Σ E per layer
    weighte: np.ndarray             # tungsten-weighted energy Σ E·W/X0
    weighte_per_layer: np.ndarray   # (n_events, n_layers): weighted E per layer
    bar_x: np.ndarray               # energy-weighted transverse barycenter x
    bar_y: np.ndarray               # energy-weighted transverse barycenter y
    bar_r: np.ndarray               # sqrt(bar_x² + bar_y²)
    moliere: np.ndarray             # 90% transverse containment radius (showers)
    transverse_rms: np.ndarray      # energy-weighted RMS hit radius
    is_shower: np.ndarray           # bool: longitudinal shower flag
    shower_start: np.ndarray        # first shower layer (NaN if not a shower)
    shower_max: np.ndarray          # peak layer
    shower_end: np.ndarray          # last shower layer
    shower_start_10: np.ndarray     # first layer above 10% of the peak
    shower_end_10: np.ndarray       # last layer above 10% of the peak
    shower_length: np.ndarray       # number of layers above threshold
    first_layer: np.ndarray         # lowest layer with a hit
    last_layer: np.ndarray          # highest layer with a hit
    n_layers_hit: np.ndarray        # number of layers with at least one hit
    e_over_nhit: np.ndarray         # energy / nhit (hit energy density)

    # --- provenance (EDM4hep path only) ------------------------------------
    #: Original podio frame index of each event (set by :meth:`from_edm4hep`;
    #: ``None`` for the ``from_root`` path). Lets ``--save-tree`` map a cut mask
    #: back to the frames to copy into a filtered EDM4hep file. Sliced by
    #: :meth:`select` like any other per-event array.
    source_index: np.ndarray = None

    #: Dataclass fields that are not per-event metric arrays.
    _META_FIELDS = ("label", "source_index")

    def __len__(self) -> int:
        return len(self.nhit)

    # ------------------------------------------------------------- loading --
    @classmethod
    def from_root(cls, path: str, label: str, config,
                  create_tree: bool = False, cache_dir: str = None) -> "EventData":
        """Read the ``ecal`` tree and compute every per-event metric.

        Events with no hits or with non-positive total energy are skipped (as
        before). The original variables (``nhit``, ``zbary``, ``energy``,
        ``hits_per_layer``, ``mip_likeness``) are preserved exactly; the rest are
        the discrimination metrics from :mod:`metrics`. ``config`` is a
        :class:`~siwecal_validation.config.PlotConfig` (layers, tree name,
        tungsten map and shower thresholds).

        Caching: a ``<stem>.valcache.root`` stores every derived variable (see
        :mod:`vars_cache`). It is written next to the input by default, or into
        ``cache_dir`` if given. When it exists and matches the current config it
        is read back instead of recomputing (the fast path); ``create_tree=True``
        forces a fresh recomputation and rewrite.
        """
        cache = vars_cache.cache_path_for(path, cache_dir=cache_dir)
        if not create_tree and vars_cache.is_valid(cache, config):
            return vars_cache.read(cache, label, config)

        n_layers = config.n_layers
        w_over_x0 = config.w_over_x0()
        profile_kind = config.shower_profile
        thr, peak_min = config.shower_e_threshold, config.shower_max_min
        start_frac = config.shower_start_frac
        containment = config.moliere_containment

        root_file = ROOT.TFile(path, "READ")
        tree = root_file.Get(config.tree_name)
        if not tree:
            keys = [k.GetName() for k in root_file.GetListOfKeys()]
            root_file.Close()
            raise RuntimeError(f"Tree '{config.tree_name}' not found in {path}. "
                               f"Available: {keys}")

        # Channels with no MIP calibration are flagged ``hit_ismasked`` by the
        # event builder and excluded from every metric. Older ``ecal`` files
        # predate the branch; there we keep all hits (no masking available).
        has_ismasked = bool(tree.GetListOfBranches().FindObject("hit_ismasked"))

        # One list per EventData array field (filled per accepted event).
        cols = {f.name: [] for f in fields(cls) if f.name not in cls._META_FIELDS}

        # Extra MIP-cut threshold columns: {threshold: {scalar_name: list}}.
        from .vars_cache import MIP_CUT_THRESHOLDS, _scalar_derived_names
        _scalar_names = _scalar_derived_names()
        extra_cols = {t: {n: [] for n in _scalar_names}
                      for t in MIP_CUT_THRESHOLDS}

        for entry_index in range(tree.GetEntries()):
            tree.GetEntry(entry_index)
            n_channels = int(tree.nhit_chan)
            if n_channels == 0:
                continue
            slab = np.frombuffer(tree.hit_slab, np.int32, count=n_channels).copy()
            energy = np.frombuffer(tree.hit_energy, np.float32,
                                   count=n_channels).copy()
            x = np.frombuffer(tree.hit_x, np.float32, count=n_channels).copy()
            y = np.frombuffer(tree.hit_y, np.float32, count=n_channels).copy()

            # Drop masked (uncalibrated) channels before any metric is computed.
            if has_ismasked:
                ismasked = np.frombuffer(tree.hit_ismasked, np.int32,
                                         count=n_channels).astype(bool)
                keep = ~ismasked
                slab, energy, x, y = slab[keep], energy[keep], x[keep], y[keep]
            n_channels = int(slab.size)
            if n_channels == 0:
                continue
            energy_sum = float(energy.sum())
            if not (energy_sum > 0):
                continue

            # Per-layer profiles.
            hits_layer = metrics.hits_per_layer(slab, n_layers)
            e_layer = metrics.energy_per_layer(slab, energy, n_layers)
            w_layer = metrics.weighte_per_layer(slab, energy, w_over_x0, n_layers)

            # Transverse moments (barycenter, RMS, Molière) use only hits with
            # positive energy: in real data pedestal subtraction yields negative
            # hit energies (noise) which would otherwise bias the centroid and
            # break the energy-ordered Molière cumulant. The longitudinal sums
            # above keep the signed (calorimetric) energy.
            pos = energy > 0
            px, py = x[pos], y[pos]
            pw = metrics.hit_weights(slab[pos], energy[pos], w_over_x0)
            bx, by, br = metrics.barycenter_xy(px, py, pw)

            # Longitudinal shower descriptors from the configured profile.
            profile = {"nhit": hits_layer, "sume": e_layer,
                       "weighte": w_layer}[profile_kind]
            sh = metrics.shower_features(profile, thr, peak_min, start_frac)

            # Molière radius only for showers (mirrors the template).
            moliere = (metrics.moliere_radius(px, py, pw, bx, by, containment)
                       if sh.is_shower else 0.0)
            first, last, n_layers_hit = metrics.layer_extent(hits_layer)

            cols["nhit"].append(n_channels)
            cols["zbary"].append(metrics.zbary(slab, energy))
            cols["energy"].append(energy_sum)
            cols["hits_per_layer"].append(hits_layer)
            cols["mip_likeness"].append(metrics.mip_likeness(hits_layer, n_layers))
            cols["energy_per_layer"].append(e_layer)
            cols["weighte"].append(metrics.weighte_total(slab, energy, w_over_x0))
            cols["weighte_per_layer"].append(w_layer)
            cols["bar_x"].append(bx)
            cols["bar_y"].append(by)
            cols["bar_r"].append(br)
            cols["moliere"].append(moliere)
            cols["transverse_rms"].append(
                metrics.transverse_rms(px, py, pw, bx, by))
            cols["is_shower"].append(sh.is_shower)
            cols["shower_start"].append(sh.start_layer)
            cols["shower_max"].append(sh.max_layer)
            cols["shower_end"].append(sh.end_layer)
            cols["shower_start_10"].append(sh.start_layer_10)
            cols["shower_end_10"].append(sh.end_layer_10)
            cols["shower_length"].append(sh.length)
            cols["first_layer"].append(first)
            cols["last_layer"].append(last)
            cols["n_layers_hit"].append(n_layers_hit)
            cols["e_over_nhit"].append(energy_sum / n_channels)

            # --- extra MIP-threshold metrics (same pass, filtered hits) ---
            for mip_thr, ecols in extra_cols.items():
                _append_scalar_metrics(
                    ecols, slab, energy, x, y,
                    w_over_x0, n_layers, profile_kind,
                    thr, peak_min, start_frac, containment, mip_thr)

        root_file.Close()

        arrays = {name: np.array(values,
                                 dtype=bool if name == "is_shower" else float)
                  for name, values in cols.items()}
        data = cls(label=label, **arrays)

        extra_scalars = {
            t: {n: np.array(v, dtype=bool if n == "is_shower" else float)
                for n, v in ecols.items()}
            for t, ecols in extra_cols.items()
        }

        # Persist the augmented cache for fast re-runs (best-effort: a read-only
        # input directory just means no cache, not a failure).
        if create_tree or not vars_cache.is_valid(cache, config):
            try:
                vars_cache.write(path, cache, config, data,
                                 extra_scalars=extra_scalars)
                print(f"[Cache] wrote derived-variable tree -> {cache}")
            except OSError as error:
                print(f"WARNING: could not write metrics cache {cache}: {error}")
        return data

    # ----------------------------------------------- loading (EDM4hep) ------
    #: EventData fields stored as 2-D per-layer profiles (length ``n_layers``).
    _PER_LAYER_FIELDS = ("hits_per_layer", "energy_per_layer", "weighte_per_layer")

    @classmethod
    def from_edm4hep(cls, path: str, label: str, config,
                     mip_thresholds=(0.5, 1.0)) -> "EventData":
        """Build an :class:`EventData` from a ``k4SiWEcalReco`` EDM4hep file.

        The per-event discrimination variables are read straight from each
        ``Cluster``'s ``shapeParameters`` (computed in C++ by
        ``EcalPidTransformer``); nothing is recomputed here. Events with no hits
        or non-positive energy are dropped, matching :meth:`from_root`.

        ``mip_thresholds`` must match the file's ``shapeParameters`` layout
        (empty when the PID stage ran in physics mode without the MIP-variant
        blocks); it only affects the width check, not the cut variables read.
        """
        from siwecal_common.edm4hep_pid import PidFileReader

        reader = PidFileReader(path, n_layers=config.n_layers,
                               mip_thresholds=mip_thresholds)
        cols = reader.scalar_columns()
        ids = reader.identifiers()
        n_layers = config.n_layers

        valid = (ids["nhit_chan"] > 0) & (cols["energy"] > 0)

        scalar_fields = [f.name for f in fields(cls)
                         if f.name not in cls._META_FIELDS
                         and f.name not in cls._PER_LAYER_FIELDS]
        arrays = {}
        for name in scalar_fields:
            if name not in cols:
                raise RuntimeError(
                    f"Variable '{name}' not found in {path}. Available shape "
                    f"parameters: {reader.shape_names[:5]}...")
            values = cols[name][valid]
            arrays[name] = values.astype(bool) if name == "is_shower" else values.astype(float)

        for name in cls._PER_LAYER_FIELDS:
            profile = np.stack([cols[f"{name}_{i}"] for i in range(n_layers)], axis=1)
            arrays[name] = profile[valid].astype(float)

        # Original podio frame index of each surviving event, so --save-tree can
        # copy exactly the cut-passing frames into a filtered EDM4hep file.
        arrays["source_index"] = np.nonzero(valid)[0].astype(np.int64)

        return cls(label=label, **arrays)

    # ----------------------------------------------------------- filtering --
    def select(self, cutset) -> "EventData":
        """Return a new EventData containing only events passing ``cutset``.

        Every numpy-array field is sliced by the same boolean mask, so new
        per-event variables are filtered automatically.
        """
        if cutset.is_empty:
            return self
        keep = cutset.mask(self)
        sliced = {f.name: getattr(self, f.name)[keep]
                  for f in fields(self)
                  if isinstance(getattr(self, f.name), np.ndarray)}
        return replace(self, **sliced)
