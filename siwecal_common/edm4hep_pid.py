"""
Reader for the EDM4hep output of ``k4SiWEcalReco`` (the Gaudi PID stage).

The PID job writes one podio event (Frame) per physics event with:

* ``ECalHits``    -- an ``edm4hep::CalorimeterHitCollection`` (per-hit energy,
  position, cellID), the EDM4hep replacement for the old per-hit ``hit_*``
  branches of the ``ecal`` tree.
* ``ECalPid``     -- an ``edm4hep::ClusterCollection`` with exactly one Cluster
  per event whose ``shapeParameters`` hold every per-event discrimination
  variable (the old ``*.valcache.root`` derived branches). The names are stored
  once in the ``metadata`` frame parameter ``ECalPid_shapeParameterNames``.
* ``EventHeader`` -- run / event numbers (+ bcid in ``timeStamp``, spill in
  ``weight``).

This is the single source of truth both ``siwecal_validation`` and
``event_viewer`` use to read EDM4hep, so neither recomputes the metrics: they are
read straight from the Cluster (computed in C++ by ``EcalPidTransformer``).
"""

from __future__ import annotations

from typing import Dict, List, Optional

import numpy as np

SHAPE_PARAM_META = "ECalPid_shapeParameterNames"
IDENTIFIER_COLUMNS = ("run", "event", "spill", "bcid", "nhit_chan")

# Canonical shapeParameters layout -- mirror of k4SiWEcalReco/EcalShowerVars.h
# (scalarNames / perLayerNames) and EcalPidTransformer's MIP-cut prefixes. Kept
# in Python so opening a file needs no podio call (the slow part); the actual
# width is asserted against the file in PidFileReader._build_scalars.
SCALAR_NAMES = (
    "nhit", "zbary", "energy", "mip_likeness", "weighte", "bar_x", "bar_y", "bar_r",
    "moliere", "transverse_rms", "is_shower", "shower_start", "shower_max", "shower_end",
    "shower_start_10", "shower_end_10", "shower_length", "first_layer", "last_layer",
    "n_layers_hit", "e_over_nhit")
PER_LAYER_NAMES = ("hits_per_layer", "energy_per_layer", "weighte_per_layer")
_MIP_PREFIX = {0.5: "mip05", 1.0: "mip1"}


def canonical_shape_names(n_layers: int = 15, mip_thresholds=(0.5, 1.0)) -> List[str]:
    """The ordered shapeParameter names written by ``EcalPidTransformer``."""
    names = list(SCALAR_NAMES)
    for block in PER_LAYER_NAMES:
        names += [f"{block}_{i}" for i in range(n_layers)]
    for thr in mip_thresholds:
        names += [f"{_MIP_PREFIX[thr]}_{s}" for s in SCALAR_NAMES]
    return names
# Per-hit fields exposed under the legacy ``hit_*`` names (so downstream is
# unchanged). slab/x/y/z/energy come from the CalorimeterHit; chip/chan/sca and
# hg/lg from the parallel podio UserDataCollections written by EcalToEDM4hep.
PERHIT_FIELDS = ("hit_slab", "hit_x", "hit_y", "hit_z", "hit_energy",
                 "hit_chip", "hit_chan", "hit_sca", "hit_hg", "hit_lg")

# Legacy hit field -> (collection name, numpy dtype) for the UserDataCollections.
_USERDATA_HITS = {
    "hit_chip": ("ECalHitChip", np.int32),
    "hit_chan": ("ECalHitChan", np.int32),
    "hit_sca": ("ECalHitSca", np.int32),
    "hit_hg": ("ECalHitHG", np.float32),
    "hit_lg": ("ECalHitLG", np.float32),
}


class PidFileReader:
    """Lazy reader over one EDM4hep PID file (``ECalHits`` + ``ECalPid``)."""

    def __init__(self, path: str, events_coll: str = "ECalPid",
                 hits_coll: str = "ECalHits", header_coll: str = "EventHeader",
                 n_layers: int = 15, mip_thresholds=(0.5, 1.0)):
        self.path = path
        self._events_coll = events_coll
        self._hits_coll = hits_coll
        self._header_coll = header_coll

        # Names from the Python schema -> opening a file does no podio work.
        self.shape_names: List[str] = canonical_shape_names(n_layers, mip_thresholds)
        self._name_idx = {n: i for i, n in enumerate(self.shape_names)}

        self._reader = None
        self._frames_cache = None
        self._scalars: Optional[np.ndarray] = None   # (n_events, n_shape_params)
        self._ids: Optional[Dict[str, np.ndarray]] = None
        self._hits: Optional[Dict[str, np.ndarray]] = None  # jagged object arrays

    @property
    def _frames(self):
        """podio frames, created lazily (only needed for per-hit access)."""
        if self._frames_cache is None:
            import podio.root_io as rio  # heavy: defer until a hit read is needed
            self._reader = rio.Reader(self.path)
            self._frames_cache = self._reader.get("events")
        return self._frames_cache

    # ----------------------------------------------------------- metadata ----
    @property
    def n_events(self) -> int:
        if self._scalars is None:
            self._build_scalars()
        return self._scalars.shape[0]

    def has_param(self, name: str) -> bool:
        return name in self._name_idx

    # ------------------------------------------------- per-event scalars -----
    def _build_scalars(self) -> None:
        """Vectorised read of the per-event scalars with uproot.

        The podio TTree stores the Cluster ``shapeParameters`` of all events in a
        single flat branch (one cluster per event -> a regular ``(n, n_params)``
        block) and the ``EventHeader`` fields as length-1 jagged branches.
        Reading them columnar with uproot is far faster than looping the podio
        frames in Python (which is what made the viewer slow to load).
        """
        import awkward as ak
        import uproot

        tree = uproot.open(self.path)["events"]
        ev = self._events_coll
        sp = tree[f"_{ev}_shapeParameters"].array()
        mat = ak.to_numpy(sp)                      # regular: 1 cluster/event
        if mat.shape[1] != len(self.shape_names):  # safety against a layout change
            raise RuntimeError(
                f"shapeParameters width {mat.shape[1]} != "
                f"{len(self.shape_names)} names in {self.path}")

        hd = self._header_coll
        first = lambda field: ak.to_numpy(  # noqa: E731
            ak.firsts(tree[f"{hd}/{hd}.{field}"].array()))
        self._scalars = mat.astype(float)
        self._ids = {
            "run": first("runNumber").astype(np.int64),
            "event": first("eventNumber").astype(np.int64),
            "bcid": first("timeStamp").astype(np.int64),
            "spill": first("weight").astype(np.int64),
            # nhit (cluster) == number of ECalHits == legacy nhit_chan.
            "nhit_chan": mat[:, self._name_idx["nhit"]].astype(np.int64),
        }

    def scalar(self, name: str) -> np.ndarray:
        """Per-event array of one shape-parameter (e.g. ``moliere``)."""
        if self._scalars is None:
            self._build_scalars()
        return self._scalars[:, self._name_idx[name]]

    def identifiers(self) -> Dict[str, np.ndarray]:
        if self._ids is None:
            self._build_scalars()
        return self._ids

    def scalar_columns(self, names: Optional[List[str]] = None) -> Dict[str, np.ndarray]:
        """Dict ``name -> per-event array`` for the requested shape-parameters
        (default: all of them), plus the identifier columns."""
        if self._scalars is None:
            self._build_scalars()
        names = names if names is not None else self.shape_names
        out = {n: self._scalars[:, self._name_idx[n]] for n in names if n in self._name_idx}
        out.update(self.identifiers())
        return out

    # ------------------------------------------------------- per-hit ---------
    def read_hits(self, index: int) -> Dict[str, np.ndarray]:
        """Per-hit arrays for a single event, under legacy ``hit_*`` names."""
        frame = self._frames[index]
        hits = frame.get(self._hits_coll)
        m = len(hits)
        slab = np.empty(m, np.int32)
        x = np.empty(m, np.float32)
        y = np.empty(m, np.float32)
        z = np.empty(m, np.float32)
        energy = np.empty(m, np.float32)
        for j, hit in enumerate(hits):
            pos = hit.getPosition()
            x[j], y[j], z[j] = pos.x, pos.y, pos.z
            energy[j] = hit.getEnergy()
            slab[j] = hit.getType()   # layer stored natively in CalorimeterHit.type
        out = {"hit_slab": slab, "hit_x": x, "hit_y": y, "hit_z": z, "hit_energy": energy}
        # Parallel per-hit UserDataCollections (same order as the hits).
        for field, (coll_name, dtype) in _USERDATA_HITS.items():
            coll = frame.get(coll_name)
            out[field] = np.fromiter(coll, dtype=dtype, count=len(coll))
        return out

    def all_hits(self) -> Dict[str, np.ndarray]:
        """All events' per-hit arrays as jagged object arrays (cached)."""
        if self._hits is None:
            n = len(self._frames)
            store = {f: np.empty(n, dtype=object) for f in PERHIT_FIELDS}
            for i in range(n):
                hits = self.read_hits(i)
                for f in PERHIT_FIELDS:
                    store[f][i] = hits[f]
            self._hits = store
        return self._hits

    def close(self) -> None:
        # podio Reader has no explicit close; drop references.
        self._frames = None
