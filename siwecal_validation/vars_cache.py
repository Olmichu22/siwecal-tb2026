"""
On-disk cache of the per-event derived variables (an augmented TTree).

The validation recomputes every per-event metric (Molière, shower shape,
barycenter, ...) on each run, which dominates the time when only the *cuts*
change. To avoid that, the first time a file is processed we write, next to it, a
``<stem>.valcache.root`` containing the **full original tree plus all derived
variables**; later runs read it back (vectorised, no recomputation) and only
apply cuts.

The cache holds **all valid events** (those with hits and positive energy, i.e.
exactly the :class:`EventData` set) with **no selection cuts** applied — cuts only
affect the plots. The same writer, with a ``row_mask``, also backs ``--save-tree``
(an augmented tree of the cut-passing events, written to the output directory).

Staleness: the derived variables depend on a few :class:`PlotConfig` parameters
(tungsten map, shower thresholds, ...). A fingerprint of those is stored in the
file; if it no longer matches, the cache is considered invalid and regenerated.

MIP-cut thresholds: the cache also stores pre-computed scalar metrics for
``hit_energy >= 0.5`` and ``hit_energy >= 1.0`` as extra branches with prefixes
``mip05_`` and ``mip1_``.  The event viewer reads these instead of recomputing
on the fly, making the interactive slider instantaneous.

Implementation note: writing uses PyROOT ``CloneTree`` so the original tree
(including its variable-length per-hit branches and their ``nhit_chan`` counter)
is reproduced faithfully — uproot's automatic counter naming clashes with the
existing ``nhit_slab``/``nhit_chip`` branches. Reading the derived branches back
uses uproot (vectorised ``arrays(..., library="np")``), which is the fast path.
"""

import json
import os
from dataclasses import fields

import numpy as np
import ROOT
import uproot

# Bump to invalidate every existing cache (e.g. when the schema changes).
CACHE_FORMAT_VERSION = 2
FINGERPRINT_KEY = "valcache_fingerprint"

# Pre-computed MIP-cut threshold levels stored in the cache.
MIP_CUT_THRESHOLDS = [0.5, 1.0]

# EventData fields that are 2-D per-layer arrays; not duplicated per threshold.
_PER_LAYER_FIELDS = {"hits_per_layer", "energy_per_layer", "weighte_per_layer"}


def _threshold_prefix(threshold: float) -> str:
    """Branch-name-safe prefix for a given hit_energy threshold.

    0.5 → ``"mip05"``, 1.0 → ``"mip1"``.
    """
    s = f"{threshold:.1f}".replace(".", "").rstrip("0")
    return f"mip{s}"


def _derived_names() -> list:
    """Names of the per-event variables stored in the cache (all but ``label``)."""
    from .event_data import EventData          # lazy: avoid import cycle
    return [f.name for f in fields(EventData) if f.name != "label"]


def _scalar_derived_names() -> list:
    """Scalar (1-D) derived field names — the ones replicated per MIP threshold."""
    return [n for n in _derived_names() if n not in _PER_LAYER_FIELDS]


def cache_path_for(input_path: str, cache_dir: str = None) -> str:
    """Path of the ``<stem>.valcache.root`` cache for ``input_path``.

    By default the cache sits next to the input file. Pass ``cache_dir`` to
    redirect every cache to a single directory instead -- useful when the data
    directory is read-only (the directory is created if needed).
    """
    stem = os.path.splitext(os.path.basename(input_path))[0]
    directory = cache_dir if cache_dir else os.path.dirname(input_path)
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)
    return os.path.join(directory, f"{stem}.valcache.root")


def config_fingerprint(config) -> str:
    """Stable JSON of the config parameters that affect the derived variables."""
    payload = {
        "version": CACHE_FORMAT_VERSION,
        "w_thicknesses": [float(w) for w in config.w_thicknesses],
        "w_x0_mm": float(config.w_x0_mm),
        "shower_profile": config.shower_profile,
        "shower_e_threshold": float(config.shower_e_threshold),
        "shower_max_min": float(config.shower_max_min),
        "shower_start_frac": float(config.shower_start_frac),
        "moliere_containment": float(config.moliere_containment),
        "n_layers": int(config.n_layers),
        "mip_cut_thresholds": MIP_CUT_THRESHOLDS,
    }
    return json.dumps(payload, sort_keys=True)


def is_valid(cache_path: str, config) -> bool:
    """True if the cache exists and its fingerprint matches the current config."""
    if not os.path.exists(cache_path):
        return False
    handle = ROOT.TFile(cache_path, "READ")
    if handle.IsZombie():
        return False
    stored = handle.Get(FINGERPRINT_KEY)
    matches = bool(stored) and stored.GetTitle() == config_fingerprint(config)
    handle.Close()
    return matches


def write(input_path: str, out_path: str, config, event_data,
          row_mask=None, extra_scalars=None) -> None:
    """Write an augmented tree (original branches + derived variables).

    Clones the original ``config.tree_name`` tree and appends one branch per
    :class:`EventData` field (scalars as ``/D``, ``is_shower`` as ``/O``, the
    per-layer profiles as ``[n_layers]/D`` arrays).  Additionally appends
    prefixed scalar branches for each threshold in ``extra_scalars``
    (``{0.5: {name: array}, 1.0: {name: array}}``).

    ``event_data`` must hold the valid events of ``input_path`` in order;
    ``row_mask`` (over those rows) optionally keeps a subset — used by
    ``--save-tree``. Written atomically.
    """
    derived = _derived_names()
    scalar_derived = _scalar_derived_names()
    n_layers = config.n_layers
    row_mask = None if row_mask is None else np.asarray(row_mask)
    extra_scalars = extra_scalars or {}

    fin = ROOT.TFile(input_path, "READ")
    in_tree = fin.Get(config.tree_name)
    if not in_tree:
        fin.Close()
        raise RuntimeError(f"Tree '{config.tree_name}' not found in {input_path}")

    clash = set(derived) & {b.GetName() for b in in_tree.GetListOfBranches()}
    if clash:
        fin.Close()
        raise ValueError(f"derived branch name(s) collide with original tree "
                         f"branches: {sorted(clash)}")

    per_layer = {n for n in derived if getattr(event_data, n).ndim == 2}

    tmp_path = out_path + ".tmp"
    fout = ROOT.TFile(tmp_path, "RECREATE")
    out_tree = in_tree.CloneTree(0)              # empty clone, shares in_tree buffers

    # Derived-branch buffers (kept alive for the whole fill loop).
    buffers = {}
    for name in derived:
        if name == "is_shower":
            buf = np.zeros(1, dtype=np.bool_)
            out_tree.Branch(name, buf, f"{name}/O")
        elif name in per_layer:
            buf = np.zeros(n_layers, dtype=np.float64)
            out_tree.Branch(name, buf, f"{name}[{n_layers}]/D")
        else:
            buf = np.zeros(1, dtype=np.float64)
            out_tree.Branch(name, buf, f"{name}/D")
        buffers[name] = buf

    # Extra MIP-threshold scalar branches.
    extra_buffers = {}   # (threshold, name) -> numpy buffer
    for thr, scalar_dict in extra_scalars.items():
        prefix = _threshold_prefix(thr)
        for name in scalar_derived:
            full = f"{prefix}_{name}"
            if name == "is_shower":
                buf = np.zeros(1, dtype=np.bool_)
                out_tree.Branch(full, buf, f"{full}/O")
            else:
                buf = np.zeros(1, dtype=np.float64)
                out_tree.Branch(full, buf, f"{full}/D")
            extra_buffers[(thr, name)] = buf

    valid_index = 0
    for entry in range(in_tree.GetEntries()):
        in_tree.GetEntry(entry)
        n_channels = int(in_tree.nhit_chan)
        if n_channels == 0:
            continue
        energy = np.frombuffer(in_tree.hit_energy, np.float32, count=n_channels)
        if energy.sum() <= 0:
            continue
        if row_mask is None or row_mask[valid_index]:
            for name in derived:
                value = getattr(event_data, name)[valid_index]
                if name in per_layer:
                    buffers[name][:] = value
                else:
                    buffers[name][0] = value
            for (thr, name), buf in extra_buffers.items():
                buf[0] = extra_scalars[thr][name][valid_index]
            out_tree.Fill()                      # copies the original branches too
        valid_index += 1
    fin.Close()

    if valid_index != len(event_data):
        fout.Close()
        os.remove(tmp_path)
        raise ValueError(
            f"valid-event count mismatch for {input_path}: tree {valid_index} "
            f"vs EventData {len(event_data)}; the input may have changed.")

    fout.cd()
    out_tree.Write()
    ROOT.TNamed(FINGERPRINT_KEY, config_fingerprint(config)).Write()
    fout.Close()
    os.replace(tmp_path, out_path)


def read(cache_path: str, label: str, config):
    """Reconstruct :class:`EventData` from the cache, reading only derived branches."""
    from .event_data import EventData          # lazy: avoid import cycle
    derived = _derived_names()
    with uproot.open(cache_path) as handle:
        arrays = handle[config.tree_name].arrays(derived, library="np")
    fields_dict = {name: np.asarray(arrays[name]) for name in derived}
    return EventData(label=label, **fields_dict)
