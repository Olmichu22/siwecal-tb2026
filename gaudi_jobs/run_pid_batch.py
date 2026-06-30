#!/usr/bin/env python
"""
Batch driver for the Gaudi PID stage: resolve inputs exactly like
``siwecal_validation`` (--run / --file / --all / --point / --cfg) and run
``options/run_pid.py`` once per ``ecal_<run>.root``.

For each input the Gaudi stage writes a full EDM4hep PID file to a temporary
path; this driver then applies the event-selection cuts (the same ``CutSet`` as
``siwecal_validation``) and writes the requested final, cut-passing output(s):

* ``ecal_<label>.edm4hep.root`` -- filtered EDM4hep PID file.
* ``ecal_<label>.valtree.root`` -- a plain valcache-schema TTree.

Run modes (``--validation``):

* default (physics): a hit-level MIP cut ``>= 0.5`` is baked in per cell and the
  ``mip05_/mip1_`` variant blocks are not computed.
* ``--validation`` (visualizer): no hit cut, and the MIP-variant blocks are
  computed in the Cluster (current behavior, for the online event viewer).

Cuts: every ``CutSet`` variable is a ``--<var>-min/--<var>-max`` flag, applied to
all runs; in ``--all``/``--point`` a per-energy ``cuts:`` block in the YAML
overrides them. All cuts default to off except total per-event energy > 0, which
is always enforced (events with no hits or non-positive energy are dropped).

Run after sourcing key4hep and exporting LD_LIBRARY_PATH / PYTHONPATH (see the
top-level README). Example::

    python gaudi_jobs/run_pid_batch.py --all --format both
    python gaudi_jobs/run_pid_batch.py --run TB2026CERN_run_000007 --validation
"""
import argparse
import os
import subprocess
import sys

import yaml

from siwecal_common import paths
from siwecal_common import edm4hep_pid
from siwecal_validation.config import BASE_PATH, DEFAULT_CONFIG, PlotConfig
from siwecal_validation.cli import cutset_from_args
from siwecal_validation.event_data import EventData
from siwecal_validation.selection import CutSet, _CUT_SPEC

_STEERING = os.path.join(os.path.dirname(__file__), "..", "gaudi_source", "options", "run_pid.py")


def _label(events_path):
    stem = os.path.splitext(os.path.basename(events_path))[0]
    return stem[len("ecal_"):] if stem.startswith("ecal_") else stem


def _resolve(events_path, base_path):
    if os.path.isabs(events_path):
        return events_path
    candidate = os.path.join(base_path, events_path)
    return candidate if os.path.exists(candidate) else paths.resolve_input(events_path)


def _jobs(args):
    """Yield (label, ecal_path, entry) triples from the CLI selection.

    ``entry`` is the YAML mapping for the sample (with an optional per-energy
    ``cuts:`` block and ``hit_mip_cut``), or ``{}`` for --run/--file.
    """
    if args.all or args.point is not None:
        with open(args.cfg) as handle:
            cfg = yaml.safe_load(handle) or {}
        data_map = cfg.get("event_data", {})
        if args.point is not None:
            data_map = {k: v for k, v in data_map.items() if k.startswith(f"P{args.point}_")}
        if not data_map:
            raise SystemExit(f"ERROR: no matching 'event_data' entries in {args.cfg}")
        base_path = cfg.get("main_path", BASE_PATH)
        for entry in data_map.values():
            path = _resolve(entry["path"], base_path)
            yield _label(path), path, entry
    elif args.file:
        yield (args.run or _label(args.file)), args.file, {}
    else:
        run = args.run or "TB2026CERN_run_000007"
        yield run, paths.resolve_input(os.path.join(run, f"ecal_{run}.root")), {}


def _post_process(temp, out_dir, label, effective_cut, mip_thresholds, fmt):
    """Apply ``effective_cut`` to the full EDM4hep ``temp`` file and write the
    requested final output(s). Returns the number of cut-passing events."""
    config = PlotConfig()
    mip = tuple(mip_thresholds)
    data = EventData.from_edm4hep(temp, label, config, mip_thresholds=mip)
    keep = effective_cut.mask(data)
    frame_indices = data.source_index[keep]
    n_kept = int(keep.sum())

    if fmt in ("edm4hep", "both"):
        out = os.path.join(out_dir, f"ecal_{label}.edm4hep.root")
        edm4hep_pid.write_filtered(temp, out, frame_indices)
        print(f"[Output] {n_kept} event(s) -> {out}")
    if fmt in ("valtree", "both"):
        out = os.path.join(out_dir, f"ecal_{label}.valtree.root")
        reader = edm4hep_pid.PidFileReader(temp, n_layers=config.n_layers,
                                           mip_thresholds=mip)
        edm4hep_pid.write_valtree(reader, out, config, frame_indices,
                                  mip_thresholds=mip)
        print(f"[Output] {n_kept} event(s) -> {out}")
    return n_kept


def main(argv=None):
    p = argparse.ArgumentParser(description="Batch Gaudi PID stage (ecal -> EDM4hep / valtree).")
    p.add_argument("--run", default=None, help="Single run name")
    p.add_argument("--file", default=None, help="Explicit ecal_<run>.root path")
    p.add_argument("--cfg", default=DEFAULT_CONFIG, help="data_reference YAML")
    bulk = p.add_mutually_exclusive_group()
    bulk.add_argument("--all", action="store_true", help="Every 'event_data' entry of the YAML")
    bulk.add_argument("--point", type=int, choices=range(1, 6), help="Only beam point P<N>_*")
    p.add_argument("--outdir", default=None,
                   help="Output directory. Default: settings.yml 'pid_dir', else "
                        "next to each input ecal file.")
    p.add_argument("--tree", default="ecal", help="Input TTree name")
    p.add_argument("--format", choices=("edm4hep", "valtree", "both"), default="edm4hep",
                   help="Output format(s) to write (default: edm4hep).")
    p.add_argument("--validation", action="store_true",
                   help="Visualizer mode: compute the mip05_/mip1_ variant blocks "
                        "and disable the hit-level MIP cut (current behavior). "
                        "Without it the hit cut defaults to >= 0.5 MIP and the "
                        "variant blocks are omitted.")
    p.add_argument("--hit-mip-cut", type=float, default=None, dest="hit_mip_cut",
                   help="Drop hits below this MIP energy (override; default: 0.5 "
                        "in physics mode, off with --validation).")
    # Event-selection cuts: one --<var>-min/--<var>-max per CutSet variable,
    # generated from the spec so the CLI never drifts from selection.py. All
    # default off; total energy > 0 is always enforced by the EDM4hep loader.
    cuts = p.add_argument_group(
        "selection cuts",
        "min/max bounds applied to all samples (per-energy YAML cuts override)")
    for field_min, field_max, _name, _attr in _CUT_SPEC:
        cuts.add_argument(f"--{field_min.replace('_', '-')}", type=float, default=None)
        cuts.add_argument(f"--{field_max.replace('_', '-')}", type=float, default=None)
    shower = cuts.add_mutually_exclusive_group()
    shower.add_argument("--is-shower", dest="is_shower", action="store_const",
                        const=True, default=None, help="keep only shower events")
    shower.add_argument("--no-shower", dest="is_shower", action="store_const",
                        const=False, help="keep only non-shower events")
    args = p.parse_args(argv)

    general_cut = cutset_from_args(args)
    mip_thresholds = [0.5, 1.0] if args.validation else []
    default_hit_cut = -1.0 if args.validation else 0.5

    out_base = args.outdir or paths.pid_dir()   # None -> next to each input
    if out_base:
        os.makedirs(out_base, exist_ok=True)
    failures = 0
    for label, ecal_path, entry in _jobs(args):
        if not os.path.exists(ecal_path):
            print(f"SKIP {label}: not found ({ecal_path})")
            failures += 1
            continue
        out_dir = out_base or os.path.dirname(ecal_path)
        effective_cut = general_cut.merge(CutSet.from_mapping(entry.get("cuts")))
        hit_cut = (args.hit_mip_cut if args.hit_mip_cut is not None
                   else float(entry.get("hit_mip_cut", default_hit_cut)))
        temp = os.path.join(out_dir, f".{label}.pid.tmp.root")

        print(f"\n=== {label}: {ecal_path} (hit_mip_cut={hit_cut}, "
              f"mip_blocks={mip_thresholds or 'none'}, format={args.format}) ===")
        env = {**os.environ, "ECAL_FILE": ecal_path, "ECAL_PID_OUT": temp,
               "ECAL_TREE": args.tree, "ECAL_HIT_MIP_CUT": str(hit_cut),
               "ECAL_MIP_THRESHOLDS": ",".join(str(t) for t in mip_thresholds)}
        result = subprocess.run(["k4run", _STEERING], env=env)
        if result.returncode != 0:
            print(f"ERROR: k4run failed for {label}")
            failures += 1
            continue
        try:
            _post_process(temp, out_dir, label, effective_cut, mip_thresholds, args.format)
        except (OSError, RuntimeError, ValueError) as error:
            print(f"ERROR: post-processing failed for {label}: {error}")
            failures += 1
        finally:
            if os.path.exists(temp):
                os.remove(temp)
    print(f"\n[Done] {failures} failure(s).")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
