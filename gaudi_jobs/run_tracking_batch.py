#!/usr/bin/env python
"""
Batch driver for the Gaudi ACTS tracking stage: resolve inputs the same way
``run_pid_batch.py`` does (--run / --file / --all) and run
``options/run_tracking.py`` once per ``ecal_<label>.edm4hep.root`` (the PID
batch driver's own output -- it keeps the ``ECalHits`` collection tracking
needs, see ``siwecal_common/edm4hep_pid.py::write_filtered``).

Writes ``tracks_<label>.edm4hep.root`` (``ACTSTracks`` + ``EMShowers``) next to
each input by default, mirroring ``run_pid_batch.py``'s ``pid_dir()``-based
output layout.

Run after ``run_pid_batch.py`` has produced the EDM4hep PID file(s), and after
sourcing key4hep / exporting LD_LIBRARY_PATH / PYTHONPATH (see the top-level
README). Example::

    python gaudi_jobs/run_pid_batch.py --run TB2026CERN_run_000013
    python gaudi_jobs/run_tracking_batch.py --run TB2026CERN_run_000013

    python gaudi_jobs/run_tracking_batch.py --all
"""
import argparse
import os
import subprocess
import sys

import yaml

from siwecal_common import paths
from siwecal_validation.config import BASE_PATH, DEFAULT_CONFIG

_STEERING = os.path.join(os.path.dirname(__file__), "..", "gaudi_source", "options", "run_tracking.py")


def _label(events_path):
    stem = os.path.splitext(os.path.splitext(os.path.basename(events_path))[0])[0]
    return stem[len("ecal_"):] if stem.startswith("ecal_") else stem


def _default_pid_path(run, base_dir):
    return os.path.join(base_dir, f"ecal_{run}.edm4hep.root")


def _jobs(args):
    """Yield (label, ecal_pid_path) pairs from the CLI selection."""
    pid_base = paths.pid_dir() or paths.reconstruction_dir()
    if args.all:
        with open(args.cfg) as handle:
            cfg = yaml.safe_load(handle) or {}
        data_map = cfg.get("event_data", {})
        if not data_map:
            raise SystemExit(f"ERROR: no matching 'event_data' entries in {args.cfg}")
        for name in data_map:
            yield name, _default_pid_path(name, pid_base)
    elif args.file:
        yield (args.run or _label(args.file)), args.file
    else:
        run = args.run or "TB2026CERN_run_000007"
        yield run, _default_pid_path(run, pid_base)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Batch Gaudi ACTS tracking stage (ecal_<label>.edm4hep.root -> tracks_<label>.edm4hep.root).")
    p.add_argument("--run", default=None, help="Single run name (label of an existing ecal_<run>.edm4hep.root)")
    p.add_argument("--file", default=None, help="Explicit input EDM4hep file (must contain ECalHits)")
    p.add_argument("--cfg", default=DEFAULT_CONFIG, help="data_reference YAML, for --all")
    p.add_argument("--all", action="store_true", help="Every 'event_data' entry of the YAML")
    p.add_argument("--outdir", default=None,
                   help="Output directory. Default: settings.yml 'pid_dir', else "
                        "next to each input file.")
    p.add_argument("--input-collection", default="ECalHits",
                   help="Hit collection to track on (default: ECalHits)")
    p.add_argument("--seed-momentum", type=float, default=3.0, help="Beam momentum [GeV] (default: 3.0)")
    p.add_argument("--debug", action="store_true", help="TRACKING_LOGLEVEL=DEBUG")
    args = p.parse_args(argv)

    out_base = args.outdir or paths.pid_dir()   # None -> next to each input
    if out_base:
        os.makedirs(out_base, exist_ok=True)

    failures = 0
    for label, in_path in _jobs(args):
        if not os.path.exists(in_path):
            print(f"SKIP {label}: not found ({in_path}) -- run run_pid_batch.py first?")
            failures += 1
            continue
        out_dir = out_base or os.path.dirname(in_path)
        out_path = os.path.join(out_dir, f"tracks_{label}.edm4hep.root")

        print(f"\n=== {label}: {in_path} -> {out_path} ===")
        env = {**os.environ,
               "TRACKING_INPUT_FILE": in_path,
               "TRACKING_OUTPUT_FILE": out_path,
               "TRACKING_INPUT_COLLECTION": args.input_collection,
               "SEED_MOMENTUM": str(args.seed_momentum),
               "TRACKING_LOGLEVEL": "DEBUG" if args.debug else "INFO"}
        result = subprocess.run(["k4run", _STEERING], env=env)
        if result.returncode != 0:
            print(f"ERROR: k4run failed for {label}")
            failures += 1
            continue
        print(f"[Output] -> {out_path}")

    print(f"\n[Done] {failures} failure(s).")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
