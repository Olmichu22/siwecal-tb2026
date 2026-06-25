#!/usr/bin/env python
"""
Batch driver for the Gaudi PID stage: resolve inputs exactly like
``siwecal_validation`` (--run / --file / --all / --point / --cfg) and run
``options/run_pid.py`` once per ``ecal_<run>.root``, writing one EDM4hep file each.

Run after sourcing key4hep and exporting LD_LIBRARY_PATH / PYTHONPATH (see the
top-level README). Example::

    python k4SiWEcalReco/run_pid_batch.py --all --outdir /tmp/pid
    python k4SiWEcalReco/run_pid_batch.py --run TB2026CERN_run_000007
"""
import argparse
import os
import subprocess
import sys

import yaml

from siwecal_common import paths
from siwecal_validation.config import BASE_PATH, DEFAULT_CONFIG

_STEERING = os.path.join(os.path.dirname(__file__), "options", "run_pid.py")


def _label(events_path):
    stem = os.path.splitext(os.path.basename(events_path))[0]
    return stem[len("ecal_"):] if stem.startswith("ecal_") else stem


def _resolve(events_path, base_path):
    if os.path.isabs(events_path):
        return events_path
    candidate = os.path.join(base_path, events_path)
    return candidate if os.path.exists(candidate) else paths.resolve_input(events_path)


def _jobs(args):
    """Yield (label, ecal_path) pairs from the CLI selection."""
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
            yield _label(path), path
    elif args.file:
        yield (args.run or _label(args.file)), args.file
    else:
        run = args.run or "TB2026CERN_run_000007"
        yield run, paths.resolve_input(os.path.join(run, f"ecal_{run}.root"))


def main(argv=None):
    p = argparse.ArgumentParser(description="Batch Gaudi PID stage (ecal -> EDM4hep).")
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
    args = p.parse_args(argv)

    out_base = args.outdir or paths.pid_dir()   # None -> next to each input
    if out_base:
        os.makedirs(out_base, exist_ok=True)
    failures = 0
    for label, ecal_path in _jobs(args):
        if not os.path.exists(ecal_path):
            print(f"SKIP {label}: not found ({ecal_path})")
            failures += 1
            continue
        out_dir = out_base or os.path.dirname(ecal_path)
        out = os.path.join(out_dir, f"ecal_pid_{label}.root")
        print(f"\n=== {label}: {ecal_path} -> {out} ===")
        env = {**os.environ, "ECAL_FILE": ecal_path, "ECAL_PID_OUT": out, "ECAL_TREE": args.tree}
        result = subprocess.run(["k4run", _STEERING], env=env)
        if result.returncode != 0:
            print(f"ERROR: k4run failed for {label}")
            failures += 1
    print(f"\n[Done] {failures} failure(s).")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
