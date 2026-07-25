#!/usr/bin/env python
"""
Classify every eudaq run's decode route (raw vs LCIO) for
generate_convert_eudaq_dag.py's --route-map, by actually test-decoding one
raw chunk per run and checking entries -- NOT by guessing from the run name
or which raw filename pattern matched. Neither is reliable: verified against
real data (2026-07-20), TB2026CERN_eudaq_run_000316 uses plain _raw.bin data
and decodes fine via the raw path, while _000183 and _000146 use the bare
".bin" naming and decode to 0 entries via the raw path EVEN WITH
RAW2ROOT_EUDAQ=1 -- they only work via LCIO (gaudi_jobs/decode_lcio_runs.py).
A naming-based guess broke 159 of 161 runs the first time this was tried.

For each run: test-decode its first raw chunk with the filename-implied
EudaqFormat flag (bare "<run>.bin" -> True, "<run>_raw.bin" -> False) into a
throwaway scratch file and check entries via ROOT.
  entries > 0  -> route "raw", recording the eudaq_format flag that worked.
  entries == 0 -> route "lcio" if Data/eudaq/ROC_run_<N>_tp.slcio exists,
                  else route "failed" (needs manual attention).

Needs key4hep sourced (k4run on PATH) before running -- this is NOT a
lightweight generator like generate_convert_eudaq_dag.py; it actually runs
Gaudi. Runs are classified in parallel (--workers) since each is dominated by
k4run's fixed startup overhead, not the (small, single-chunk) decode itself.

Usage::

    source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "$(cat .key4hep-release)"
    export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
    export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH

    python gaudi_jobs/condor/classify_eudaq_routes.py \\
        --out gaudi_jobs/condor/generated/convert_eudaq/route_map.json
"""
import argparse
import concurrent.futures
import json
import os
import subprocess
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_REPO, "gaudi_jobs", "calibration"))
sys.path.insert(0, os.path.join(_REPO, "gaudi_jobs"))
sys.path.insert(0, _REPO)

from calib_run_utils import DEFAULT_RAW_BASE, raw_chunk_files  # noqa: E402
from decode_lcio_runs import DEFAULT_EUDAQ_BASE, lcio_path_for  # noqa: E402

_STEERING = os.path.join(_REPO, "gaudi_source", "options", "run_raw2root.py")

_ENTRIES_CHECK = (
    "import ROOT,sys; ROOT.gErrorIgnoreLevel=ROOT.kFatal; "
    "f=ROOT.TFile.Open(sys.argv[1]); t=f.Get('siwecaldecoded') if f and not f.IsZombie() else None; "
    "print(t.GetEntries() if t else -1)"
)


def _entries_of(path):
    r = subprocess.run(["python3", "-c", _ENTRIES_CHECK, path], capture_output=True, text=True, timeout=60)
    try:
        return int(r.stdout.strip())
    except ValueError:
        return -1


def _eudaq_run_number(run):
    return int("".join(ch for ch in run.split("_")[-1] if ch.isdigit()))


def classify(run, raw_base, eudaq_base, scratch_dir):
    raw_dir = os.path.join(raw_base, run)
    try:
        chunks = raw_chunk_files(raw_dir, run)
    except SystemExit as e:
        return run, {"route": "failed", "reason": str(e)}

    first = chunks[0]
    eudaq_format = not os.path.basename(first).startswith(f"{run}_raw.bin")
    out = os.path.join(scratch_dir, f"{run}.root")
    run_settings = os.path.join(raw_dir, "Run_Settings.txt")
    env = {**os.environ, "RAW_FILES": first, "RAW2ROOT_OUT": out, "RAW2ROOT_RUN_SETTINGS_FILE": run_settings}
    if eudaq_format:
        env["RAW2ROOT_EUDAQ"] = "1"
    try:
        r = subprocess.run(["k4run", _STEERING], env=env, capture_output=True, text=True, timeout=180)
    except subprocess.TimeoutExpired:
        return run, {"route": "failed", "reason": "raw test-decode timed out"}
    if r.returncode != 0:
        return run, {"route": "failed", "reason": f"k4run exit {r.returncode}: {r.stderr[-300:]}"}

    n = _entries_of(out)
    try:
        os.remove(out)
    except OSError:
        pass
    if n > 0:
        return run, {"route": "raw", "eudaq_format": eudaq_format, "n_chunks": len(chunks), "test_entries": n}

    lcio_in = lcio_path_for(_eudaq_run_number(run), eudaq_base)
    if os.path.exists(lcio_in):
        return run, {"route": "lcio", "lcio_in": lcio_in, "n_chunks": len(chunks)}
    return run, {"route": "failed", "reason": "raw gives 0 entries and no LCIO file", "n_chunks": len(chunks)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", default=None,
                   help="Comma-separated run names to classify. Default: every 'TB2026CERN_eudaq_run_*' "
                        "folder under --raw-base.")
    p.add_argument("--raw-base", default=DEFAULT_RAW_BASE, help=f"Default: {DEFAULT_RAW_BASE}")
    p.add_argument("--eudaq-base", default=DEFAULT_EUDAQ_BASE, help=f"Default: {DEFAULT_EUDAQ_BASE}")
    p.add_argument("--out", required=True, help="Where to write the route_map.json.")
    p.add_argument("--workers", type=int, default=10,
                   help="Parallel k4run test-decodes (default 10). Each is dominated by k4run's fixed "
                        "startup overhead (~10-15s), not the actual (single-chunk) decode.")
    args = p.parse_args(argv)

    if args.runs:
        runs = [r.strip() for r in args.runs.split(",") if r.strip()]
    else:
        runs = sorted(
            name for name in os.listdir(args.raw_base)
            if "eudaq" in name and os.path.isdir(os.path.join(args.raw_base, name))
        )
    print(f"{len(runs)} run(s) to classify", file=sys.stderr)

    with tempfile.TemporaryDirectory(prefix="classify_eudaq_") as scratch_dir:
        result = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
            futures = {ex.submit(classify, run, args.raw_base, args.eudaq_base, scratch_dir): run for run in runs}
            for i, fut in enumerate(concurrent.futures.as_completed(futures), 1):
                run, info = fut.result()
                result[run] = info
                print(f"[{i}/{len(runs)}] {run}: {info}", file=sys.stderr)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)

    n_raw = sum(1 for v in result.values() if v["route"] == "raw")
    n_lcio = sum(1 for v in result.values() if v["route"] == "lcio")
    n_failed = sum(1 for v in result.values() if v["route"] == "failed")
    print(f"\nraw={n_raw} lcio={n_lcio} failed={n_failed}  -> {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
