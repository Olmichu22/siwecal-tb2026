#!/usr/bin/env python
"""
Convert a range or list of runs to `siwecaldecoded` ROOT files (raw2root ONLY).

One `k4run EcalRawDecoder` process per run decodes all of that run's raw chunk
files into a single `<converted-dir>/<run>/<run>.root`. No event building, no
PID -- just the decoded tree. This is the convert-only counterpart of
`run_full_pipeline_batch.py` (which additionally does event building + PID for
one run), for when you just want the decoded trees of several runs at once.

Since the decoder streams acquisitions to disk (constant memory, see
EcalRawDecoder.cpp), decoding a whole run in one process is fine regardless of
how many chunks it has. For fanning a single huge run out across the batch farm
per-chunk instead, use gaudi_jobs/calibration/condor/generate_convert_jobs.py.

SAFETY: the raw data directory is READ-ONLY -- this script only reads chunks
from it and writes decoded output to --converted-dir; it never deletes anything.

Usage::

    source setup.sh
    export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
    export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH

    # a contiguous run range (inclusive)
    python gaudi_jobs/run_convert_batch.py --runs 12-20

    # ranges, bare names and explicit folder paths can be mixed
    python gaudi_jobs/run_convert_batch.py --runs 12-20,60,TB2026CERN_eudaq_run_000146

    # see what would run without decoding anything
    python gaudi_jobs/run_convert_batch.py --runs 12-20 --dry-run
"""
import argparse
import os
import re
import subprocess
import sys

# calib_run_utils lives in gaudi_jobs/calibration/; reuse its raw-chunk
# discovery (handles the plain vs eudaq `_raw.bin`/`.bin` naming quirk) and the
# shared read-only raw base, rather than re-deriving them here.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration"))
from calib_run_utils import DEFAULT_RAW_BASE, raw_chunk_files  # noqa: E402

_OPTIONS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "gaudi_source", "options")
_RAW2ROOT_STEERING = os.path.join(_OPTIONS_DIR, "run_raw2root.py")

_DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"

_RANGE_RE = re.compile(r"^(\d+)-(\d+)$")


def _assert_not_under(path, readonly_root):
    real_path = os.path.realpath(path)
    real_root = os.path.realpath(readonly_root)
    if real_path == real_root or real_path.startswith(real_root + os.sep):
        raise SystemExit(f"REFUSING to write inside the read-only raw data directory: {path} is under {readonly_root}")


def _resolve_run(token, raw_base):
    """Resolve one --runs token to ``(run_name, raw_dir)``.

    A token is a full folder path (contains a separator), a bare run name, or a
    plain integer run number -- which is resolved to whichever of the plain
    ``TB2026CERN_run_<N>`` or ``TB2026CERN_eudaq_run_<N>`` folder actually
    exists on disk (the prefix is just the trigger mode, not a different run).
    """
    if os.sep in token:
        raw_dir = token.rstrip(os.sep)
        return os.path.basename(raw_dir), raw_dir
    if token.isdigit():
        n = int(token)
        for tmpl in ("TB2026CERN_run_%06d", "TB2026CERN_eudaq_run_%06d"):
            name = tmpl % n
            raw_dir = os.path.join(raw_base, name)
            if os.path.isdir(raw_dir):
                return name, raw_dir
        raise SystemExit(f"ERROR: no raw folder for run {n} under {raw_base} "
                         f"(tried TB2026CERN_run_{n:06d} and TB2026CERN_eudaq_run_{n:06d})")
    return token, os.path.join(raw_base, token)


def _parse_runs(spec, raw_base):
    """Expand a ``--runs`` spec into an ordered ``[(run_name, raw_dir), ...]``.

    Comma-separated tokens, each a contiguous integer range (``12-20``,
    inclusive), a bare run name, an integer, or a full folder path.
    """
    resolved = []
    for token in (t.strip() for t in spec.split(",")):
        if not token:
            continue
        m = _RANGE_RE.match(token)
        if m:
            lo, hi = int(m.group(1)), int(m.group(2))
            if hi < lo:
                raise SystemExit(f"ERROR: inverted run range '{token}'")
            for n in range(lo, hi + 1):
                resolved.append(_resolve_run(str(n), raw_base))
        else:
            resolved.append(_resolve_run(token, raw_base))
    if not resolved:
        raise SystemExit("ERROR: --runs is empty")
    return resolved


def main(argv=None):
    p = argparse.ArgumentParser(description="Convert a range/list of runs to siwecaldecoded ROOT (raw2root only).")
    p.add_argument("--runs", required=True,
                   help="Comma-separated runs: integer ranges (e.g. 12-20), bare run names, run numbers, or full "
                        "folder paths, freely mixed.")
    p.add_argument("--raw-base", default=DEFAULT_RAW_BASE,
                   help=f"Base directory to resolve run numbers/names against (READ-ONLY). Default: {DEFAULT_RAW_BASE}")
    p.add_argument("--converted-dir", default=_DEFAULT_CONVERTED_DIR,
                   help=f"Output base directory; each run writes <converted-dir>/<run>/<run>.root. "
                        f"Default: {_DEFAULT_CONVERTED_DIR}")
    p.add_argument("--tree", default="siwecaldecoded", help="Output TTree name (default siwecaldecoded)")
    p.add_argument("--max-cycle-jump", type=int, default=10, help="EcalRawDecoder MaxReadoutCycleJump (default 10)")
    p.add_argument("--force", action="store_true",
                   help="Re-decode even if the output <run>.root already exists (default: skip existing)")
    p.add_argument("--dry-run", action="store_true", help="Print what would be decoded without running k4run")
    args = p.parse_args(argv)

    run_folders = _parse_runs(args.runs, args.raw_base)
    _assert_not_under(args.converted_dir, args.raw_base)

    print(f"[convert] {len(run_folders)} run(s) -> {args.converted_dir}/<run>/<run>.root")
    failures = []
    for run_name, raw_dir in run_folders:
        out_dir = os.path.join(args.converted_dir, run_name)
        out_file = os.path.join(out_dir, f"{run_name}.root")
        try:
            chunks = raw_chunk_files(raw_dir, run_name)
        except SystemExit as exc:
            print(f"[skip] {run_name}: {exc}")
            failures.append(run_name)
            continue

        if os.path.exists(out_file) and not args.force:
            print(f"[skip] {run_name}: {out_file} already exists (use --force to re-decode)")
            continue

        print(f"[convert] {run_name}: {len(chunks)} chunk(s) -> {out_file}")
        if args.dry_run:
            continue

        os.makedirs(out_dir, exist_ok=True)
        # Pass the chunk list via a file (RAW_FILES_LIST), not RAW_FILES: a run
        # with hundreds/thousands of chunks would blow past the execve argument
        # size limit as one comma-joined env var.
        chunklist = os.path.join(out_dir, "chunklist.txt")
        with open(chunklist, "w") as fh:
            fh.write("\n".join(chunks) + "\n")

        env = {
            **os.environ,
            "RAW_FILES_LIST": chunklist,
            "RAW2ROOT_OUT": out_file,
            "RAW2ROOT_TREE": args.tree,
            "RAW2ROOT_MAX_CYCLE_JUMP": str(args.max_cycle_jump),
            "RAW2ROOT_RUN_SETTINGS_FILE": os.path.join(raw_dir, "Run_Settings.txt"),
        }
        result = subprocess.run(["k4run", _RAW2ROOT_STEERING], env=env)
        if result.returncode != 0:
            print(f"[FAIL] {run_name}: k4run raw2root exited {result.returncode}")
            failures.append(run_name)

    if failures:
        print(f"\n[Done with errors] {len(failures)} run(s) failed: {', '.join(failures)}")
        return 1
    print(f"\n[Done] {len(run_folders)} run(s) processed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
