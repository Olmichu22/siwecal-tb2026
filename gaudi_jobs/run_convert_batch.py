#!/usr/bin/env python
"""
Convert a range or list of runs to `siwecaldecoded` ROOT files (raw2root ONLY).

One `k4run EcalRawDecoder` process per RAW CHUNK, decoded into
`<converted-dir>/<run>/chunks/chunk_NNNN.root`. No event building, no PID --
just the decoded trees. This is the convert-only counterpart of
`run_full_pipeline_batch.py`, for when you just want several runs decoded.

It used to decode a whole run in one process, into a single `<run>.root`. That
silently dropped ~75% of the acquisitions on some runs; see gaudi_jobs/decode_chunks.py
for the measurements. Every decode now ends with a health check that fails loudly
if the acquisitions do not add up, so a truncated run can no longer be quietly
reconstructed.

The chunks ARE the decoded data -- not an intermediate to merge away. The
calibration Fill reads them and the event builder chains them.

For fanning a run out across the batch farm instead of decoding locally, use
gaudi_jobs/condor/generate_reco_dag.py (or calibration/condor/generate_convert_jobs.py).

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
import sys

# calib_run_utils lives in gaudi_jobs/calibration/; reuse its raw-chunk
# discovery (handles the plain vs eudaq `_raw.bin`/`.bin` naming quirk) and the
# shared read-only raw base, rather than re-deriving them here.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibration"))
from calib_run_utils import DEFAULT_RAW_BASE  # noqa: E402
from decode_chunks import assert_healthy, assert_not_under, chunks_dir, decode_run  # noqa: E402

_DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"

_RANGE_RE = re.compile(r"^(\d+)-(\d+)$")


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
                   help=f"Output base directory; each run writes <converted-dir>/<run>/chunks/chunk_NNNN.root. "
                        f"Default: {_DEFAULT_CONVERTED_DIR}")
    p.add_argument("--tree", default="siwecaldecoded", help="Output TTree name (default siwecaldecoded)")
    p.add_argument("--max-cycle-jump", type=int, default=10, help="EcalRawDecoder MaxReadoutCycleJump (default 10)")
    p.add_argument("--force", action="store_true",
                   help="Re-decode chunks whose output already exists (default: skip them, so an "
                        "interrupted decode is resumed by simply re-running)")
    p.add_argument("--dry-run", action="store_true", help="Print what would be decoded without running k4run")
    args = p.parse_args(argv)

    run_folders = _parse_runs(args.runs, args.raw_base)
    assert_not_under(args.converted_dir, args.raw_base)

    print(f"[convert] {len(run_folders)} run(s) -> {args.converted_dir}/<run>/chunks/")
    failures = []
    for run_name, raw_dir in run_folders:
        try:
            chunks = decode_run(run_name, raw_dir, args.converted_dir,
                                force=args.force, dry_run=args.dry_run,
                                tree=args.tree, max_cycle_jump=args.max_cycle_jump)
            if not args.dry_run:
                assert_healthy(chunks, run_name, tree=args.tree)
        except SystemExit as exc:
            print(f"[FAIL] {run_name}: {exc}")
            failures.append(run_name)
            continue

    if failures:
        print(f"\n[Done with errors] {len(failures)} run(s) failed: {', '.join(failures)}")
        return 1
    print(f"\n[Done] {len(run_folders)} run(s) decoded to {chunks_dir(args.converted_dir, '<run>')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
