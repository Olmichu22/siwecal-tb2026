#!/usr/bin/env python
"""
Generate a Condor/DAGMan pipeline that CONVERTs every eudaq run (raw ->
siwecaldecoded chunks), routing each run to the raw decoder (EcalRawDecoder)
or the LCIO decoder (EcalLcioDecoder, gaudi_jobs/decode_lcio_runs.py) per a
pre-computed --route-map.

Why a route map instead of a naming heuristic: neither "eudaq" in the run
name nor which raw filename pattern matched ("<run>_raw.bin" vs "<run>.bin")
reliably predicts which decoder actually produces a non-empty tree. Verified
against real data (2026-07-20): TB2026CERN_eudaq_run_000316 uses plain
_raw.bin data and decodes fine via the raw path; TB2026CERN_eudaq_run_000183
and _000146 use the bare ".bin" naming and decode to 0 entries via the raw
path EVEN WITH RAW2ROOT_EUDAQ=1 -- they only work via LCIO. A first attempt
at this script keyed routing off the run name and broke 159/161 runs. Build
the route map with the classification helper (test-decode one chunk per run,
check entries, fall back to LCIO if an .slcio file exists) before running
this generator -- see gaudi_jobs/decode_chunks.py for the same per-run
eudaq_format detection used on the serial/local decode path.

Convert-only counterpart of generate_reco_dag.py, restricted to eudaq runs.
Reuses the same grouped-CONVERT machinery as generate_reco_dag.py (one Condor
job decodes --chunks-per-job raw chunks, each in its own k4run) and the same
idempotent skip-if-already-decoded check (entries > 0, not just "tree
exists": production already holds 0-entry chunks from before this was
understood).

Only WRITES files -- never calls condor_submit itself.

Route map schema (JSON): {run_name: {"route": "raw", "eudaq_format": bool}}
or {"route": "lcio", "lcio_in": path} or {"route": "failed", "reason": str}.

Usage::

    python gaudi_jobs/condor/generate_convert_eudaq_dag.py \\
        --route-map gaudi_jobs/condor/generated/convert_eudaq/route_map.json \\
        --out-dir gaudi_jobs/condor/generated/convert_eudaq

    condor_submit_dag gaudi_jobs/condor/generated/convert_eudaq/convert_eudaq.dag
"""
import argparse
import json
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_REPO, "gaudi_jobs", "calibration"))
sys.path.insert(0, os.path.join(_REPO, "gaudi_jobs", "calibration", "condor"))
sys.path.insert(0, os.path.join(_REPO, "gaudi_jobs"))
sys.path.insert(0, _REPO)

from calib_run_utils import DEFAULT_RAW_BASE, raw_chunk_files
from condor_common import DEFAULT_CONVERTED_DIR, OPTIONS_DIR, chunks_dir, condor_log_dir, \
    env_wrapper_preamble, mkdirs_line, write_executable, write_text

_CLEANUP = os.path.join(_REPO, "gaudi_jobs", "calibration", "condor", "cleanup_job_logs.sh")

_VALID_TREE_CHECK = (
    'python3 -c "import ROOT,sys; ROOT.gErrorIgnoreLevel=ROOT.kFatal; '
    "f=ROOT.TFile.Open(sys.argv[1]); t=f.Get('siwecaldecoded') if f and not f.IsZombie() else None; "
    "sys.exit(0 if t and t.GetEntries() > 0 else 1)\" "
)
# Entries > 0, not just "the tree exists": production already holds eudaq
# chunks decoded before this was understood -- valid ROOT files, 0-entry
# trees. is_valid_chunk() elsewhere treats those as fine (structural check
# only); this pipeline exists specifically to redo them, so its own
# idempotency check must not wave them through as "already decoded".


def _convert_sh(out_dir):
    # Same grouped-per-chunk approach as generate_reco_dag.py's CONVERT stage
    # (a job decodes a GROUP of chunks, each still in its own k4run process --
    # see that script for why: ~14s of fixed per-job overhead vs ~1.5s of
    # actual decoding). EUDAQ_FLAG is per-RUN (all rows of one run's grouplist
    # carry the same value, from the route map's "eudaq_format"), not
    # hardcoded, since it differs across runs that are both routed "raw".
    content = f"""#!/bin/bash
set -eo pipefail
# convert.sh <grouplist> <run_settings> <eudaq_flag>
#   grouplist: text file, one "<raw_chunk> <out_decoded>" pair per line
#   eudaq_flag: 1 to set RAW2ROOT_EUDAQ (route-map "eudaq_format"), else 0
GROUPLIST="$1"; RUN_SETTINGS="$2"; EUDAQ_FLAG="$3"

""" + env_wrapper_preamble() + f"""
N=0; SKIPPED=0
while read -r CHUNK OUT_DECODED; do
  [ -z "$CHUNK" ] && continue
  # Idempotent: a chunk that already opens with its tree is left alone, so a
  # partial submit or a DAG rescue only re-does what is missing.
  if [ -s "$OUT_DECODED" ] && {_VALID_TREE_CHECK}"$OUT_DECODED" >/dev/null 2>&1; then
    SKIPPED=$((SKIPPED+1))
    continue
  fi
""" + "  " + mkdirs_line('$(dirname "$OUT_DECODED")') + f"""
  RAW_FILES="$CHUNK" RAW2ROOT_OUT="$OUT_DECODED" RAW2ROOT_RUN_SETTINGS_FILE="$RUN_SETTINGS" \\
  RAW2ROOT_EUDAQ="$EUDAQ_FLAG" \\
    k4run "{OPTIONS_DIR}/run_raw2root.py"
  if ! {_VALID_TREE_CHECK}"$OUT_DECODED" >/dev/null 2>&1; then
    echo "ERROR: $OUT_DECODED was written but has no siwecaldecoded tree" >&2
    exit 1
  fi
  N=$((N+1))
done < "$GROUPLIST"
echo "decoded $N chunk(s), skipped $SKIPPED already done"
"""
    write_executable(os.path.join(out_dir, "convert.sh"), content)


def _convert_sub(out_dir, log_dir, run, request_memory, job_flavour):
    content = f"""universe                = vanilla
executable              = {out_dir}/convert.sh
arguments               = "$(grouplist) $(run_settings) $(eudaq_flag)"
log                     = {log_dir}/convert_{run}_$(Process).log
output                  = {log_dir}/convert_{run}_$(Process).out
error                   = {log_dir}/convert_{run}_$(Process).err
request_cpus            = 1
request_memory          = ifthenelse(NumJobStarts <= 0, {request_memory}, {request_memory} * (NumJobStarts + 1))
request_disk            = 4000M
should_transfer_files   = NO
getenv                  = False
periodic_release        = (JobStatus == 5) && (HoldReasonCode == 34) && (NumJobStarts < 5)
+JobFlavour             = "{job_flavour}"
queue grouplist,run_settings,eudaq_flag from {out_dir}/grouplist_{run}.txt
"""
    write_text(os.path.join(out_dir, f"convert_{run}.sub"), content)


def _convert_lcio_sh(out_dir):
    # LCIO twin of _convert_sh: a run whose raw .bin data the raw decoder
    # cannot parse at all (see gaudi_jobs/decode_lcio_runs.py) decodes instead
    # from its pre-existing LCIO file via EcalLcioDecoder. One job per run,
    # one chunk_0000.root -- an LCIO run is a single .slcio file, not chunked.
    content = f"""#!/bin/bash
set -eo pipefail
# convert_lcio.sh <lcio_in> <out_decoded> <run_settings>
LCIO_IN="$1"; OUT_DECODED="$2"; RUN_SETTINGS="$3"

""" + env_wrapper_preamble() + f"""
if [ -s "$OUT_DECODED" ] && {_VALID_TREE_CHECK}"$OUT_DECODED" >/dev/null 2>&1; then
  echo "already decoded -> $OUT_DECODED"
  exit 0
fi
""" + mkdirs_line('$(dirname "$OUT_DECODED")') + f"""
LCIO_IN="$LCIO_IN" LCIO_OUT="$OUT_DECODED" LCIO_RUN_SETTINGS_FILE="$RUN_SETTINGS" \\
  k4run "{OPTIONS_DIR}/run_lcio_decode.py"
if ! {_VALID_TREE_CHECK}"$OUT_DECODED" >/dev/null 2>&1; then
  echo "ERROR: $OUT_DECODED was written but has no siwecaldecoded tree" >&2
  exit 1
fi
"""
    write_executable(os.path.join(out_dir, "convert_lcio.sh"), content)


def _convert_lcio_sub(out_dir, log_dir, run, request_memory, job_flavour):
    content = f"""universe                = vanilla
executable              = {out_dir}/convert_lcio.sh
arguments               = "$(lcio_in) $(out_decoded) $(run_settings)"
log                     = {log_dir}/convert_lcio_{run}.log
output                  = {log_dir}/convert_lcio_{run}.out
error                   = {log_dir}/convert_lcio_{run}.err
request_cpus            = 1
request_memory          = ifthenelse(NumJobStarts <= 0, {request_memory}, {request_memory} * (NumJobStarts + 1))
request_disk            = 8000M
should_transfer_files   = NO
getenv                  = False
periodic_release        = (JobStatus == 5) && (HoldReasonCode == 34) && (NumJobStarts < 5)
+JobFlavour             = "{job_flavour}"
queue
"""
    write_text(os.path.join(out_dir, f"convert_lcio_{run}.sub"), content)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Generate a Condor/DAGMan pipeline that CONVERTs every eudaq run "
                    "(raw -> siwecaldecoded chunks), per a pre-computed --route-map.")
    p.add_argument("--route-map", required=True,
                   help="JSON file: {run: {route: raw|lcio|failed, eudaq_format, lcio_in, ...}}. "
                        "Build with the classification helper (test-decode + entries check) -- "
                        "do NOT guess routing from the run name or raw filename pattern, neither "
                        "reliably predicts which decoder actually works (see module docstring).")
    p.add_argument("--runs", default=None,
                   help="Comma-separated run name(s) to restrict to. Default: every run in --route-map.")
    p.add_argument("--raw-base", default=DEFAULT_RAW_BASE,
                   help=f"Base dir the raw runs live in. Default: {DEFAULT_RAW_BASE}")
    p.add_argument("--converted-dir", default=DEFAULT_CONVERTED_DIR,
                   help=f"Where the decoded <RUN>/chunks/ go. Default: {DEFAULT_CONVERTED_DIR}")
    p.add_argument("--out-dir", required=True, help="Directory to write the DAG, subs, wrappers and logs into.")
    p.add_argument("--chunks-per-job", type=int, default=20,
                   help="Raw chunks decoded by ONE Condor job on the raw path, each in its own k4run "
                        "(default 20, same rationale as generate_reco_dag.py).")
    p.add_argument("--convert-request-memory", type=int, default=2000,
                   help="request_memory in MB for a raw-path CONVERT job (default 2000, same measured basis "
                        "as generate_reco_dag.py).")
    p.add_argument("--lcio-request-memory", type=int, default=3000,
                   help="request_memory in MB for an LCIO-path CONVERT job (default 3000 -- LCIO runs decode "
                        "tens of thousands of events into one multi-GB file in a single process).")
    p.add_argument("--convert-job-flavour", default="microcentury", help="+JobFlavour for CONVERT (default 1h).")
    p.add_argument("--lcio-job-flavour", default="workday",
                   help="+JobFlavour for the LCIO path (default 8h -- a single-process decode of tens of "
                        "thousands of events, not a small chunk).")
    p.add_argument("--keep-job-logs", action="store_true",
                   help="Keep every job's .out/.err. By default they are deleted when the node succeeds, to "
                        "protect the AFS logs/ directory's entry limit (see calibration/condor/README.md).")
    args = p.parse_args(argv)

    with open(args.route_map) as fh:
        route_map = json.load(fh)

    if args.runs:
        runs = [r.strip() for r in args.runs.split(",") if r.strip()]
        missing = [r for r in runs if r not in route_map]
        if missing:
            print(f"ERROR: not in --route-map: {', '.join(missing)}", file=sys.stderr)
            return 1
    else:
        runs = sorted(route_map)

    out_dir = os.path.abspath(args.out_dir)
    log_dir = condor_log_dir(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)
    groups_dir = os.path.join(out_dir, "groups")
    os.makedirs(groups_dir, exist_ok=True)

    dag = []
    n_raw = n_lcio = n_skipped = 0
    for run in runs:
        info = route_map[run]
        route = info.get("route")
        raw_dir = os.path.join(args.raw_base, run)
        run_settings = os.path.join(raw_dir, "Run_Settings.txt")
        cdir = chunks_dir(args.converted_dir, run)

        if route == "raw":
            raw_chunks = raw_chunk_files(raw_dir, run)
            eudaq_flag = "1" if info["eudaq_format"] else "0"
            os.makedirs(cdir, exist_ok=True)
            pairs = [(c, os.path.join(cdir, f"chunk_{i:04d}.root")) for i, c in enumerate(raw_chunks)]
            rows = []
            for g, start in enumerate(range(0, len(pairs), args.chunks_per_job)):
                group = pairs[start:start + args.chunks_per_job]
                gpath = os.path.join(groups_dir, f"{run}_g{g:04d}.txt")
                write_text(gpath, "\n".join(f"{raw} {out}" for raw, out in group) + "\n")
                rows.append(f"{gpath},{run_settings},{eudaq_flag}")
            write_text(os.path.join(out_dir, f"grouplist_{run}.txt"), "\n".join(rows) + "\n")
            _convert_sub(out_dir, log_dir, run, args.convert_request_memory, args.convert_job_flavour)

            dag.append(f"JOB convert_{run} {out_dir}/convert_{run}.sub")
            dag.append(f"RETRY convert_{run} 3")
            if not args.keep_job_logs:
                dag.append(f"SCRIPT POST convert_{run} {_CLEANUP} $RETURN {os.path.join(log_dir, f'convert_{run}')}")
            dag.append("")
            n_raw += 1
            print(f"[{run}] {len(raw_chunks)} raw chunk(s) in {len(rows)} job(s), "
                  f"eudaq_format={info['eudaq_format']} -> {cdir}/  (raw path)")
            continue

        if route == "lcio":
            lcio_in = info["lcio_in"]
            os.makedirs(cdir, exist_ok=True)
            out_decoded = os.path.join(cdir, "chunk_0000.root")
            _convert_lcio_sub(out_dir, log_dir, run, args.lcio_request_memory, args.lcio_job_flavour)
            dag.append(f"JOB convert_lcio_{run} {out_dir}/convert_lcio_{run}.sub")
            dag.append(f'VARS convert_lcio_{run} lcio_in="{lcio_in}" out_decoded="{out_decoded}" '
                       f'run_settings="{run_settings}"')
            dag.append(f"RETRY convert_lcio_{run} 2")
            if not args.keep_job_logs:
                dag.append(f"SCRIPT POST convert_lcio_{run} {_CLEANUP} $RETURN "
                           f"{os.path.join(log_dir, f'convert_lcio_{run}')}")
            dag.append("")
            n_lcio += 1
            print(f"[{run}] LCIO path: {lcio_in} -> {out_decoded}  (LCIO path)")
            continue

        n_skipped += 1
        print(f"[{run}] SKIPPED -- route-map says {route!r} ({info.get('reason', 'no reason given')})",
              file=sys.stderr)

    _convert_sh(out_dir)
    if n_lcio:
        _convert_lcio_sh(out_dir)

    dag_path = os.path.join(out_dir, "convert_eudaq.dag")
    write_text(dag_path, "\n".join(dag) + "\n")
    print(f"\n{n_raw} run(s) raw path, {n_lcio} run(s) LCIO path, {n_skipped} run(s) skipped (unresolved).")
    print(f"Submit with:\n  condor_submit_dag {dag_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
