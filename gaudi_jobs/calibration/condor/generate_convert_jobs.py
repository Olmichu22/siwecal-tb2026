#!/usr/bin/env python
"""
Generate the CONVERT stage of the Condor calibration pipeline: one
independent HTCondor job per raw DAQ chunk file, each running EcalRawDecoder
on just that chunk. This is the ONLY stage of the whole pipeline that opens
a file under the read-only raw data area -- see gaudi_jobs/calibration/
README.md for the full two-stage design (CONVERT here, CALIBRATION in
generate_dag.py).

No dependencies between chunks, so no DAGMan is needed: this writes a single
multi-queue .sub file (``queue chunk,out_decoded,run_settings from
chunklist.txt``) plus its executable wrapper. Only WRITES files -- never
calls condor_submit itself; submitting potentially thousands of jobs to the
shared CERN batch farm is left as a deliberate, explicit step for the user.

Usage::

    python gaudi_jobs/calibration/condor/generate_convert_jobs.py \\
        --runs TB2026CERN_run_000060,TB2026CERN_run_000061 \\
        --out-dir gaudi_jobs/calibration/condor/generated/th220

    condor_submit gaudi_jobs/calibration/condor/generated/th220/convert.sub
"""
import argparse
import os
import sys

# calib_run_utils lives one directory up (gaudi_jobs/calibration/); must be
# on sys.path BEFORE importing it -- condor_common.py is a sibling of this
# script so it doesn't need this.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from calib_run_utils import parse_run_folder_list, raw_chunk_files
from condor_common import DEFAULT_FILL_SCRATCH_DIR, OPTIONS_DIR, condor_log_dir, decoded_dir, env_wrapper_preamble, \
    mkdirs_line, write_executable, write_text


def _convert_sh(out_dir):
    content = f"""#!/bin/bash
set -eo pipefail
# convert.sh <chunk> <out_decoded> <run_settings>
CHUNK="$1"; OUT_DECODED="$2"; RUN_SETTINGS="$3"
""" + env_wrapper_preamble() + f"""
""" + mkdirs_line('$(dirname "$OUT_DECODED")') + f"""RAW_FILES="$CHUNK" RAW2ROOT_OUT="$OUT_DECODED" RAW2ROOT_RUN_SETTINGS_FILE="$RUN_SETTINGS" \\
  k4run "{OPTIONS_DIR}/run_raw2root.py"
"""
    write_executable(os.path.join(out_dir, "convert.sh"), content)


def _convert_sub(out_dir, log_dir, request_memory, job_flavour):
    content = f"""universe                = vanilla
executable              = {out_dir}/convert.sh
arguments               = "$(chunk) $(out_decoded) $(run_settings)"
log                     = {log_dir}/convert_$(Process).log
output                  = {log_dir}/convert_$(Process).out
error                   = {log_dir}/convert_$(Process).err
request_cpus            = 1
request_memory          = {request_memory}M
request_disk            = 4000M
should_transfer_files   = NO
getenv                  = False
+JobFlavour             = "{job_flavour}"
queue chunk,out_decoded,run_settings from {out_dir}/chunklist.txt
"""
    write_text(os.path.join(out_dir, "convert.sub"), content)


def main(argv=None):
    p = argparse.ArgumentParser(description="Generate the CONVERT (raw2root per chunk) Condor job set.")
    p.add_argument("--runs", required=True,
                   help="Comma-separated run folder path(s) or bare run name(s), same format as "
                        "run_calibration_batch.py --runs. No auto-discovery by threshold.")
    p.add_argument("--raw-base", default=None,
                   help="Base directory to resolve bare run names in --runs against. "
                        "Default: calib_run_utils.DEFAULT_RAW_BASE")
    p.add_argument("--fill-scratch-dir", default=DEFAULT_FILL_SCRATCH_DIR,
                   help=f"EOS scratch area (decoded/ and hist/ live under it). Default: {DEFAULT_FILL_SCRATCH_DIR}")
    p.add_argument("--out-dir", required=True,
                   help="Directory to write convert.sub/convert.sh/chunklist.txt/logs/ into.")
    p.add_argument("--request-memory", type=int, default=7000,
                   help="request_memory in MB (default 7000). EcalRawDecoder buffers a whole chunk's decoded "
                        "acquisitions in memory before flushing to disk (~5.5MB per acquisition cycle), and chunk "
                        "density varies a lot by run -- observed real Condor OOM holds at 3.5-5.3GB actual usage "
                        "for a single chunk of a dense run (~1000 acquisitions/chunk), well above the low-density "
                        "~500MB seen in small test runs. Condor's cgroup can also round a low request up to a "
                        "pool minimum (observed: a 2000M request was enforced as a 3000M limit), so don't rely on "
                        "a low request being honored exactly either.")
    p.add_argument("--job-flavour", default="microcentury",
                   help="HTCondor +JobFlavour (default microcentury, 1h -- decoding one ~11-18MB chunk should take "
                        "minutes).")
    args = p.parse_args(argv)

    kwargs = {"default_base": args.raw_base} if args.raw_base else {}
    run_folders = parse_run_folder_list(args.runs, **kwargs)

    decoded_base = decoded_dir(args.fill_scratch_dir)
    # Condor logs stay on AFS (CERN standard schedds reject /eos log paths);
    # each convert job's .out is only a handful of lines, so they don't pile up.
    log_dir = condor_log_dir(args.out_dir)
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    lines = []
    n_chunks = 0
    for run_name, raw_dir in run_folders:
        chunks = raw_chunk_files(raw_dir, run_name)
        run_settings = os.path.join(raw_dir, "Run_Settings.txt")
        # Pre-create the per-run output directory now (os.makedirs works
        # fine on the EOS FUSE mount; shell `mkdir -p` on the worker node
        # does not -- see condor_common.py) so convert.sh never needs to
        # create more than the leaf itself.
        os.makedirs(os.path.join(decoded_base, run_name), exist_ok=True)
        for idx, chunk in enumerate(chunks):
            out_decoded = os.path.join(decoded_base, run_name, f"chunk_{idx:04d}.root")
            lines.append(f"{chunk},{out_decoded},{run_settings}")
            n_chunks += 1
        print(f"[convert] {run_name}: {len(chunks)} chunk(s) -> {os.path.join(decoded_base, run_name)}/")

    write_text(os.path.join(args.out_dir, "chunklist.txt"), "\n".join(lines) + "\n")
    _convert_sh(args.out_dir)
    _convert_sub(args.out_dir, log_dir, args.request_memory, args.job_flavour)

    print(f"\n[convert] {len(run_folders)} run(s), {n_chunks} chunk(s) total -> {n_chunks} Condor job(s)")
    print(f"[convert] Wrote {args.out_dir}/convert.sub (+ convert.sh, chunklist.txt, logs/)")
    print(f"\nTo submit:\n  condor_submit {args.out_dir}/convert.sub")
    return 0


if __name__ == "__main__":
    sys.exit(main())
