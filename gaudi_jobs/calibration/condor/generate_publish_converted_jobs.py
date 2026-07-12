#!/usr/bin/env python
"""
Generate the PUBLISH stage: one HTCondor job per run that hadds that run's
per-chunk decoded files into the single converted file the rest of the
project expects, ``<converted_dir>/<RUN>/<RUN>.root``.

Why this stage exists. CONVERT (generate_convert_jobs.py) decodes each raw
chunk in its own job, so a calibration run ends up as N per-chunk ROOT files
under the Condor scratch area's ``decoded/``. That is a parallelisation
artefact, not a converted run: every other converted run in the experiment --
the ones in rundata_converted_old/, the ones written by gaudi_jobs/run000012
and run000013, and the one run_calibration_batch.py writes -- is ONE directory
per run holding ONE <RUN>.root. Without this stage the calibration runs are
the odd ones out, invisible to anything that looks in the converted area.

Nothing downstream in the calibration pipeline needs the merged file (Fill
reads the chunks directly, which is the whole point of chunking them), so this
stage is about making the converted runs findable and uniform, not about
feeding the fits.

No dependencies between runs, so no DAGMan: a single multi-queue .sub file
(``queue run,out_file,decoded_run_dir from runlist.txt``) plus its wrapper.
Only WRITES files -- never calls condor_submit itself.

Usage::

    python gaudi_jobs/calibration/condor/generate_publish_converted_jobs.py \\
        --runs TB2026CERN_run_000060,TB2026CERN_run_000061 \\
        --out-dir gaudi_jobs/calibration/condor/generated/publish_th220

    condor_submit gaudi_jobs/calibration/condor/generated/publish_th220/publish.sub
"""
import argparse
import glob
import os
import sys

from condor_common import DEFAULT_CONVERTED_DIR, DEFAULT_FILL_SCRATCH_DIR, condor_log_dir, converted_run_file, \
    decoded_dir, env_wrapper_preamble, mkdirs_line, write_executable, write_text


def _parse_runs(spec):
    """Bare run names from a comma-separated --runs value.

    Deliberately NOT calib_run_utils.parse_run_folder_list: that resolves each
    name against the raw data area, and this stage never touches raw -- it only
    reads what CONVERT already decoded.
    """
    runs = [r.strip().rstrip(os.sep).split(os.sep)[-1] for r in spec.split(",") if r.strip()]
    if not runs:
        raise SystemExit("ERROR: --runs is empty")
    return runs


def _publish_sh(out_dir):
    # hadd streams TTrees sequentially, so it runs fine straight off EOS here
    # (~30 MB/s measured) -- unlike the histogram merges, whose pathology was
    # hundreds of thousands of tiny keys, not volume. No local staging needed.
    #
    # Idempotent: an existing, valid output is left alone, so a partially
    # completed submit can simply be resubmitted.
    # The skip-check must come AFTER the preamble: it needs key4hep's ROOT, and
    # a bare system python3 would just fail the import and re-hadd every time.
    content = f"""#!/bin/bash
set -eo pipefail
# publish.sh <run> <out_file> <decoded_run_dir>
RUN="$1"; OUT_FILE="$2"; DECODED_DIR="$3"

""" + env_wrapper_preamble() + f"""
if [ -f "$OUT_FILE" ] && python3 -c "
import ROOT, sys
f = ROOT.TFile.Open(sys.argv[1])
sys.exit(0 if f and not f.IsZombie() and f.Get('siwecaldecoded') else 1)" "$OUT_FILE" 2>/dev/null; then
  echo "$RUN: already published, skipping"
  exit 0
fi

""" + mkdirs_line('$(dirname "$OUT_FILE")') + f"""
hadd -f -v 0 "$OUT_FILE" "$DECODED_DIR"/*.root
echo "$RUN: published -> $OUT_FILE"
"""
    write_executable(os.path.join(out_dir, "publish.sh"), content)


def _publish_sub(out_dir, log_dir, request_memory, job_flavour):
    content = f"""universe                = vanilla
executable              = {out_dir}/publish.sh
arguments               = "$(run) $(out_file) $(decoded_run_dir)"
log                     = {log_dir}/publish_$(run).log
output                  = {log_dir}/publish_$(run).out
error                   = {log_dir}/publish_$(run).err
request_cpus            = 1
request_memory          = {request_memory}M
request_disk            = 4000M
should_transfer_files   = NO
getenv                  = False
+JobFlavour             = "{job_flavour}"
queue run,out_file,decoded_run_dir from {out_dir}/runlist.txt
"""
    write_text(os.path.join(out_dir, "publish.sub"), content)


def main(argv=None):
    p = argparse.ArgumentParser(description="Generate the PUBLISH (decoded chunks -> <RUN>.root) Condor job set.")
    p.add_argument("--runs", required=True,
                   help="Comma-separated run folder names (e.g. TB2026CERN_run_000060,TB2026CERN_eudaq_run_000146).")
    p.add_argument("--fill-scratch-dir", default=DEFAULT_FILL_SCRATCH_DIR,
                   help=f"EOS scratch area holding decoded/. Default: {DEFAULT_FILL_SCRATCH_DIR}")
    p.add_argument("--converted-dir", default=DEFAULT_CONVERTED_DIR,
                   help=f"Where to publish <RUN>/<RUN>.root. Default: {DEFAULT_CONVERTED_DIR}")
    p.add_argument("--out-dir", required=True,
                   help="Directory to write publish.sub/publish.sh/runlist.txt/logs/ into.")
    p.add_argument("--request-memory", type=int, default=4000,
                   help="request_memory in MB (default 4000). hadd streams, so memory is flat in the number of chunks.")
    p.add_argument("--job-flavour", default="workday",
                   help="HTCondor +JobFlavour (default workday, 8h -- the largest run here is ~57 GB at ~30 MB/s).")
    args = p.parse_args(argv)

    runs = _parse_runs(args.runs)
    decoded_base = decoded_dir(args.fill_scratch_dir)
    out_dir = os.path.abspath(args.out_dir)
    log_dir = condor_log_dir(out_dir)
    os.makedirs(log_dir, exist_ok=True)

    rows = []
    for run in runs:
        run_decoded = os.path.join(decoded_base, run)
        chunks = sorted(glob.glob(os.path.join(run_decoded, "*.root")))
        if not chunks:
            print(f"ERROR: no decoded chunks for {run} in {run_decoded} -- run CONVERT first.", file=sys.stderr)
            return 1
        rows.append(f"{run},{converted_run_file(args.converted_dir, run)},{run_decoded}")
        print(f"  {run}: {len(chunks)} chunk(s) -> {converted_run_file(args.converted_dir, run)}")

    write_text(os.path.join(out_dir, "runlist.txt"), "\n".join(rows) + "\n")
    _publish_sh(out_dir)
    _publish_sub(out_dir, log_dir, args.request_memory, args.job_flavour)

    print(f"\n{len(rows)} job(s) written. Submit with:\n  condor_submit {os.path.join(out_dir, 'publish.sub')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
