#!/usr/bin/env python
"""
Generate the full raw -> EDM4hep Condor pipeline for any number of runs, as one
DAG: CONVERT (one job per raw chunk, all in parallel) then RECO (one job per
run: event building + PID).

Why a DAG and not the serial drivers (gaudi_jobs/run_full_pipeline_*.py). Those
decode every raw chunk of a run inside ONE k4run process and write a single
multi-GB siwecaldecoded file, which the event builder then reads. That file is
pure overhead:

  * it doubles the decoded data on EOS (~10 GB for one physics run, ~390 GB
    across the calibration runs) for something no analysis ever opens;
  * writing it is serial, so a run takes ~20 min of wall clock that 300 parallel
    chunk jobs do in a couple of minutes;
  * closing a 10 GB file over XRootD times out often enough to be a real
    operational problem (it stalled run_000013 three times in one night), and a
    failed Close leaves a 0-byte file that fails the next stage;
  * past ~100 GB it is not even possible: ROOT's TTree::fgMaxTreeSize rolls the
    tree over into a second file, which no merger can produce as one file.

So: decode chunks in parallel, keep them, and let the event builder chain them
(EcalEventBuilder's InputFiles). The per-run OUTPUT is still exactly one
ecal_<run>.root and one ecal_<run>.edm4hep.root -- only the intermediate merge
disappears.

Layout, per run, under --converted-dir:

    <RUN>/chunks/chunk_0000.root ...   decoded (CONVERT)
    <RUN>/ecal_<RUN>.root              event-built (RECO)
    <RUN>/ecal_<RUN>.edm4hep.root      PID/EDM4hep (RECO)

Calibration tables are resolved per run from its own ThresholdDAC (read from the
run's Run_Settings.txt), exactly as the serial drivers do, so runs taken at
different thresholds can be submitted together.

Only WRITES files -- never calls condor_submit itself.

Usage::

    python gaudi_jobs/condor/generate_reco_dag.py \\
        --runs TB2026CERN_run_000012,TB2026CERN_run_000013 \\
        --out-dir gaudi_jobs/condor/generated/physics

    condor_submit_dag gaudi_jobs/condor/generated/physics/reco.dag
"""
import argparse
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
sys.path.insert(0, os.path.join(_REPO, "gaudi_jobs", "calibration"))
sys.path.insert(0, os.path.join(_REPO, "gaudi_jobs", "calibration", "condor"))
sys.path.insert(0, _REPO)

from calib_run_utils import parse_run_folder_list, raw_chunk_files
from condor_common import DEFAULT_CONVERTED_DIR, OPTIONS_DIR, chunks_dir, condor_log_dir, env_wrapper_preamble, \
    mkdirs_line, write_executable, write_text
from siwecal_common import paths
from siwecal_eventbuilder.cli import resolve_gaudi_calib_files
from siwecal_eventbuilder.run_settings import read_threshold_dac

_CLEANUP = os.path.join(_REPO, "gaudi_jobs", "calibration", "condor", "cleanup_job_logs.sh")


def _convert_sh(out_dir):
    content = f"""#!/bin/bash
set -eo pipefail
# convert.sh <chunk> <out_decoded> <run_settings>
CHUNK="$1"; OUT_DECODED="$2"; RUN_SETTINGS="$3"

# Idempotent: a chunk already decoded is left alone, so a partially completed
# submit can simply be resubmitted (and a DAG rescue re-runs only what is missing).
if [ -s "$OUT_DECODED" ]; then
  echo "already decoded, skipping: $OUT_DECODED"
  exit 0
fi

""" + env_wrapper_preamble() + f"""
""" + mkdirs_line('$(dirname "$OUT_DECODED")') + f"""
RAW_FILES="$CHUNK" RAW2ROOT_OUT="$OUT_DECODED" RAW2ROOT_RUN_SETTINGS_FILE="$RUN_SETTINGS" \\
  k4run "{OPTIONS_DIR}/run_raw2root.py"
"""
    write_executable(os.path.join(out_dir, "convert.sh"), content)


def _reco_sh(out_dir):
    # One job per run: event building (chaining that run's chunks) then PID.
    # Two k4run processes, not one: EcalToEDM4hep/EcalPidTransformer are
    # k4FWCore components whose EvtMax must be fixed before the process starts,
    # and that count is only known once the event builder has run. run_pid.py
    # reads it from the ecal file itself, so no driver is needed in between.
    content = f"""#!/bin/bash
set -eo pipefail
# reco.sh <run> <chunks_glob> <ecal_out> <pid_out> <ped> <mip> <ped_lg> <mip_lg> <padmap> <padmap_ovr> <slab_z> <run_id>
RUN="$1"; CHUNKS="$2"; ECAL_OUT="$3"; PID_OUT="$4"
PED="$5"; MIP="$6"; PED_LG="$7"; MIP_LG="$8"
PADMAP="$9"; PADMAP_OVR="${{10}}"; SLAB_Z="${{11}}"; RUN_ID="${{12}}"

""" + env_wrapper_preamble() + f"""
""" + mkdirs_line('$(dirname "$ECAL_OUT")') + f"""
echo "[reco] $RUN: event building from $CHUNKS"
EVBLD_INPUT_FILES="$CHUNKS" \\
EVBLD_OUTPUT="$ECAL_OUT" \\
EVBLD_RUN_ID="$RUN_ID" \\
EVBLD_PEDESTAL_FILE="$PED" \\
EVBLD_MIP_FILE="$MIP" \\
EVBLD_PEDESTAL_FILE_LOWGAIN="$PED_LG" \\
EVBLD_MIP_FILE_LOWGAIN="$MIP_LG" \\
EVBLD_PADMAP_DEFAULT="$PADMAP" \\
EVBLD_PADMAP_SLAB_OVERRIDES="$PADMAP_OVR" \\
EVBLD_SLAB_Z_FILE="$SLAB_Z" \\
  k4run "{OPTIONS_DIR}/run_event_builder.py"

echo "[reco] $RUN: PID/EDM4hep -> $PID_OUT"
ECAL_FILE="$ECAL_OUT" ECAL_PID_OUT="$PID_OUT" \\
  k4run "{OPTIONS_DIR}/run_pid.py"

echo "[reco] $RUN: done"
"""
    write_executable(os.path.join(out_dir, "reco.sh"), content)


def _convert_sub(out_dir, log_dir, run, request_memory, job_flavour):
    content = f"""universe                = vanilla
executable              = {out_dir}/convert.sh
arguments               = "$(chunk) $(out_decoded) $(run_settings)"
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
queue chunk,out_decoded,run_settings from {out_dir}/chunklist_{run}.txt
"""
    write_text(os.path.join(out_dir, f"convert_{run}.sub"), content)


def _reco_sub(out_dir, log_dir, request_memory, job_flavour):
    content = f"""universe                = vanilla
executable              = {out_dir}/reco.sh
arguments               = "$(run) '$(chunks)' $(ecal_out) $(pid_out) $(ped) $(mip) $(ped_lg) $(mip_lg) $(padmap) $(padmap_ovr) $(slab_z) $(run_id)"
log                     = {log_dir}/reco_$(run).log
output                  = {log_dir}/reco_$(run).out
error                   = {log_dir}/reco_$(run).err
request_cpus            = 1
request_memory          = ifthenelse(NumJobStarts <= 0, {request_memory}, {request_memory} * (NumJobStarts + 1))
request_disk            = 8000M
should_transfer_files   = NO
getenv                  = False
periodic_release        = (JobStatus == 5) && (HoldReasonCode == 34) && (NumJobStarts < 5)
+JobFlavour             = "{job_flavour}"
queue
"""
    write_text(os.path.join(out_dir, "reco.sub"), content)


def main(argv=None):
    p = argparse.ArgumentParser(description="Generate the CONVERT+RECO Condor DAG (raw -> ecal -> EDM4hep).")
    p.add_argument("--runs", required=True,
                   help="Comma-separated run folder path(s) or bare run name(s).")
    p.add_argument("--raw-base", default=None, help="Base dir to resolve bare run names against.")
    p.add_argument("--converted-dir", default=DEFAULT_CONVERTED_DIR,
                   help=f"Where <RUN>/chunks/ and <RUN>/ecal_*.root go. Default: {DEFAULT_CONVERTED_DIR}")
    p.add_argument("--out-dir", required=True, help="Directory to write the DAG, subs, wrappers and logs into.")
    p.add_argument("--convert-request-memory", type=int, default=7000,
                   help="request_memory in MB for a CONVERT (one-chunk decode) job. Default 7000; see "
                        "generate_convert_jobs.py for why it is this high.")
    p.add_argument("--reco-request-memory", type=int, default=8000,
                   help="request_memory in MB for a RECO job (event building holds a run's events in memory).")
    p.add_argument("--convert-job-flavour", default="microcentury", help="+JobFlavour for CONVERT (default 1h).")
    p.add_argument("--reco-job-flavour", default="workday", help="+JobFlavour for RECO (default 8h).")
    p.add_argument("--keep-job-logs", action="store_true",
                   help="Keep every job's .out/.err. By default they are deleted when the node succeeds, to "
                        "protect the AFS logs/ directory's entry limit (see calibration/condor/README.md).")
    args = p.parse_args(argv)

    kwargs = {"default_base": args.raw_base} if args.raw_base else {}
    run_folders = parse_run_folder_list(args.runs, **kwargs)

    out_dir = os.path.abspath(args.out_dir)
    log_dir = condor_log_dir(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(log_dir, exist_ok=True)

    mappings = paths.geometry_dir()
    padmap = os.path.join(mappings, "fev10_rotate_chip_channel_x_y_mapping.txt")
    padmap_slab12 = os.path.join(mappings, "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
    slab_z = os.path.join(_REPO, "event_display", "conversion", "slab_z_positions.yml")

    dag = []
    for run, raw_dir in run_folders:
        raw_chunks = raw_chunk_files(raw_dir, run)
        if not raw_chunks:
            print(f"ERROR: no raw chunks for {run} in {raw_dir}", file=sys.stderr)
            return 1
        run_settings = os.path.join(raw_dir, "Run_Settings.txt")

        th = read_threshold_dac(run_settings)
        if th < 0:
            print(f"ERROR: cannot read ThresholdDAC from {run_settings}", file=sys.stderr)
            return 1
        ped, mip, ped_lg, mip_lg = resolve_gaudi_calib_files(th)

        cdir = chunks_dir(args.converted_dir, run)
        os.makedirs(cdir, exist_ok=True)
        rows = [f"{c},{os.path.join(cdir, f'chunk_{i:04d}.root')},{run_settings}"
                for i, c in enumerate(raw_chunks)]
        write_text(os.path.join(out_dir, f"chunklist_{run}.txt"), "\n".join(rows) + "\n")
        _convert_sub(out_dir, log_dir, run, args.convert_request_memory, args.convert_job_flavour)

        digits = "".join(ch for ch in run.split("_")[-1] if ch.isdigit())
        run_id = int(digits) if digits else -1
        ecal_out = os.path.join(args.converted_dir, run, f"ecal_{run}.root")
        pid_out = os.path.join(args.converted_dir, run, f"ecal_{run}.edm4hep.root")

        dag.append(f"JOB convert_{run} {out_dir}/convert_{run}.sub")
        dag.append(f"RETRY convert_{run} 3")
        dag.append(f"JOB reco_{run} {out_dir}/reco.sub")
        dag.append(
            f'VARS reco_{run} run="{run}" chunks="{os.path.join(cdir, "chunk_*.root")}" '
            f'ecal_out="{ecal_out}" pid_out="{pid_out}" ped="{ped}" mip="{mip}" '
            f'ped_lg="{ped_lg}" mip_lg="{mip_lg}" padmap="{padmap}" '
            f'padmap_ovr="12:{padmap_slab12}" slab_z="{slab_z}" run_id="{run_id}"')
        dag.append(f"RETRY reco_{run} 2")
        dag.append(f"PARENT convert_{run} CHILD reco_{run}")
        if not args.keep_job_logs:
            dag.append(f"SCRIPT POST reco_{run} {_CLEANUP} $RETURN {os.path.join(log_dir, f'reco_{run}')}")
        dag.append("")

        print(f"[{run}] th={th}  {len(raw_chunks)} chunk(s) -> {cdir}/")
        print(f"          calib: {os.path.basename(ped)} / {os.path.basename(mip)}"
              f"{'' if ped_lg else '   (no low-gain tables: no saturation recovery)'}")

    _convert_sh(out_dir)
    _reco_sh(out_dir)
    _reco_sub(out_dir, log_dir, args.reco_request_memory, args.reco_job_flavour)
    dag_path = os.path.join(out_dir, "reco.dag")
    write_text(dag_path, "\n".join(dag) + "\n")
    print(f"\n{len(run_folders)} run(s). Submit with:\n  condor_submit_dag {dag_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
