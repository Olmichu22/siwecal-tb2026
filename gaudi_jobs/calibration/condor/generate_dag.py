#!/usr/bin/env python
"""
Generate the CALIBRATION stage of the Condor pipeline: an HTCondor DAGMan
pipeline (fill per already-decoded chunk -> merge per run -> merge per
threshold -> fit) producing the same Pedestal_*.txt/MIP_*.txt tables
run_calibration_batch.py would, but scaled out across the batch farm.

This stage NEVER opens a file under the read-only raw data area -- it only
reads chunks already decoded by generate_convert_jobs.py's CONVERT stage
(validated to exist before anything is written; see _check_decoded_chunks
below). See gaudi_jobs/calibration/README.md for the full two-stage design.

Only WRITES files (.dag/.sub/.sh, hadd file lists, logs/) -- never calls
condor_submit_dag itself; submitting is left as a deliberate, explicit step.

Usage::

    python gaudi_jobs/calibration/condor/generate_dag.py \\
        --runs TB2026CERN_run_000060,TB2026CERN_run_000061 \\
        --fill-scratch-dir /eos/.../calib_fill_scratch \\
        --dag-dir gaudi_jobs/calibration/condor/generated/th220

    condor_submit_dag gaudi_jobs/calibration/condor/generated/th220/calibration_th220.dag
"""
import argparse
import os
import sys

# calib_run_utils/run_calibration_batch live one directory up
# (gaudi_jobs/calibration/); must be on sys.path BEFORE importing them --
# condor_common.py is a sibling of this script so it doesn't need this.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from siwecal_common import paths

from calib_run_utils import parse_run_folder_list, raw_chunk_files
from condor_common import DEFAULT_FILL_SCRATCH_DIR, OPTIONS_DIR, decoded_dir, env_wrapper_preamble, \
    hist_dir, mkdirs_line, write_executable, write_text
from run_calibration_batch import _check_threshold_consistency, _combined_label  # noqa: E402

_CHECK_SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "check_merged_histograms.sh")

# A single merge_run job hadd-ing thousands of per-chunk histogram files
# serially is what made the th220/th230 merges take the better part of a
# week: `hadd -n 100` re-merges the already-growing target against each new
# batch of 100, so cost grows with the number of batches processed so far,
# and it all runs on one CPU. Above this size, split into groups of
# _MERGE_GROUP_SIZE first -- each group is a single hadd batch (no
# regrowth), and the groups run as separate Condor jobs in parallel. The
# per-run merge then only has to combine a handful of already-merged group
# files instead of thousands of raw ones.
_MERGE_GROUP_SIZE = 100


def _expected_chunks(run_folders):
    """{run_name: [decoded_chunk_path, ...]} using the SAME enumeration order
    generate_convert_jobs.py used (sorted raw_chunk_files, 0-indexed) -- this
    is the correspondence between a raw chunk and its decoded output path."""
    expected = {}
    for run_name, raw_dir in run_folders:
        chunks = raw_chunk_files(raw_dir, run_name)
        expected[run_name] = [f"chunk_{idx:04d}.root" for idx in range(len(chunks))]
    return expected


def _check_decoded_chunks(run_folders, decoded_base):
    """Aborts with a clear message if generate_convert_jobs.py's CONVERT
    stage hasn't (fully) run yet for the requested runs."""
    expected = _expected_chunks(run_folders)
    missing = []
    for run_name, chunk_names in expected.items():
        for chunk_name in chunk_names:
            path = os.path.join(decoded_base, run_name, chunk_name)
            if not os.path.exists(path):
                missing.append(path)
    if missing:
        preview = "\n".join(f"  {m}" for m in missing[:20])
        more = f"\n  ... and {len(missing) - 20} more" if len(missing) > 20 else ""
        raise SystemExit(
            f"ERROR: {len(missing)} decoded chunk(s) missing under {decoded_base} -- run CONVERT first:\n"
            f"  python gaudi_jobs/calibration/condor/generate_convert_jobs.py --runs ... --out-dir ...\n"
            f"  condor_submit <out-dir>/convert.sub\n"
            f"Missing:\n{preview}{more}"
        )
    return expected


def _fill_sh(out_dir):
    content = f"""#!/bin/bash
set -eo pipefail
# fill.sh <decoded_chunk> <out_hist>
DECODED_CHUNK="$1"; OUT_HIST="$2"
""" + env_wrapper_preamble() + f"""
""" + mkdirs_line('$(dirname "$OUT_HIST")') + f"""CALIB_INPUT_FILES="$DECODED_CHUNK" CALIB_MODE=Fill CALIB_OUTPUT_HISTOGRAM_FILE="$OUT_HIST" \\
  k4run "{OPTIONS_DIR}/run_pedestal_mip.py"
"""
    write_executable(os.path.join(out_dir, "fill.sh"), content)


def _merge_sh(out_dir):
    content = """#!/bin/bash
set -eo pipefail
# merge.sh <filelist.txt> <out.root>
FILELIST="$1"; OUT="$2"
""" + env_wrapper_preamble() + """
""" + mkdirs_line('$(dirname "$OUT")') + """hadd -f -n 100 "$OUT" "@$FILELIST"
"""
    write_executable(os.path.join(out_dir, "merge.sh"), content)


def _fit_sh(out_dir):
    content = f"""#!/bin/bash
set -eo pipefail
# fit.sh <mode: FitPedestal|FitMip> <gain> <in_hist> <out_path> <diag_path (may be empty)>
MODE="$1"; GAIN="$2"; IN_HIST="$3"; OUT_PATH="$4"; DIAG_PATH="${{5:-}}"
""" + env_wrapper_preamble() + f"""
""" + mkdirs_line('$(dirname "$OUT_PATH")') + f"""if [ "$MODE" = "FitMip" ]; then export CALIB_OUTPUT_MIP_FILE="$OUT_PATH"; else export CALIB_OUTPUT_PEDESTAL_FILE="$OUT_PATH"; fi
CALIB_INPUT_HISTOGRAM_FILE="$IN_HIST" CALIB_MODE="$MODE" CALIB_GAIN="$GAIN" CALIB_DIAGNOSTICS_FILE="$DIAG_PATH" \\
  k4run "{OPTIONS_DIR}/run_pedestal_mip.py"
"""
    write_executable(os.path.join(out_dir, "fit.sh"), content)


def _write_subs(out_dir, fill_mem, fill_flavour, merge_mem, merge_flavour, fit_mem, fit_flavour):
    common = """universe                = vanilla
should_transfer_files   = NO
getenv                  = False
log                     = {out_dir}/logs/$(JOB).log
output                  = {out_dir}/logs/$(JOB).out
error                   = {out_dir}/logs/$(JOB).err
"""
    fill_sub = common.format(out_dir=out_dir) + f"""executable              = {out_dir}/fill.sh
arguments               = "$(in_decoded) $(out_hist)"
request_cpus            = 1
request_memory          = {fill_mem}M
request_disk            = 2000M
+JobFlavour             = "{fill_flavour}"
queue
"""
    merge_sub = common.format(out_dir=out_dir) + f"""executable              = {out_dir}/merge.sh
arguments               = "$(filelist) $(out)"
request_cpus            = 1
request_memory          = {merge_mem}M
request_disk            = 8000M
+JobFlavour             = "{merge_flavour}"
queue
"""
    fit_sub = common.format(out_dir=out_dir) + f"""executable              = {out_dir}/fit.sh
arguments               = "$(mode) $(gain) $(in_hist) $(out_path) $(diag_path)"
request_cpus            = 1
request_memory          = {fit_mem}M
request_disk            = 1000M
+JobFlavour             = "{fit_flavour}"
queue
"""
    write_text(os.path.join(out_dir, "fill.sub"), fill_sub)
    write_text(os.path.join(out_dir, "merge.sub"), merge_sub)
    write_text(os.path.join(out_dir, "fit.sub"), fit_sub)


def _emit_parent_child(dag_lines, parent_names, child_name, batch_size=200):
    # Batched into multiple PARENT/CHILD lines (rather than one huge line)
    # for readability/robustness -- a run with thousands of chunks would
    # otherwise produce one absurdly long line. Must stay >= the largest
    # single fan-in we ever emit (currently _MERGE_GROUP_SIZE=100 fill
    # parents per merge_group): splitting one child's parents across
    # multiple PARENT/CHILD lines triggers DAGMan's edge-adjustment pass,
    # which asserts ("GetStatus() != STATUS_DONE") when any of those split
    # parents is pre-marked DONE (see gaudi_jobs/calibration/README.md
    # troubleshooting note) -- observed on condor_dagman 24.12.17.
    for i in range(0, len(parent_names), batch_size):
        batch = parent_names[i:i + batch_size]
        dag_lines.append(f"PARENT {' '.join(batch)} CHILD {child_name}")


def main(argv=None):
    p = argparse.ArgumentParser(description="Generate the CALIBRATION (fill->merge->merge->fit) DAGMan pipeline.")
    p.add_argument("--runs", required=True,
                   help="Comma-separated run folder path(s) or bare run name(s), same format as "
                        "run_calibration_batch.py --runs. No auto-discovery by threshold.")
    p.add_argument("--raw-base", default=None, help="Base directory to resolve bare run names in --runs against.")
    p.add_argument("--th", default=None, help="Override the output th<N> label instead of auto-deriving it.")
    p.add_argument("--fill-scratch-dir", default=DEFAULT_FILL_SCRATCH_DIR,
                   help=f"EOS scratch area (decoded/ and hist/ live under it). Default: {DEFAULT_FILL_SCRATCH_DIR}")
    p.add_argument("--outdir", default=None,
                   help="Calibration table output base directory. Default: settings.yml calib_dir()/MuonCalib_gaudi")
    p.add_argument("--gain", choices=("high", "low", "both"), default="high")
    p.add_argument("--dag-dir", required=True,
                   help="Directory to write the .dag/.sub/.sh/file-lists/logs/ into.")
    p.add_argument("--diagnostics", action="store_true",
                    help="Also write a <output>.diagnostics.root cross-check file per Fit job.")
    p.add_argument("--fill-request-memory", type=int, default=5000,
                   help="Fill job request_memory in MB (default 5000: ~3.3GB for the 4 histogram grids + margin).")
    p.add_argument("--fill-job-flavour", default="microcentury")
    p.add_argument("--merge-request-memory", type=int, default=4000)
    p.add_argument("--merge-job-flavour", default="workday")
    p.add_argument("--fit-request-memory", type=int, default=4000)
    p.add_argument("--fit-job-flavour", default="workday")
    args = p.parse_args(argv)
    # Always absolute: a relative --dag-dir would leak into VARS filelist=
    # paths (chunklist_<run>.txt / chunklist_th<N>.txt) that merge.sh passes
    # straight to hadd on the Condor worker. Condor's own executable/log/
    # output/error resolution copes with a relative dag_dir fine (resolved
    # against the submission Iwd), but a plain `hadd @relative/path` inside
    # the job's shell depends on the worker's actual cwd matching the
    # submission directory -- which failed in practice for th230 ("could not
    # open indirect file"), forcing a manual local hadd. Absolute paths sidestep
    # the question entirely.
    args.dag_dir = os.path.abspath(args.dag_dir)

    kwargs = {"default_base": args.raw_base} if args.raw_base else {}
    run_folders = parse_run_folder_list(args.runs, **kwargs)

    decoded_base = decoded_dir(args.fill_scratch_dir)
    expected_chunks = _check_decoded_chunks(run_folders, decoded_base)

    th = _check_threshold_consistency(run_folders, args.th)
    label = _combined_label(run_folders, th)

    hist_base = hist_dir(args.fill_scratch_dir)
    out_base = args.outdir or os.path.join(paths.calib_dir(), "MuonCalib_gaudi")
    gains = ["high", "low"] if args.gain == "both" else [args.gain]

    os.makedirs(args.dag_dir, exist_ok=True)
    os.makedirs(os.path.join(args.dag_dir, "logs"), exist_ok=True)

    _fill_sh(args.dag_dir)
    _merge_sh(args.dag_dir)
    _fit_sh(args.dag_dir)
    _write_subs(args.dag_dir, args.fill_request_memory, args.fill_job_flavour, args.merge_request_memory,
                args.merge_job_flavour, args.fit_request_memory, args.fit_job_flavour)

    dag_lines = []
    merge_run_names = []
    for run_name, chunk_names in expected_chunks.items():
        fill_job_names = []
        chunk_paths = []
        run_hist_dir = os.path.join(hist_base, run_name)
        # Pre-create now: os.makedirs works on the EOS FUSE mount, shell
        # `mkdir -p` on the worker node does not (see condor_common.py).
        os.makedirs(run_hist_dir, exist_ok=True)
        for idx, chunk_name in enumerate(chunk_names):
            job_name = f"fill_{run_name}_{idx:04d}"
            in_decoded = os.path.join(decoded_base, run_name, chunk_name)
            out_hist = os.path.join(run_hist_dir, chunk_name)
            dag_lines.append(f"JOB {job_name} {args.dag_dir}/fill.sub")
            dag_lines.append(f'VARS {job_name} in_decoded="{in_decoded}" out_hist="{out_hist}"')
            dag_lines.append(f"RETRY {job_name} 2")
            dag_lines.append(f"CATEGORY {job_name} FillJobs")
            fill_job_names.append(job_name)
            chunk_paths.append(out_hist)

        merge_run_name = f"merge_run_{run_name}"
        merge_run_out = os.path.join(run_hist_dir, "merged_run.root")

        if len(chunk_paths) <= _MERGE_GROUP_SIZE:
            run_filelist = os.path.join(args.dag_dir, f"chunklist_{run_name}.txt")
            write_text(run_filelist, "\n".join(chunk_paths) + "\n")
            merge_run_parents = fill_job_names
        else:
            group_out_paths = []
            merge_run_parents = []
            for gidx, start in enumerate(range(0, len(chunk_paths), _MERGE_GROUP_SIZE)):
                group_paths = chunk_paths[start:start + _MERGE_GROUP_SIZE]
                group_fill_jobs = fill_job_names[start:start + _MERGE_GROUP_SIZE]
                group_name = f"merge_group_{run_name}_{gidx:04d}"
                group_out = os.path.join(run_hist_dir, f"merged_group_{gidx:04d}.root")
                group_filelist = os.path.join(args.dag_dir, f"chunklist_{run_name}_group{gidx:04d}.txt")
                write_text(group_filelist, "\n".join(group_paths) + "\n")
                dag_lines.append(f"JOB {group_name} {args.dag_dir}/merge.sub")
                dag_lines.append(f'VARS {group_name} filelist="{group_filelist}" out="{group_out}"')
                dag_lines.append(f"RETRY {group_name} 2")
                dag_lines.append(f"CATEGORY {group_name} MergeGroupJobs")
                _emit_parent_child(dag_lines, group_fill_jobs, group_name)
                group_out_paths.append(group_out)
                merge_run_parents.append(group_name)

            run_filelist = os.path.join(args.dag_dir, f"chunklist_{run_name}.txt")
            write_text(run_filelist, "\n".join(group_out_paths) + "\n")

        dag_lines.append(f"JOB {merge_run_name} {args.dag_dir}/merge.sub")
        dag_lines.append(f'VARS {merge_run_name} filelist="{run_filelist}" out="{merge_run_out}"')
        dag_lines.append(f"RETRY {merge_run_name} 2")
        _emit_parent_child(dag_lines, merge_run_parents, merge_run_name)
        merge_run_names.append(merge_run_name)

    merge_th_name = f"merge_threshold_th{th}"
    th_hist_dir = os.path.join(hist_base, f"th{th}")
    os.makedirs(th_hist_dir, exist_ok=True)
    merge_th_out = os.path.join(th_hist_dir, f"merged_th{th}.root")
    th_filelist = os.path.join(args.dag_dir, f"chunklist_th{th}.txt")
    write_text(th_filelist, "\n".join(
        os.path.join(hist_base, run_name, "merged_run.root") for run_name in expected_chunks) + "\n")
    dag_lines.append(f"JOB {merge_th_name} {args.dag_dir}/merge.sub")
    dag_lines.append(f'VARS {merge_th_name} filelist="{th_filelist}" out="{merge_th_out}"')
    dag_lines.append(f"RETRY {merge_th_name} 2")
    dag_lines.append(f"SCRIPT POST {merge_th_name} {_CHECK_SCRIPT} $RETURN {merge_th_out}")
    _emit_parent_child(dag_lines, merge_run_names, merge_th_name)

    fit_job_names = []
    os.makedirs(os.path.join(out_base, "pedestals", f"th{th}"), exist_ok=True)
    os.makedirs(os.path.join(out_base, "mips", f"th{th}"), exist_ok=True)
    for gain in gains:
        ped_out = os.path.join(out_base, "pedestals", f"th{th}", f"Pedestal_{label}_{gain}gain.txt")
        mip_out = os.path.join(out_base, "mips", f"th{th}", f"MIP_pedestalsubmode1_{label}_{gain}gain.txt")
        for mode, out_path in (("FitPedestal", ped_out), ("FitMip", mip_out)):
            job_name = f"fit_{mode.lower()}_{gain}"
            diag_path = os.path.splitext(out_path)[0] + ".diagnostics.root" if args.diagnostics else ""
            dag_lines.append(f"JOB {job_name} {args.dag_dir}/fit.sub")
            dag_lines.append(f'VARS {job_name} mode="{mode}" gain="{gain}" in_hist="{merge_th_out}" '
                              f'out_path="{out_path}" diag_path="{diag_path}"')
            fit_job_names.append(job_name)
    dag_lines.append(f"PARENT {merge_th_name} CHILD {' '.join(fit_job_names)}")

    dag_lines.append("MAXJOBS FillJobs 200")
    dag_lines.append("MAXJOBS MergeGroupJobs 100")

    dag_path = os.path.join(args.dag_dir, f"calibration_th{th}.dag")
    write_text(dag_path, "\n".join(dag_lines) + "\n")

    n_chunks = sum(len(c) for c in expected_chunks.values())
    print(f"[calibration] {len(run_folders)} run(s), {n_chunks} chunk(s) -> th{th}, label={label}")
    print(f"[calibration] Wrote {dag_path} (+ fill/merge/fit .sub/.sh, chunklists, logs/)")
    print(f"\nTo submit:\n  condor_submit_dag {dag_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
