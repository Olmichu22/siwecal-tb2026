#!/usr/bin/env python
"""
Generic/batch version of gaudi_jobs/run000013/run_full_pipeline_000013.py:
raw2root -> event building -> PID/EDM4hep for any run, entirely in Gaudi/C++
(except the thin Python subprocess orchestration between the two k4run
processes -- see that script's module docstring for why two processes are
needed).

Calibration files are resolved via
siwecal_eventbuilder.cli.resolve_muon_calib_files(th) (the same lookup the
existing Python pipeline already uses), reused here rather than duplicated.

SAFETY: the raw data directory is READ-ONLY -- this script never writes
there; decoded/reconstructed outputs go to --converted-dir, and this script
never deletes anything.

Example::

    python gaudi_jobs/run_full_pipeline_batch.py --run TB2026CERN_run_000013 --th 230
"""
import argparse
import glob
import os
import subprocess
import sys

import ROOT

from siwecal_common import paths
from siwecal_eventbuilder.cli import resolve_muon_calib_files

_OPTIONS_DIR = os.path.join(os.path.dirname(__file__), "..", "gaudi_source", "options")
_RAW2ROOT_AND_EVBLD_STEERING = os.path.join(_OPTIONS_DIR, "run_raw2root_and_eventbuilder.py")
_PID_STEERING = os.path.join(_OPTIONS_DIR, "run_pid.py")

# Read-only source of raw DAQ chunk files (per-run subdirectories). NEVER
# write or delete anything here. Matches run_calibration_batch.py.
_DEFAULT_RAW_BASE = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata"
_DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test"


def _assert_not_under(path, readonly_root):
    real_path = os.path.realpath(path)
    real_root = os.path.realpath(readonly_root)
    if real_path == real_root or real_path.startswith(real_root + os.sep):
        raise SystemExit(f"REFUSING to write inside the read-only raw data directory: {path} is under {readonly_root}")


def main(argv=None):
    p = argparse.ArgumentParser(description="Full raw2root -> event building -> PID/EDM4hep pipeline for one run.")
    p.add_argument("--run", required=True, help="Run name (e.g. TB2026CERN_run_000013)")
    p.add_argument("--th", required=True, help="Threshold label to resolve calibration files (e.g. 230)")
    p.add_argument("--run-id", type=int, default=None, help="Numeric run id for the ecal tree's `run` branch "
                                                              "(default: parsed from the run name)")
    p.add_argument("--raw-dir", default=None,
                   help=f"Directory with <run>_raw.bin_NNNN chunk files (READ-ONLY). Default: {_DEFAULT_RAW_BASE}/<run>/")
    p.add_argument("--converted-dir", default=None,
                   help=f"Output base directory (separate from --raw-dir). Default: {_DEFAULT_CONVERTED_DIR}/<run>/")
    p.add_argument("--pedestal-file", default=None, help="Override the auto-resolved pedestal calibration file")
    p.add_argument("--mip-file", default=None, help="Override the auto-resolved MIP calibration file")
    p.add_argument("--no-calibration", action="store_true", help="Raw-ADC mode: pedestal=0, mip=1")
    p.add_argument("--padmap-default", default=None,
                   help="Default pad map file. Default: mappings/fev10_rotate_chip_channel_x_y_mapping.txt")
    p.add_argument("--padmap-slab12", default=None,
                   help="Slab-12 pad map override. Default: mappings/fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
    p.add_argument("--slab-z-file", default=None, help="Default: mappings/slab_z_positions.yml")
    p.add_argument("--hit-mip-cut", type=float, default=0.5, help="Hit-level MIP cut for the PID stage (default 0.5)")
    args = p.parse_args(argv)

    run = args.run
    raw_dir = args.raw_dir or os.path.join(_DEFAULT_RAW_BASE, run)
    converted_dir = args.converted_dir or os.path.join(_DEFAULT_CONVERTED_DIR, run)
    _assert_not_under(converted_dir, raw_dir)
    os.makedirs(converted_dir, exist_ok=True)

    siwecaldecoded_out = os.path.join(converted_dir, f"{run}.root")
    ecal_out = os.path.join(converted_dir, f"ecal_{run}.root")
    pid_out = os.path.join(converted_dir, f"ecal_{run}.edm4hep.root")

    mappings_dir = paths.geometry_dir()
    padmap_default = args.padmap_default or os.path.join(mappings_dir, "fev10_rotate_chip_channel_x_y_mapping.txt")
    padmap_slab12 = args.padmap_slab12 or os.path.join(mappings_dir, "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
    slab_z_file = args.slab_z_file or os.path.join(mappings_dir, "slab_z_positions.yml")

    if args.no_calibration:
        pedestal_file, mip_file = "", ""
    else:
        pedestal_file = args.pedestal_file
        mip_file = args.mip_file
        if not pedestal_file or not mip_file:
            auto_pedestal, auto_mip = resolve_muon_calib_files(args.th)
            pedestal_file = pedestal_file or auto_pedestal
            mip_file = mip_file or auto_mip

    run_id = args.run_id
    if run_id is None:
        digits = "".join(c for c in run.split("_")[-1] if c.isdigit())
        run_id = int(digits) if digits else -1

    # Stage 1: raw2root + event building, one k4run process.
    raw_files = sorted(glob.glob(os.path.join(raw_dir, f"{run}_raw.bin*")))
    if not raw_files:
        raise SystemExit(f"ERROR: no raw chunk files found in {raw_dir}")
    print(f"[stage1] {run}: {len(raw_files)} raw chunk file(s) -> {siwecaldecoded_out} -> {ecal_out}")
    env = {
        **os.environ,
        "RAW_FILES": ",".join(raw_files),
        "RAW2ROOT_OUT": siwecaldecoded_out,
        "RAW2ROOT_RUN_SETTINGS_FILE": os.path.join(raw_dir, "Run_Settings.txt"),
        "EVBLD_OUTPUT": ecal_out,
        "EVBLD_RUN_ID": str(run_id),
        "EVBLD_NO_CALIBRATION": "1" if args.no_calibration else "0",
        "EVBLD_PEDESTAL_FILE": pedestal_file,
        "EVBLD_MIP_FILE": mip_file,
        "EVBLD_PADMAP_DEFAULT": padmap_default,
        "EVBLD_PADMAP_SLAB_OVERRIDES": f"12:{padmap_slab12}",
        "EVBLD_SLAB_Z_FILE": slab_z_file,
    }
    result = subprocess.run(["k4run", _RAW2ROOT_AND_EVBLD_STEERING], env=env)
    if result.returncode != 0:
        raise SystemExit(f"ERROR: stage1 (raw2root + event building) k4run failed for {run}")

    # Stage 2: PID/EDM4hep -- separate k4run process, see module docstring.
    f = ROOT.TFile.Open(ecal_out)
    n_events = int(f.Get("ecal").GetEntries())
    f.Close()
    print(f"[stage2] {run}: {n_events} reconstructed event(s) -> {pid_out}")
    env = {
        **os.environ,
        "ECAL_FILE": ecal_out,
        "ECAL_TREE": "ecal",
        "ECAL_PID_OUT": pid_out,
        "ECAL_HIT_MIP_CUT": str(args.hit_mip_cut),
        "ECAL_MIP_THRESHOLDS": "",
    }
    result = subprocess.run(["k4run", _PID_STEERING], env=env)
    if result.returncode != 0:
        raise SystemExit(f"ERROR: stage2 (PID/EDM4hep) k4run failed for {run}")

    print(f"\n[Done] {pid_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
