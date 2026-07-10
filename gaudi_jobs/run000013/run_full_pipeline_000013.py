#!/usr/bin/env python
"""
Full raw-to-EDM4hep pipeline example for TB2026CERN_run_000013 (52 GeV
electrons, beam point P1_52GeV -- see configs/data/data_reference_base.yml),
entirely in Gaudi/C++: raw2root -> event building -> PID/EDM4hep.

This is a plain Python script (run with ``python``, NOT ``k4run`` -- unlike
``steer_run000013.py``, which IS a k4run steering file), because the full
chain needs two separate k4run processes:

1. ``options/run_raw2root_and_eventbuilder.py``
   (``TopAlg=[EcalRawDecoder, EcalEventBuilder]``): raw binary chunks ->
   siwecaldecoded -> ``ecal_TB2026CERN_run_000013.root``.
2. ``options/run_pid.py`` (``TopAlg=[EcalToEDM4hep, EcalPidTransformer]``):
   the reconstructed ``ecal`` tree -> EDM4hep.

These can't be one single TopAlg/k4run process: EcalToEDM4hep/
EcalPidTransformer are k4FWCore Producer/Transformer components that need
``EvtMax`` (the number of ``ecal`` tree entries) fixed before the process
starts, but that count is only known after EcalEventBuilder has run --
inside the SAME process for a combined TopAlg. Splitting into two k4run
calls, with this script reading the entry count in between, resolves it (the
same pattern already used by gaudi_jobs/run_calibration_batch.py's two-stage
chaining and gaudi_jobs/run_pid_batch.py).

SAFETY: the raw data directory is READ-ONLY -- this script never writes
there. The pedestal/MIP calibration tables are resolved from MuonCalib_gaudi
for the run's own ThresholdDAC (read from its Run_Settings.txt), and come from
a DIFFERENT data-taking run than the electron physics run being converted
(000013), matching how pedestal/MIP calibration is generated in practice (see
gaudi_jobs/run_calibration_batch.py).

Usage::

    source setup.sh   # or: source /cvmfs/sw.hsf.org/key4hep/setup.sh -r <release>
    export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
    export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
    python gaudi_jobs/run000013/run_full_pipeline_000013.py
"""
import glob
import os
import subprocess
import sys

import ROOT

from siwecal_common import paths
from siwecal_eventbuilder.cli import resolve_gaudi_calib_files
from siwecal_eventbuilder.run_settings import read_threshold_dac

_RUN = "TB2026CERN_run_000013"

_OPTIONS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "gaudi_source", "options")
_RAW2ROOT_AND_EVBLD_STEERING = os.path.join(_OPTIONS_DIR, "run_raw2root_and_eventbuilder.py")
_PID_STEERING = os.path.join(_OPTIONS_DIR, "run_pid.py")

_RAW_DIR = os.path.join("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata", _RUN)  # READ-ONLY
_RUN_SETTINGS_FILE = os.path.join(_RAW_DIR, "Run_Settings.txt")
_CONVERTED_DIR = os.path.join("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test", _RUN)
_SIWECALDECODED_OUT = os.path.join(_CONVERTED_DIR, f"{_RUN}.root")
_ECAL_OUT = os.path.join(_CONVERTED_DIR, f"ecal_{_RUN}.root")
_PID_OUT = os.path.join(_CONVERTED_DIR, f"ecal_{_RUN}.edm4hep.root")

_MAPPINGS_DIR = paths.geometry_dir()
_PADMAP_DEFAULT = os.path.join(_MAPPINGS_DIR, "fev10_rotate_chip_channel_x_y_mapping.txt")
_PADMAP_SLAB12 = os.path.join(_MAPPINGS_DIR, "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
_SLAB_Z_FILE = os.path.join(_MAPPINGS_DIR, "slab_z_positions.yml")


def main():
    os.makedirs(_CONVERTED_DIR, exist_ok=True)

    threshold_dac = read_threshold_dac(_RUN_SETTINGS_FILE)
    if threshold_dac < 0:
        raise SystemExit(f"ERROR: cannot read ThresholdDAC from {_RUN_SETTINGS_FILE}")
    pedestal_file, mip_file = resolve_gaudi_calib_files(threshold_dac)
    print(f"[calib] {_RUN}: ThresholdDAC={threshold_dac} (from {_RUN_SETTINGS_FILE})")
    print(f"[calib] {_RUN}: pedestal={pedestal_file}")
    print(f"[calib] {_RUN}: mip={mip_file}")

    # Stage 1: raw2root + event building, one k4run process.
    raw_files = sorted(glob.glob(os.path.join(_RAW_DIR, f"{_RUN}_raw.bin*")))
    if not raw_files:
        raise SystemExit(f"ERROR: no raw chunk files found in {_RAW_DIR}")
    print(f"[stage1] {_RUN}: {len(raw_files)} raw chunk file(s) -> {_SIWECALDECODED_OUT} -> {_ECAL_OUT}")
    env = {
        **os.environ,
        "RAW_FILES": ",".join(raw_files),
        "RAW2ROOT_OUT": _SIWECALDECODED_OUT,
        "RAW2ROOT_RUN_SETTINGS_FILE": _RUN_SETTINGS_FILE,
        "EVBLD_OUTPUT": _ECAL_OUT,
        "EVBLD_RUN_ID": "13",
        "EVBLD_PEDESTAL_FILE": pedestal_file,
        "EVBLD_MIP_FILE": mip_file,
        "EVBLD_PADMAP_DEFAULT": _PADMAP_DEFAULT,
        "EVBLD_PADMAP_SLAB_OVERRIDES": f"12:{_PADMAP_SLAB12}",
        "EVBLD_SLAB_Z_FILE": _SLAB_Z_FILE,
    }
    result = subprocess.run(["k4run", _RAW2ROOT_AND_EVBLD_STEERING], env=env)
    if result.returncode != 0:
        raise SystemExit("ERROR: stage1 (raw2root + event building) k4run failed")

    # Stage 2: PID/EDM4hep. EvtMax must be known before k4run starts -- read
    # the ecal tree's entry count now, in this driver script, between the two
    # k4run processes (see module docstring).
    f = ROOT.TFile.Open(_ECAL_OUT)
    n_events = int(f.Get("ecal").GetEntries())
    f.Close()
    print(f"[stage2] {_RUN}: {n_events} reconstructed event(s) -> {_PID_OUT}")
    env = {
        **os.environ,
        "ECAL_FILE": _ECAL_OUT,
        "ECAL_TREE": "ecal",
        "ECAL_PID_OUT": _PID_OUT,
        "ECAL_HIT_MIP_CUT": "0.5",
        "ECAL_MIP_THRESHOLDS": "0.5,1.0",
    }
    result = subprocess.run(["k4run", _PID_STEERING], env=env)
    if result.returncode != 0:
        raise SystemExit("ERROR: stage2 (PID/EDM4hep) k4run failed")

    print(f"\n[Done] {_PID_OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
