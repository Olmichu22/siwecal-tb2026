#!/usr/bin/env python
"""
PIPELINE 2/2 -- event building + PID/EDM4hep + validation for
TB2026CERN_run_000012 (th230), continuing from the ALREADY-CONVERTED
siwecaldecoded file produced by pipeline 1 (convert_run000012.py).

Independent of pipeline 1: this script never touches the raw *_raw.bin*
chunks or runs EcalRawDecoder -- it only reads the siwecaldecoded ROOT file
convert_run000012.py already wrote, and errors out clearly if it isn't there
yet. Run pipeline 1 first (once); pipeline 2 can then be re-run on its own
(e.g. to try different calibration files) without re-converting.

Three steps, run as a plain Python driver (not a single k4run steering file)
because EcalToEDM4hep/EcalPidTransformer need EvtMax fixed before their k4run
process starts, but that count is only known after the event builder has run
-- same reasoning as gaudi_jobs/run000013/run_full_pipeline_000013.py:

1. ``options/run_event_builder.py``   (EcalEventBuilder):
   siwecaldecoded -> ecal_TB2026CERN_run_000012.root
2. ``options/run_pid.py``             (EcalToEDM4hep + EcalPidTransformer):
   ecal tree -> ecal_TB2026CERN_run_000012.edm4hep.root
3. ``python -m siwecal_validation``: validation plots from the PID output.

Usage::

    source setup.sh
    export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
    export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
    python gaudi_jobs/run000012/process_run000012.py
"""
import os
import subprocess
import sys

import ROOT

from siwecal_common import paths

_RUN = "TB2026CERN_run_000012"
_TH = "230"

_OPTIONS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "gaudi_source", "options")
_EVENT_BUILDER_STEERING = os.path.join(_OPTIONS_DIR, "run_event_builder.py")
_PID_STEERING = os.path.join(_OPTIONS_DIR, "run_pid.py")

_CONVERTED_DIR = os.path.join("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test", _RUN)
_SIWECALDECODED_IN = os.path.join(_CONVERTED_DIR, f"{_RUN}.root")  # produced by convert_run000012.py
_ECAL_OUT = os.path.join(_CONVERTED_DIR, f"ecal_{_RUN}.root")
_PID_OUT = os.path.join(_CONVERTED_DIR, f"ecal_{_RUN}.edm4hep.root")

# th230 calibration from TB2026CERN_run_000004 alone (the only muon/MIP
# calibration run at this threshold -- eudaq_255/eudaq_256 are electron runs
# and must not be used), MaxNhit=2/NSlabsHit=5, chip-fit + global MIP
# fallback, instead of the old MuonCalib_it2_corrected/run_000004 reference
# tables.
_CALIB_BASE = os.path.join(paths.calib_dir(), "MuonCalib_gaudi")
_PEDESTAL_FILE = os.path.join(_CALIB_BASE, "pedestals", f"th{_TH}", "Pedestal_TB2026CERN_run_000004_highgain.txt")
_MIP_FILE = os.path.join(_CALIB_BASE, "mips", f"th{_TH}", "MIP_pedestalsubmode1_TB2026CERN_run_000004_highgain.txt")

_MAPPINGS_DIR = paths.geometry_dir()
_PADMAP_DEFAULT = os.path.join(_MAPPINGS_DIR, "fev10_rotate_chip_channel_x_y_mapping.txt")
_PADMAP_SLAB12 = os.path.join(_MAPPINGS_DIR, "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt")
_SLAB_Z_FILE = os.path.join(_MAPPINGS_DIR, "slab_z_positions.yml")


def main():
    if not os.path.exists(_SIWECALDECODED_IN):
        raise SystemExit(
            f"ERROR: {_SIWECALDECODED_IN} not found.\n"
            f"Run pipeline 1 first: k4run gaudi_jobs/run000012/convert_run000012.py"
        )

    # Step 1: event building, one k4run process (EvtMax=1, all work in
    # EcalEventBuilder::initialize()).
    print(f"[step1] {_RUN}: {_SIWECALDECODED_IN} -> {_ECAL_OUT}")
    env = {
        **os.environ,
        "EVBLD_INPUT": _SIWECALDECODED_IN,
        "EVBLD_TREE": "siwecaldecoded",
        "EVBLD_OUTPUT": _ECAL_OUT,
        "EVBLD_RUN_ID": "12",
        "EVBLD_PEDESTAL_FILE": _PEDESTAL_FILE,
        "EVBLD_MIP_FILE": _MIP_FILE,
        "EVBLD_PADMAP_DEFAULT": _PADMAP_DEFAULT,
        "EVBLD_PADMAP_SLAB_OVERRIDES": f"12:{_PADMAP_SLAB12}",
        "EVBLD_SLAB_Z_FILE": _SLAB_Z_FILE,
    }
    result = subprocess.run(["k4run", _EVENT_BUILDER_STEERING], env=env)
    if result.returncode != 0:
        raise SystemExit("ERROR: step1 (event building) k4run failed")

    # Step 2: PID/EDM4hep. EvtMax must be known before k4run starts -- read
    # the ecal tree's entry count now, between the two k4run processes.
    f = ROOT.TFile.Open(_ECAL_OUT)
    n_events = int(f.Get("ecal").GetEntries())
    f.Close()
    print(f"[step2] {_RUN}: {n_events} reconstructed event(s) -> {_PID_OUT}")
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
        raise SystemExit("ERROR: step2 (PID/EDM4hep) k4run failed")

    # Step 3: validation plots, reading metrics from the PID output above
    # (found next to _ECAL_OUT since it's not under settings.yml's pid_dir).
    print(f"[step3] {_RUN}: validation plots from {_PID_OUT}")
    result = subprocess.run([
        sys.executable, "-m", "siwecal_validation",
        "--file", _ECAL_OUT, "--run", _RUN,
    ])
    if result.returncode != 0:
        raise SystemExit("ERROR: step3 (validation) failed")

    print(f"\n[Done] {_PID_OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
