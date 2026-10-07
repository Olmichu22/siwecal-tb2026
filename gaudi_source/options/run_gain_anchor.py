#
# Gaudi/k4FWCore steering: siwecaldecoded chunks of an electron run -> LG->HG anchor file, via the
# GainAnchorCalibrator component.
#
#   ANCHOR_INPUT_FILES=/path/chunks/chunk_0000.root,/path/chunks/chunk_0001.root \
#   ANCHOR_PEDESTAL_FILE=calibration/MuonCalib_gaudi_fixed/pedestals/th210/Pedestal_TB2026CERN_run_000th210_highgain.txt \
#   ANCHOR_PEDESTAL_FILE_LOWGAIN=calibration/MuonCalib_gaudi_fixed/pedestals/th210/Pedestal_TB2026CERN_run_000th210_lowgain.txt \
#   ANCHOR_OUTPUT_FILE=/path/anchor/th210/gain_anchor_th210.txt \
#       k4run gaudi_source/options/run_gain_anchor.py
#
# ANCHOR_INPUT_FILES also takes a glob (sorted). Optional:
#   ANCHOR_DIAGNOSTICS_FILE   ROOT file with the 2D histogram, slice medians, fitted points, k and c
#   ANCHOR_MAX_ACQ            acquisitions read (default 8000, 0 = all)
#   ANCHOR_FID_LO / _HI       fiducial band on pedestal-subtracted adc_high (default 150 / 1200)
#   ANCHOR_MIN_SLICE_ENTRIES  (default 50)
#   ANCHOR_MIN_HITBITS / ANCHOR_MAX_HITBITS   occupancy window, hit bits per chip-SCA (default 0 / 64 = all)
#   ANCHOR_SWITCH             written to the file as 'switch' (default 1500)
#   ANCHOR_SOURCE / ANCHOR_TH  run name and th label for the file header
#
# All work happens in GainAnchorCalibrator::initialize(); execute() is a no-op, EvtMax=1.
# The batch driver is gaudi_jobs/calibration/run_anchor_calibration.py.
import glob
import os

import Gaudi.Configuration as _gaudi_cfg
from Configurables import EventDataSvc
from Configurables import GainAnchorCalibrator
from k4FWCore import ApplicationMgr

_output_level = getattr(_gaudi_cfg, os.environ.get("ANCHOR_OUTPUT_LEVEL", "INFO"), _gaudi_cfg.INFO)

raw = os.environ.get("ANCHOR_INPUT_FILES", "")
if not raw:
    raise SystemExit("Set ANCHOR_INPUT_FILES (comma-separated siwecaldecoded chunks, or a glob)")
if any(ch in raw for ch in "*?["):
    input_files = sorted(glob.glob(raw))
    if not input_files:
        raise SystemExit(f"ANCHOR_INPUT_FILES glob matched nothing: {raw}")
else:
    input_files = [f for f in raw.split(",") if f.strip()]

for var in ("ANCHOR_PEDESTAL_FILE", "ANCHOR_PEDESTAL_FILE_LOWGAIN", "ANCHOR_OUTPUT_FILE"):
    if not os.environ.get(var):
        raise SystemExit(f"Set {var}")

anchor = GainAnchorCalibrator(
    "GainAnchorCalibrator",
    InputFiles=input_files,
    TreeName=os.environ.get("ANCHOR_TREE", "siwecaldecoded"),
    PedestalFile=os.environ["ANCHOR_PEDESTAL_FILE"],
    PedestalFileLowGain=os.environ["ANCHOR_PEDESTAL_FILE_LOWGAIN"],
    OutputFile=os.environ["ANCHOR_OUTPUT_FILE"],
    DiagnosticsFile=os.environ.get("ANCHOR_DIAGNOSTICS_FILE", ""),
    MaxAcquisitions=int(os.environ.get("ANCHOR_MAX_ACQ", "8000")),
    FiducialLow=float(os.environ.get("ANCHOR_FID_LO", "150")),
    FiducialHigh=float(os.environ.get("ANCHOR_FID_HI", "1200")),
    MinEntriesPerSlice=int(os.environ.get("ANCHOR_MIN_SLICE_ENTRIES", "50")),
    MinHitBits=int(os.environ.get("ANCHOR_MIN_HITBITS", "0")),
    MaxHitBits=int(os.environ.get("ANCHOR_MAX_HITBITS", "64")),
    SwitchAdc=float(os.environ.get("ANCHOR_SWITCH", "1500")),
    SourceLabel=os.environ.get("ANCHOR_SOURCE", ""),
    ThresholdLabel=os.environ.get("ANCHOR_TH", ""),
    OutputLevel=_output_level,
)

ApplicationMgr(
    TopAlg=[anchor],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=_output_level,
)
