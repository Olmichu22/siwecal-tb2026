#
# LG->HG anchor of the th210 period, measured on TB2026CERN_eudaq_run_000166 (99 GeV electrons) with
# GainAnchorCalibrator -- the run and settings the deployed th210 anchor
# (/eos/experiment/drdcalo/siw-ecal/TB2026-06/Calib/MuonCalib/anchor/th210/gain_anchor_th210.txt) is produced with.
#
# Run it (from the repo root):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/run000166/steer_anchor_run000166.py
#
# Output, in this folder unless overridden: gain_anchor_th210.txt (k, c, switch, source_run, statistics) and
# gain_anchor_fit_000166_th210.root (2D histogram, slice medians, fitted points). Expected: k = 0.09243,
# c = 2.345 on 6,489,694 fired channels.
# Overrides: ANCHOR_OUTPUT_FILE, ANCHOR_DIAGNOSTICS_FILE, ANCHOR_MIN_HITBITS / ANCHOR_MAX_HITBITS (occupancy window).
#
# All three periods, written where the reco DAG reads them: gaudi_jobs/calibration/run_anchor_calibration.py
#
import glob
import os

from Configurables import EventDataSvc, GainAnchorCalibrator
from k4FWCore import ApplicationMgr
from Gaudi.Configuration import INFO

_RUN = "TB2026CERN_eudaq_run_000166"
_HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(_HERE))
CHUNKS = sorted(glob.glob(
    f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi/{_RUN}/chunks/chunk_*.root"))[:60]
if not CHUNKS:
    raise SystemExit(f"No decoded chunks for {_RUN}: python gaudi_jobs/decode_lcio_runs.py --runs 166")
PED = os.path.join("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Calib/MuonCalib", "pedestals", "th210")

anchor = GainAnchorCalibrator(
    "GainAnchorCalibrator",
    InputFiles=CHUNKS,
    PedestalFile=os.path.join(PED, "Pedestal_TB2026CERN_run_000th210_highgain.txt"),
    PedestalFileLowGain=os.path.join(PED, "Pedestal_TB2026CERN_run_000th210_lowgain.txt"),
    OutputFile=os.environ.get("ANCHOR_OUTPUT_FILE", os.path.join(_HERE, "gain_anchor_th210.txt")),
    DiagnosticsFile=os.environ.get("ANCHOR_DIAGNOSTICS_FILE",
                                   os.path.join(_HERE, "gain_anchor_fit_000166_th210.root")),
    MaxAcquisitions=8000,
    FiducialLow=150.0,
    FiducialHigh=1200.0,
    MinHitBits=int(os.environ.get("ANCHOR_MIN_HITBITS", "0")),
    MaxHitBits=int(os.environ.get("ANCHOR_MAX_HITBITS", "64")),
    SourceLabel=_RUN,
    ThresholdLabel="210",
    OutputLevel=INFO,
)

ApplicationMgr(TopAlg=[anchor], EvtSel="NONE", EvtMax=1, ExtSvc=[EventDataSvc("EventDataSvc")], OutputLevel=INFO)
