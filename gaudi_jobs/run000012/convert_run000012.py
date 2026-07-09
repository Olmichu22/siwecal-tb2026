#
# PIPELINE 1/2 -- raw2root ONLY for TB2026CERN_run_000012 (th230).
#
# Converts the raw SLB DAQ binary chunks into the intermediate
# `siwecaldecoded` ROOT tree, via EcalRawDecoder. Independent of pipeline 2
# (process_run000012.py): this script never runs the event builder, PID, or
# validation -- it only produces the converted file that pipeline 2 reads.
#
# Uso (desde la raiz del repo):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/run000012/convert_run000012.py
#
import glob
import os

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalRawDecoder
from k4FWCore import ApplicationMgr

_RUN = "TB2026CERN_run_000012"

_RAW_DIR = os.path.join("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata", _RUN)  # READ-ONLY
_CONVERTED_DIR = os.path.join("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test", _RUN)
_SIWECALDECODED_OUT = os.path.join(_CONVERTED_DIR, f"{_RUN}.root")

os.makedirs(_CONVERTED_DIR, exist_ok=True)

raw_files = sorted(glob.glob(os.path.join(_RAW_DIR, f"{_RUN}_raw.bin*")))
if not raw_files:
    raise SystemExit(f"ERROR: no raw chunk files found in {_RAW_DIR}")
print(f"[convert] {_RUN}: {len(raw_files)} raw chunk file(s) -> {_SIWECALDECODED_OUT}")

decoder = EcalRawDecoder(
    "EcalRawDecoder",
    InputFiles=raw_files,
    OutputFile=_SIWECALDECODED_OUT,
    TreeName="siwecaldecoded",
    MaxReadoutCycleJump=10,
    BcidThreshold=3,
    ZeroSuppression=False,
    EudaqFormat=False,
    ComputeBadBcid=True,
    ResetStatePerInputFile=True,
    RunSettingsFile=os.path.join(_RAW_DIR, "Run_Settings.txt"),
)

# All work happens in EcalRawDecoder::initialize(); EvtMax=1 is enough since
# execute() is a no-op (see EcalRawDecoder.cpp).
ApplicationMgr(TopAlg=[decoder],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
