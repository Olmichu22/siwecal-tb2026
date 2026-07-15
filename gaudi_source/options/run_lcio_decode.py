#
# Gaudi/k4FWCore steering: EUDAQ LCIO (ROC_run_<N>_tp.slcio) -> siwecaldecoded
# ROOT tree, via the EcalLcioDecoder component. The LCIO twin of run_raw2root.py:
# same output tree, so everything downstream reads it unchanged.
#
#   LCIO_IN=/eos/.../Data/eudaq/ROC_run_000147_tp.slcio \
#   LCIO_OUT=/path/run.root \
#   LCIO_RUN_SETTINGS_FILE=/eos/.../rundata/TB2026CERN_eudaq_run_000147/Run_Settings.txt \
#       k4run gaudi_source/options/run_lcio_decode.py
#
# All work happens in EcalLcioDecoder::initialize(); EvtMax=1, execute() is a
# no-op -- same reasoning as EcalRawDecoder.
import os

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalLcioDecoder
from k4FWCore import ApplicationMgr

in_file = os.environ.get("LCIO_IN", "")
if not in_file:
    raise SystemExit("Set LCIO_IN to the input .slcio file path")
out_file = os.environ.get("LCIO_OUT", "")
if not out_file:
    raise SystemExit("Set LCIO_OUT to the output siwecaldecoded ROOT file path")

decoder = EcalLcioDecoder(
    "EcalLcioDecoder",
    InputFile=in_file,
    OutputFile=out_file,
    TreeName=os.environ.get("LCIO_TREE", "siwecaldecoded"),
    CollectionName=os.environ.get("LCIO_COLLECTION", "EUDAQDataSiECAL"),
    RunSettingsFile=os.environ.get("LCIO_RUN_SETTINGS_FILE", ""),
    MaxAcq=int(os.environ.get("LCIO_MAX_ACQ", "-1")),
)

ApplicationMgr(TopAlg=[decoder],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
