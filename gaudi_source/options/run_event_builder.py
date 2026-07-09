#
# Gaudi/k4FWCore steering: siwecaldecoded ROOT tree -> reconstructed `ecal`
# tree (event building), via the EcalEventBuilder component.
#
#   EVBLD_INPUT=/path/run.root \
#   EVBLD_OUTPUT=/path/ecal_run.root \
#   EVBLD_PEDESTAL_FILE=/path/Pedestal_..._highgain.txt \
#   EVBLD_MIP_FILE=/path/MIP_pedestalsubmode1_..._highgain.txt \
#       k4run gaudi_source/options/run_event_builder.py
#
# threshold_dac is read directly from the siwecaldecoded tree's thresholdDac
# branch (written by EcalRawDecoder from Run_Settings.txt) -- this steering
# does not need to pass it explicitly unless EVBLD_THRESHOLD_DAC_OVERRIDE is
# set (e.g. for siwecaldecoded files predating that branch).
#
# All work happens in EcalEventBuilder::initialize(); EvtMax=1 is enough
# since execute() is a no-op (this is a global reduction over an entire run,
# not a per-event k4FWCore transform -- see EcalEventBuilder.cpp).
import os

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalEventBuilder
from k4FWCore import ApplicationMgr

input_file = os.environ.get("EVBLD_INPUT", "")
if not input_file:
    raise SystemExit("Set EVBLD_INPUT to the input siwecaldecoded ROOT file")

output_file = os.environ.get("EVBLD_OUTPUT", "")
if not output_file:
    raise SystemExit("Set EVBLD_OUTPUT to the output ecal ROOT file path")

no_calibration = os.environ.get("EVBLD_NO_CALIBRATION", "0") == "1"
pedestal_file = os.environ.get("EVBLD_PEDESTAL_FILE", "")
mip_file = os.environ.get("EVBLD_MIP_FILE", "")
if not no_calibration and (not pedestal_file or not mip_file):
    raise SystemExit("Set EVBLD_PEDESTAL_FILE and EVBLD_MIP_FILE, or EVBLD_NO_CALIBRATION=1 for raw-ADC mode")

pad_map_overrides = [e for e in os.environ.get("EVBLD_PADMAP_SLAB_OVERRIDES", "").split(",") if e.strip()]

builder = EcalEventBuilder(
    "EcalEventBuilder",
    InputFile=input_file,
    TreeName=os.environ.get("EVBLD_TREE", "siwecaldecoded"),
    OutputFile=output_file,
    RunId=int(os.environ.get("EVBLD_RUN_ID", "-1")),
    ThresholdDacOverride=int(os.environ.get("EVBLD_THRESHOLD_DAC_OVERRIDE", "-1")),
    NoCalibration=no_calibration,
    PedestalFile=pedestal_file,
    MipFile=mip_file,
    PadMapDefaultFile=os.environ.get("EVBLD_PADMAP_DEFAULT", ""),
    PadMapSlabOverrides=pad_map_overrides,
    SlabZFile=os.environ.get("EVBLD_SLAB_Z_FILE", ""),
    NoMapping=os.environ.get("EVBLD_NO_MAPPING", "0") == "1",
)

ApplicationMgr(TopAlg=[builder],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
