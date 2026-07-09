#
# Gaudi/k4FWCore steering: raw SLB DAQ binary chunks -> siwecaldecoded ->
# reconstructed `ecal` tree, in ONE k4run process (TopAlg=[EcalRawDecoder,
# EcalEventBuilder]).
#
# This combination is possible (unlike chaining all the way to
# EcalToEDM4hep/EcalPidTransformer in the same process) because BOTH
# EcalRawDecoder and EcalEventBuilder are plain Gaudi::Algorithm that do all
# their work in initialize() with EvtMax=1 -- neither needs to know an event
# count ahead of time. Gaudi runs TopAlg members' initialize() in order, so
# EcalRawDecoder has already written and closed its siwecaldecoded file by
# the time EcalEventBuilder::initialize() opens it.
#
# EcalToEDM4hep/EcalPidTransformer (k4FWCore Producer/Transformer) DO need
# EvtMax fixed before the process starts, which isn't known until
# EcalEventBuilder has run -- see gaudi_jobs/run000013/run_full_pipeline_000013.py
# for why that stage is a SEPARATE k4run invocation.
#
#   RAW_FILES=/path/run_raw.bin,/path/run_raw.bin_0001 \
#   RAW2ROOT_OUT=/path/run.root \
#   EVBLD_OUTPUT=/path/ecal_run.root \
#   EVBLD_PEDESTAL_FILE=/path/Pedestal_..._highgain.txt \
#   EVBLD_MIP_FILE=/path/MIP_pedestalsubmode1_..._highgain.txt \
#       k4run gaudi_source/options/run_raw2root_and_eventbuilder.py
#
import os

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalRawDecoder, EcalEventBuilder
from k4FWCore import ApplicationMgr

raw_files = os.environ.get("RAW_FILES", "")
if not raw_files:
    raise SystemExit("Set RAW_FILES to a comma-separated list of *_raw.bin_NNNN chunk files")
input_files = [f for f in raw_files.split(",") if f.strip()]

raw2root_out = os.environ.get("RAW2ROOT_OUT", "")
if not raw2root_out:
    raise SystemExit("Set RAW2ROOT_OUT to the intermediate siwecaldecoded ROOT file path")

evbld_out = os.environ.get("EVBLD_OUTPUT", "")
if not evbld_out:
    raise SystemExit("Set EVBLD_OUTPUT to the output ecal ROOT file path")

no_calibration = os.environ.get("EVBLD_NO_CALIBRATION", "0") == "1"
pedestal_file = os.environ.get("EVBLD_PEDESTAL_FILE", "")
mip_file = os.environ.get("EVBLD_MIP_FILE", "")
if not no_calibration and (not pedestal_file or not mip_file):
    raise SystemExit("Set EVBLD_PEDESTAL_FILE and EVBLD_MIP_FILE, or EVBLD_NO_CALIBRATION=1 for raw-ADC mode")

pad_map_overrides = [e for e in os.environ.get("EVBLD_PADMAP_SLAB_OVERRIDES", "").split(",") if e.strip()]

decoder = EcalRawDecoder(
    "EcalRawDecoder",
    InputFiles=input_files,
    OutputFile=raw2root_out,
    TreeName=os.environ.get("RAW2ROOT_TREE", "siwecaldecoded"),
    MaxReadoutCycleJump=int(os.environ.get("RAW2ROOT_MAX_CYCLE_JUMP", "10")),
    BcidThreshold=int(os.environ.get("RAW2ROOT_BCID_THRESHOLD", "3")),
    ZeroSuppression=os.environ.get("RAW2ROOT_ZERO_SUPPRESSION", "0") == "1",
    EudaqFormat=os.environ.get("RAW2ROOT_EUDAQ", "0") == "1",
    ComputeBadBcid=os.environ.get("RAW2ROOT_COMPUTE_BADBCID", "1") == "1",
    ResetStatePerInputFile=os.environ.get("RAW2ROOT_RESET_PER_FILE", "1") == "1",
    RunSettingsFile=os.environ.get("RAW2ROOT_RUN_SETTINGS_FILE", ""),
)

builder = EcalEventBuilder(
    "EcalEventBuilder",
    InputFile=raw2root_out,
    TreeName=os.environ.get("RAW2ROOT_TREE", "siwecaldecoded"),
    OutputFile=evbld_out,
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

ApplicationMgr(TopAlg=[decoder, builder],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
