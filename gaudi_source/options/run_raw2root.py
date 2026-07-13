#
# Gaudi/k4FWCore steering: raw SLB DAQ binary chunks -> siwecaldecoded ROOT
# tree (raw2root), via the EcalRawDecoder component.
#
#   RAW_FILES=/path/run_raw.bin_0000,/path/run_raw.bin_0001 \
#   RAW2ROOT_OUT=/path/run.root \
#       k4run gaudi_source/options/run_raw2root.py
#
# For runs with many chunks, RAW_FILES (an env var) can exceed the OS
# execve() argument/environment size limit ("Argument list too long") --
# use RAW_FILES_LIST instead, pointing at a text file with one chunk path
# per line:
#
#   RAW_FILES_LIST=/path/chunklist.txt RAW2ROOT_OUT=/path/run.root \
#       k4run gaudi_source/options/run_raw2root.py
#
# All work happens in EcalRawDecoder::initialize(); EvtMax=1 is enough since
# execute() is a no-op (see EcalRawDecoder.cpp for why this isn't a k4FWCore
# Producer: the output is a fixed-layout legacy ROOT tree, not an EDM4hep
# collection, and the number of decoded cycles isn't known until the whole
# input has been streamed).
import os

raw_files_list = os.environ.get("RAW_FILES_LIST", "")
if raw_files_list:
    with open(raw_files_list) as f:
        input_files = [line.strip() for line in f if line.strip()]
else:
    raw_files = os.environ.get("RAW_FILES", "")
    if not raw_files:
        raise SystemExit("Set RAW_FILES (comma list) or RAW_FILES_LIST (path to a newline-separated file)")
    input_files = [f for f in raw_files.split(",") if f.strip()]

# Decoding many chunks in ONE process silently drops acquisitions on some runs
# (run_000012: 30,020 of 119,211 survived; run_000013: untouched). The mechanism
# was never pinned down, so this is a warning, not a refusal -- but nothing in
# this repo should be hitting it: every driver decodes one chunk per process and
# lets the event builder chain the results. See gaudi_jobs/decode_chunks.py.
if len(input_files) > 1:
    print(f"WARNING: raw2root was handed {len(input_files)} raw files in one process. "
          f"This is the pattern that silently lost 75% of run_000012. Decode one chunk "
          f"per process instead (gaudi_jobs/decode_chunks.py), and check the result with "
          f"decode_chunks.assert_healthy().")

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalRawDecoder
from k4FWCore import ApplicationMgr

out_file = os.environ.get("RAW2ROOT_OUT", "")
if not out_file:
    raise SystemExit("Set RAW2ROOT_OUT to the output siwecaldecoded ROOT file path")

decoder = EcalRawDecoder(
    "EcalRawDecoder",
    InputFiles=input_files,
    OutputFile=out_file,
    TreeName=os.environ.get("RAW2ROOT_TREE", "siwecaldecoded"),
    MaxReadoutCycleJump=int(os.environ.get("RAW2ROOT_MAX_CYCLE_JUMP", "10")),
    BcidThreshold=int(os.environ.get("RAW2ROOT_BCID_THRESHOLD", "3")),
    ZeroSuppression=os.environ.get("RAW2ROOT_ZERO_SUPPRESSION", "0") == "1",
    EudaqFormat=os.environ.get("RAW2ROOT_EUDAQ", "0") == "1",
    ComputeBadBcid=os.environ.get("RAW2ROOT_COMPUTE_BADBCID", "1") == "1",
    ResetStatePerInputFile=os.environ.get("RAW2ROOT_RESET_PER_FILE", "1") == "1",
    RunSettingsFile=os.environ.get("RAW2ROOT_RUN_SETTINGS_FILE", ""),
)

ApplicationMgr(TopAlg=[decoder],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
