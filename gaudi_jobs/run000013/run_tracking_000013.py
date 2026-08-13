#
# ACTS tracking para TB2026CERN_run_000013, sobre el EDM4hep PID ya producido
# por EcalToEDM4hep/EcalPidTransformer (mismo INPUT que usa steer_run000013.py,
# variante "_realvalues": event builder con calibración real).
#
# Uso (desde la raíz del repo):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/run000013/run_tracking_000013.py
#
import os

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, ShowerTagger, SiPadMeasConverter, ACTSProtoTracker
from Gaudi.Configuration import DEBUG, INFO

_RUN = "TB2026CERN_run_000013"
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SLAB_Z_FILE = os.path.join(REPO_ROOT, "mappings", "slab_z_positions.yml")

# EDM4hep PID file (ECalHits + ECalPid), same input steer_run000013.py's sibling
# examples read from -- override with TRACKING_INPUT_FILE.
_INPUT = os.environ.get(
    "TRACKING_INPUT_FILE",
    f"/eos/user/o/oarquero/TB2026CERN/data/{_RUN}/ecal_pid_{_RUN}_realvalues.root")
_OUTPUT = os.environ.get(
    "TRACKING_OUTPUT_FILE",
    os.path.join(os.path.dirname(_INPUT), f"tracks_{_RUN}.edm4hep.root"))

loglevel = DEBUG if os.environ.get("TRACKING_LOGLEVEL", "INFO") == "DEBUG" else INFO

print("Input File: ", _INPUT)
print("Output File:", _OUTPUT)

iosvc = IOSvc()
iosvc.Input          = _INPUT
iosvc.Output         = _OUTPUT
iosvc.outputCommands = ["keep *"]

geosvc = ACTSGeoSvc("ACTSGeoSvc")
geosvc.SlabZFile = SLAB_Z_FILE
geosvc.OutputLevel = loglevel

shower = ShowerTagger("ShowerTagger")
shower.InputCollection      = ["ECalHits"]
shower.OutputFlags          = ["SiPadShowerFlags"]
shower.OutputShowers        = ["EMShowers"]
shower.BitField             = "system:8,slab:8,chip:16,channel:8,sca:8"
shower.ShowerNHitsThreshold = 4
shower.ShowerMinConsecutive = 2
shower.MinTrackLayers       = 4
shower.OutputLevel          = loglevel

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = "ECalHits"
sipad_conv.OutputCollection = "SiPadMeasurements"
sipad_conv.InputFlags       = "SiPadShowerFlags"
sipad_conv.BitField         = "system:8,slab:8,chip:16,channel:8,sca:8"
sipad_conv.PixelSizeX       = 5.53
sipad_conv.PixelSizeY       = 5.53
sipad_conv.OutputLevel      = loglevel

proto = ACTSProtoTracker("ACTSProtoTracker")
proto.GeoSvc           = "ACTSGeoSvc"
proto.InputSiPad       = "SiPadMeasurements"
proto.OutputCollection = "ACTSTracks"
proto.AutoSeed         = True
proto.MaxSeeds         = 3
proto.HoughBinSize     = 5.0
proto.HoughHalfSize    = 90.0
proto.HoughMinVotes    = 3
proto.SeedCompatRadius = 8.0
proto.SeedMomentum     = float(os.environ.get("SEED_MOMENTUM", "3.0"))
proto.SeedStripPitch   = 5.53
proto.Chi2CutOff       = 70.0
proto.NumMeasCutOff    = 1
proto.MaxChi2PerNdf    = 10.0

NEIGHBOUR_WINDOW = 1.5 * 5.53
proto.HoughMaxMultiplicity  = 10.0
proto.IsolationWindow       = NEIGHBOUR_WINDOW
proto.IsolationMaxNeighbors = 2
proto.HitPurgeWindow        = NEIGHBOUR_WINDOW
proto.HitPurgeMaxNeighbors  = 4
proto.SeedCleaning             = True
proto.FinalRefit               = True
proto.DuplicateOverlapFraction = 0.7
proto.OutputLevel = loglevel

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = int(os.environ.get("TRACKING_EVT_MAX", "-1")),
    TopAlg  = [shower, sipad_conv, proto],
    ExtSvc  = [iosvc, geosvc]
)
