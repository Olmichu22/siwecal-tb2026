#
# ACTS tracking over the 15 SiW-ECAL slabs -- read real reconstructed hits
# (ECalHits, from EcalToEDM4hep/run_pid.py) and write ACTSTracks + EMShowers
# back into the SAME edm4hep file, in place.
#
# Mirrors siwecal_k4sim/gaudi_jobs/pid2026_common/job4_tracking.py: same
# component names, same property names, same tuned defaults, so the two
# repos' tracking output is directly comparable. The only real difference is
# geometry input -- this repo has no DD4hep compact XML, so ACTSGeoSvc builds
# its 15 plane surfaces from mappings/slab_z_positions.yml instead (see
# k4SiWEcalReco/ACTSGeoSvc.h and docs/acts_integration.md). Unlike simulation,
# real-data ECalHits are already in the frame ACTS surfaces expect (no
# DetectorFlipper-equivalent step), so tracking here needs no cross-file
# merge at all -- just a plain in-place update of ecal_<run>.edm4hep.root.
#
# A Gaudi IOSvc reader/writer cannot write back onto the file it is still
# reading, so this job's own default output is a HIDDEN TEMP FILE next to the
# input (must end in .root -- IOSvc picks its Writer backend off the
# filename). It is NOT a product: swap it onto the input yourself
#   (mv <output> <input>)
# for a bare `k4run` invocation, or use gaudi_jobs/run_tracking_batch.py,
# which drives this job AND does that swap automatically -- the normal way
# to run tracking.
#
#   TRACKING_INPUT_FILE=/path/ecal_pid_<label>.root \
#       k4run gaudi_source/options/run_tracking.py
#   mv /path/.<label>.tracking.tmp.root /path/ecal_pid_<label>.root
#
# Do NOT run this on a file that already has ACTSTracks/EMShowers/
# SiPadMeasurements/SiPadShowerFlags (e.g. re-running tracking on its own
# output): the components write those as new Gaudi collections, which
# collide with the ones already carried through by "keep *" and abort the
# process (not a graceful error). Track a file once.
#
import os
import sys
from pathlib import Path

from k4FWCore import ApplicationMgr, IOSvc
from Configurables import ACTSGeoSvc, ShowerTagger, SiPadMeasConverter, ACTSProtoTracker
from Gaudi.Configuration import DEBUG, INFO

REPO_ROOT = Path(__file__).resolve().parents[2]
SLAB_Z_FILE = REPO_ROOT / "mappings" / "slab_z_positions.yml"

infile = os.environ.get("TRACKING_INPUT_FILE", "")
if not infile:
    raise SystemExit("Set TRACKING_INPUT_FILE to an EDM4hep file with an "
                      "ECalHits collection (e.g. ecal_pid_<label>.root from "
                      "run_pid.py).")

# Default output: a hidden temp file next to the input -- see the module
# docstring above. Not a product; swap it onto the input yourself, or use
# run_tracking_batch.py, which does that automatically.
_stem = os.path.splitext(os.path.basename(infile))[0]
_label = _stem[len("ecal_pid_"):] if _stem.startswith("ecal_pid_") else _stem
outfile = os.environ.get(
    "TRACKING_OUTPUT_FILE",
    os.path.join(os.path.dirname(infile) or ".", f".{_label}.tracking.tmp.root"),
)
inputcoll = os.environ.get("TRACKING_INPUT_COLLECTION", "ECalHits")
seed_p    = float(os.environ.get("SEED_MOMENTUM", "3.0"))
loglevel  = DEBUG if os.environ.get("TRACKING_LOGLEVEL", "INFO") == "DEBUG" else INFO

print("Input File: ", infile)
print("Output File:", outfile)

iosvc = IOSvc()
iosvc.Input          = infile
iosvc.Output         = outfile
iosvc.outputCommands = ["keep *"]

# Builds the ACTS TrackingGeometry from mappings/slab_z_positions.yml: one
# plane surface per slab, at PadPitchMm/HalfSizeX/HalfSizeY -- see
# k4SiWEcalReco/ACTSGeoSvc.h for why this differs from siwecal_k4sim's
# DD4hep-based ACTSGeoSvc (same public interface either way, so
# ACTSProtoTracker below is unchanged from the simulation).
geosvc = ACTSGeoSvc("ACTSGeoSvc")
geosvc.SlabZFile = str(SLAB_Z_FILE)
geosvc.OutputLevel = loglevel

# ---- Shower tagging ------------------------------------------------------
# A track through an electromagnetic cascade is not a physical object: the
# pads are secondaries spraying transversely, not samples of a trajectory. So
# the shower hits are identified FIRST and never reach the ACTS measurement
# pool; what the event gets instead is a reconstructed shower (EMShowers).
#
# Every slab of this detector is 1.2-2.0 X0, so a high-energy electron is
# already showering in slab 0 or 1: with fewer than MinTrackLayers slabs
# before the onset there is no incoming segment to fit, and the whole event is
# vetoed for tracking. A muon never triggers the onset and is untouched.
shower = ShowerTagger("ShowerTagger")
shower.InputCollection      = [inputcoll]
shower.OutputFlags          = ["SiPadShowerFlags"]
shower.OutputShowers        = ["EMShowers"]
shower.BitField             = "system:8,slab:8,chip:16,channel:8,sca:8"
shower.ShowerNHitsThreshold = 4    # a MIP lights 1-2 pads per slab
shower.ShowerMinConsecutive = 2    # ignore a single-slab delta-ray spike
shower.MinTrackLayers       = 4
shower.OutputLevel          = loglevel

sipad_conv = SiPadMeasConverter("SiPadMeasConverter")
sipad_conv.GeoSvc           = "ACTSGeoSvc"
sipad_conv.InputCollection  = inputcoll
sipad_conv.OutputCollection = "SiPadMeasurements"
# Shower hits are dropped here, so they are never sampled by ACTS.
sipad_conv.InputFlags       = "SiPadShowerFlags"
sipad_conv.BitField         = "system:8,slab:8,chip:16,channel:8,sca:8"
sipad_conv.PixelSizeX       = 5.53
sipad_conv.PixelSizeY       = 5.53
sipad_conv.OutputLevel      = loglevel

proto = ACTSProtoTracker("ACTSProtoTracker")
proto.GeoSvc           = "ACTSGeoSvc"
proto.InputSiPad       = "SiPadMeasurements"
proto.OutputCollection = "ACTSTracks"

# ---- Seeding (Hough transform over the 2D pad hits) ------------------------
# Same tuned values as siwecal_k4sim's job4_tracking.py -- same pad pitch,
# same physical detector.
proto.AutoSeed         = True
proto.MaxSeeds         = 3
proto.HoughBinSize     = 5.0
proto.HoughHalfSize    = 90.0   # transverse envelope, matches ACTSGeoSvc.HalfSizeX/Y
proto.HoughMinVotes    = 3
proto.SeedCompatRadius = 8.0
proto.SeedMomentum     = seed_p
proto.SeedStripPitch   = 5.53

# ---- Track finding --------------------------------------------------------
proto.Chi2CutOff       = 70.0
proto.NumMeasCutOff    = 1
proto.MaxChi2PerNdf    = 10.0

# ---- Shower rejection ----------------------------------------------------
# Both neighbour-counting windows MUST be expressed in pad pitches: any window
# <= one pitch (5.53 mm) counts zero neighbours by construction (the nearest
# other pad on a plane is exactly one pitch away), which silently disables the
# filter. 1.5 pitches reaches the 8 pads surrounding a hit.
NEIGHBOUR_WINDOW = 1.5 * 5.53

proto.HoughMaxMultiplicity  = 10.0
proto.IsolationWindow       = NEIGHBOUR_WINDOW
proto.IsolationMaxNeighbors = 2
proto.HitPurgeWindow        = NEIGHBOUR_WINDOW
proto.HitPurgeMaxNeighbors  = 4

# ---- Track quality -------------------------------------------------------
proto.SeedCleaning            = True
proto.FinalRefit              = True
proto.DuplicateOverlapFraction = 0.7

proto.OutputLevel      = loglevel

ApplicationMgr(
    EvtSel  = "NONE",
    EvtMax  = -1,
    TopAlg  = [shower, sipad_conv, proto],
    ExtSvc  = [iosvc, geosvc]
)
