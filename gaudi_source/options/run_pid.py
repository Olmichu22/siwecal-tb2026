#
# Gaudi/k4FWCore steering: ecal tree -> EDM4hep (CalorimeterHit + per-event
# Cluster with the particle-ID shower variables in shapeParameters).
#
#   ECAL_FILE=/path/ecal_run.root ECAL_PID_OUT=out.root \
#       k4run gaudi_source/options/run_pid.py
#
import os

import ROOT
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalToEDM4hep, EcalPidTransformer
from k4FWCore import ApplicationMgr, IOSvc

ecal_file = os.environ.get("ECAL_FILE", "")
tree_name = os.environ.get("ECAL_TREE", "ecal")
if not ecal_file:
    raise SystemExit("Set ECAL_FILE to the input ecal_<run>.root")

# Default output: next to the input, with the same naming run_pid_batch uses
# (ecal_<run>.root -> ecal_pid_<run>.root).
_stem = os.path.splitext(os.path.basename(ecal_file))[0]
_label = _stem[len("ecal_"):] if _stem.startswith("ecal_") else _stem
out_file = os.environ.get(
    "ECAL_PID_OUT",
    os.path.join(os.path.dirname(ecal_file), f"ecal_pid_{_label}.root"),
)
print("Output File:", out_file)
# The ecal tree is one row per physics event -> a strict 1->1 source. The event
# count is just GetEntries (no fan-out, no pre-scan needed).
_f = ROOT.TFile.Open(ecal_file)
n_events = int(_f.Get(tree_name).GetEntries())
_f.Close()

svc = IOSvc("IOSvc")
svc.Output = out_file
# No svc.Input: EcalToEDM4hep reads the (non-podio) ecal tree itself.

# Run mode. These two are DIFFERENT knobs (see the main README, "Two different
# MIP cut knobs"):
#  - ECAL_HIT_MIP_CUT: per-hit energy cut -- drop hits below this many MIPs
#    (<0 disables). Default off here; run_full_pipeline_batch.py sets 0.5.
#  - ECAL_MIP_THRESHOLDS: comma list of EXTRA thresholds at which the whole
#    shower-variable block is recomputed and stored again as mip05_/mip1_
#    ("0.5,1.0" = visualiser mode; "" = physics mode, one block only).
#    NOTE the default below ("0.5,1.0") applies ONLY to a bare `k4run
#    run_pid.py`. The batch/interactive drivers (run_pid_batch.py without
#    --validation, run_full_pipeline_batch.py) override it to "" -- physics
#    mode -- so those pipelines write one block, not three. This is intentional.
hit_mip_cut = float(os.environ.get("ECAL_HIT_MIP_CUT", "-1"))
_raw_mip = os.environ.get("ECAL_MIP_THRESHOLDS", "0.5,1.0")
mip_thresholds = [float(t) for t in _raw_mip.split(",") if t.strip()]

source = EcalToEDM4hep("EcalToEDM4hep", InputFile=ecal_file, TreeName=tree_name,
                       HitMipCut=hit_mip_cut)
pid = EcalPidTransformer("EcalPidTransformer",
                         InputCaloHits=["ECalHits"], OutputClusters=["ECalPid"],
                         MipThresholds=mip_thresholds)
# is_shower onset-criterion overrides (optional; the component defaults --
# ShowerCoreRadiusMm=30.0, ShowerOnsetMinNhit=4, ShowerOnsetMinConsecutive=2 --
# apply unless set here, same pattern as ECAL_HIT_MIP_CUT above).
if os.environ.get("ECAL_SHOWER_CORE_RADIUS_MM"):
    pid.ShowerCoreRadiusMm = float(os.environ["ECAL_SHOWER_CORE_RADIUS_MM"])
if os.environ.get("ECAL_SHOWER_ONSET_MIN_NHIT"):
    pid.ShowerOnsetMinNhit = int(os.environ["ECAL_SHOWER_ONSET_MIN_NHIT"])
if os.environ.get("ECAL_SHOWER_ONSET_MIN_CONSECUTIVE"):
    pid.ShowerOnsetMinConsecutive = int(os.environ["ECAL_SHOWER_ONSET_MIN_CONSECUTIVE"])

ApplicationMgr(TopAlg=[source, pid],
               EvtSel="NONE",
               EvtMax=n_events,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
