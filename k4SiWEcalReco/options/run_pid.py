#
# Gaudi/k4FWCore steering: ecal tree -> EDM4hep (CalorimeterHit + per-event
# Cluster with the particle-ID shower variables in shapeParameters).
#
#   ECAL_FILE=/path/ecal_run.root ECAL_PID_OUT=out.root \
#       k4run k4SiWEcalReco/options/run_pid.py
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

source = EcalToEDM4hep("EcalToEDM4hep", InputFile=ecal_file, TreeName=tree_name)
pid = EcalPidTransformer("EcalPidTransformer",
                         InputCaloHits=["ECalHits"], OutputClusters=["ECalPid"])

ApplicationMgr(TopAlg=[source, pid],
               EvtSel="NONE",
               EvtMax=n_events,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
