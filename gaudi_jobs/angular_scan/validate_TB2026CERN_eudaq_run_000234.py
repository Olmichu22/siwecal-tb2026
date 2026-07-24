#
# Gaudi/k4FWCore steering for TB2026CERN_eudaq_run_000234 (th210): the VALIDATION / EDM4hep stage.
# Reads the event-built `ecal` tree produced by steer_TB2026CERN_eudaq_run_000234.py and writes the
# EDM4hep + PID output next to it, in <run>/events/ -- NOT in the shared
# Reconstruction/ area.
#
# Run steer_TB2026CERN_eudaq_run_000234.py first (it makes <run>/events/ecal_<run>.root).
#
# Part of the angular_scan pipeline (gaudi_jobs/angular_scan/README.md).
#
# Usage (from the repo root):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/angular_scan/validate_TB2026CERN_eudaq_run_000234.py
#
import os

import ROOT
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalToEDM4hep, EcalPidTransformer
from k4FWCore import ApplicationMgr, IOSvc

_RUN = "TB2026CERN_eudaq_run_000234"
_CONVERTED = os.environ.get(
    "EVBLD_CONVERTED_DIR",
    "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi")
_EVENTS = os.path.join(_CONVERTED, _RUN, "events")
_INPUT = os.environ.get("VALIDATE_INPUT", os.path.join(_EVENTS, f"ecal_{_RUN}.root"))
_OUTPUT = os.environ.get("VALIDATE_OUTPUT", os.path.join(_EVENTS, f"ecal_{_RUN}.edm4hep.root"))
os.makedirs(os.path.dirname(_OUTPUT), exist_ok=True)

_f = ROOT.TFile.Open(_INPUT)
if not _f or _f.IsZombie():
    raise SystemExit(f"input not found (run steer_TB2026CERN_eudaq_run_000234.py first): {_INPUT}")
n_events = int(_f.Get("ecal").GetEntries())
_f.Close()

svc = IOSvc("IOSvc")
svc.Output = _OUTPUT

source = EcalToEDM4hep("EcalToEDM4hep",
                        InputFile=_INPUT,
                        TreeName="ecal",
                        HitMipCut=0.5)            # <0 to disable the per-hit MIP cut

pid = EcalPidTransformer("EcalPidTransformer",
                          InputCaloHits=["ECalHits"],
                          OutputClusters=["ECalPid"],
                          MipThresholds=[0.5, 1.0])  # validation mode: mip05_/mip1_ shape blocks

ApplicationMgr(TopAlg=[source, pid],
               EvtSel="NONE",
               EvtMax=n_events,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
