#
# Gaudi/k4FWCore steering for TB2026CERN_eudaq_run_000167 (th210, electrons --
# the TEST run twinned with run000166): the VALIDATION / EDM4hep stage. Reads the
# event-built `ecal` tree in <run>/events/ and writes the EDM4hep + PID output
# next to it -- NOT in the shared Reconstruction/ area.
#
# The event-built ecal is produced the same way as run000166 (steer with the
# th210 calibration from calibration/MuonCalib_gaudi/); this stage only adds the
# EDM4hep + PID pass on top of it.
#
# Usage (from the repo root):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/run000167/validate_run000167.py
#
import os

import ROOT
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalToEDM4hep, EcalPidTransformer
from k4FWCore import ApplicationMgr, IOSvc

_RUN = "TB2026CERN_eudaq_run_000167"
_CONVERTED = os.environ.get(
    "EVBLD_CONVERTED_DIR",
    "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi")
_EVENTS = os.path.join(_CONVERTED, _RUN, "events")
_INPUT = os.environ.get("VALIDATE_INPUT", os.path.join(_EVENTS, f"ecal_{_RUN}.root"))
_OUTPUT = os.environ.get("VALIDATE_OUTPUT", os.path.join(_EVENTS, f"ecal_{_RUN}.edm4hep.root"))
os.makedirs(os.path.dirname(_OUTPUT), exist_ok=True)

_f = ROOT.TFile.Open(_INPUT)
if not _f or _f.IsZombie():
    raise SystemExit(f"input not found (event-build run 167 first): {_INPUT}")
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
