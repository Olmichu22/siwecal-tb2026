#
# Gaudi/k4FWCore steering para TB2026CERN_run_000013 — modo físico.
# Hit MIP cut = 0.5; sin bloques de variantes mip05_/mip1_.
#
# Uso (desde la raíz del repo):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/run000013/steer_run000013.py
#
import ROOT
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalToEDM4hep, EcalPidTransformer
from k4FWCore import ApplicationMgr, IOSvc

_RUN   = "TB2026CERN_run_000013"
# Event-builder output with real calibration (MuonCalib_it2, th220)
_INPUT = f"/eos/user/o/oarquero/TB2026CERN/data/{_RUN}/ecal_{_RUN}_realvalues.root"
_OUT   = f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Reconstruction/ecal_{_RUN}.edm4hep.root"

_f = ROOT.TFile.Open(_INPUT)
n_events = int(_f.Get("ecal").GetEntries())
_f.Close()

svc = IOSvc("IOSvc")
svc.Output = _OUT

source = EcalToEDM4hep("EcalToEDM4hep",
                        InputFile=_INPUT,
                        TreeName="ecal",
                        HitMipCut=0.5)      # <0 para deshabilitar el corte de hit

pid = EcalPidTransformer("EcalPidTransformer",
                          InputCaloHits=["ECalHits"],
                          OutputClusters=["ECalPid"],
                          MipThresholds=[])  # [0.5, 1.0] para modo validación (slider)

ApplicationMgr(TopAlg=[source, pid],
               EvtSel="NONE",
               EvtMax=n_events,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
