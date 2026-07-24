#
# Gaudi/k4FWCore steering for TB2026CERN_eudaq_run_000213 (th210, angular scan):
# siwecaldecoded chunks -> event-built `ecal` tree, with the FULL th210
# calibration taken from this project's calibration/MuonCalib_gaudi/ folder:
#   pedestals/th210/  MIP mips/th210/  and the LG->HG anchor anchor/th210/.
#
# Every table is th210's own (per-threshold calibration; nothing is borrowed
# from another threshold). The anchor line (k, c) and the pedestal-subtracted
# switch are read from anchor/th210/gain_anchor_th210.txt.
#
# Part of the angular_scan pipeline (gaudi_jobs/angular_scan/README.md):
# runs 179-239, event builder + validation, from already-converted chunks
# (no CONVERT step here -- see gaudi_jobs/condor/generate_convert_eudaq_dag.py
# for how these chunks were produced).
#
# Usage (from the repo root):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/angular_scan/steer_TB2026CERN_eudaq_run_000213.py
#
import glob
import os

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalEventBuilder
from k4FWCore import ApplicationMgr

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
_RUN = "TB2026CERN_eudaq_run_000213"
_TH = "210"
_CALIB = os.path.join(_REPO, "calibration", "MuonCalib_gaudi")
_MAPPINGS = os.path.join(_REPO, "mappings")

_CONVERTED = os.environ.get(
    "EVBLD_CONVERTED_DIR",
    "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi")
_INPUT_GLOB = os.path.join(_CONVERTED, _RUN, "chunks", "chunk_*.root")
# Event-builder output goes to <run>/events/, alongside the run's <run>/chunks/,
# NOT to the shared Reconstruction/ area. Override with EVBLD_OUTPUT.
_EVENTS_DIR = os.path.join(_CONVERTED, _RUN, "events")
_OUTPUT = os.environ.get("EVBLD_OUTPUT", os.path.join(_EVENTS_DIR, f"ecal_{_RUN}.root"))
os.makedirs(os.path.dirname(_OUTPUT), exist_ok=True)


def _latest(pattern):
    hits = sorted(glob.glob(pattern))
    if not hits:
        raise SystemExit(f"calibration file not found: {pattern}")
    return hits[-1]


def _anchor(path):
    """Read k, c, switch from an anchor/thN/gain_anchor_thN.txt file."""
    k, c, switch = 0.0962, 1.45, 1500
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        tok = line.split()
        if tok[0] == "k":
            k = float(tok[1])
        elif tok[0] == "c":
            c = float(tok[1])
        elif tok[0] == "switch":
            switch = int(float(tok[1]))
    return k, c, switch


ped_hg = _latest(os.path.join(_CALIB, "pedestals", f"th{_TH}", "Pedestal_*_highgain.txt"))
mip_hg = _latest(os.path.join(_CALIB, "mips", f"th{_TH}", "MIP_pedestalsubmode1_*_highgain.txt"))
ped_lg = _latest(os.path.join(_CALIB, "pedestals", f"th{_TH}", "Pedestal_*_lowgain.txt"))
mip_lg = _latest(os.path.join(_CALIB, "mips", f"th{_TH}", "MIP_pedestalsubmode1_*_lowgain.txt"))
k, c, switch = _anchor(os.path.join(_CALIB, "anchor", f"th{_TH}", f"gain_anchor_th{_TH}.txt"))

builder = EcalEventBuilder(
    "EcalEventBuilder",
    InputFiles=sorted(glob.glob(_INPUT_GLOB)),
    TreeName="siwecaldecoded",
    OutputFile=_OUTPUT,
    PedestalFile=ped_hg,
    MipFile=mip_hg,
    PedestalFileLowGain=ped_lg,
    MipFileLowGain=mip_lg,
    AdcSaturationThreshold=switch,   # on (adc_high - ped_hg)
    GainRatio=k,
    GainIntercept=c,
    # Geometry: pad (x,y) mapping + per-slab z. Without these the hits have no
    # position and the event display comes up empty. Slab 12 is the COB and uses
    # its own pad map (fev11), every other slab the default fev10.
    PadMapDefaultFile=os.path.join(_MAPPINGS, "fev10_rotate_chip_channel_x_y_mapping.txt"),
    PadMapSlabOverrides=[f"12:{os.path.join(_MAPPINGS, 'fev11_cob_good_rotate_chip_channel_x_y_mapping.txt')}"],
    SlabZFile=os.path.join(_MAPPINGS, "slab_z_positions.yml"),
)

ApplicationMgr(TopAlg=[builder],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
