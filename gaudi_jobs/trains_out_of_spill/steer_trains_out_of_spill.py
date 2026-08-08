#
# Stage 2 of the trains_out_of_spill pipeline: event-build the extracted
# acquisitions into an `ecal` tree the event_viewer opens directly.
#
# These are the isolated above-cut acquisitions that sit BETWEEN spill trains
# on run_000060 (th220) -- see diagnostics/plot_isolated_bursts.py and section 9
# of diagnostics/pedestal_multipeak_tables.txt for what they turned out to be.
# Stage 1 (extract_acquisitions.py) has already copied them into one small
# siwecaldecoded file; this is the ordinary event builder over that file, with
# the full th220 calibration, so nothing about the reconstruction is special to
# this pipeline.
#
# Calibration is th220's own (run_000060 is a th220 run) from
# calibration/MuonCalib_gaudi/: pedestals/th220/, mips/th220/ and the LG->HG
# anchor anchor/th220/. Per-threshold, nothing borrowed -- same rule as every
# other pipeline in this repo.
#
# Usage (from the repo root):
#   source setup.sh
#   export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
#   export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
#   k4run gaudi_jobs/trains_out_of_spill/steer_trains_out_of_spill.py
#
# Or run both stages via gaudi_jobs/trains_out_of_spill/run_pipeline.sh.
#
import glob
import os

from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalEventBuilder
from k4FWCore import ApplicationMgr

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))
_TH = "220"
_CALIB = os.path.join(_REPO, "calibration", "MuonCalib_gaudi")
_MAPPINGS = os.path.join(_REPO, "mappings")

_INPUT = os.environ.get("TOS_INPUT",
                        os.path.join(_HERE, "data", "decoded_trains_out_of_spill.root"))
_OUTPUT = os.environ.get("TOS_OUTPUT",
                         os.path.join(_HERE, "data", "ecal_trains_out_of_spill.root"))
if not os.path.exists(_INPUT):
    raise SystemExit(
        f"input not found: {_INPUT}\n"
        "Run stage 1 first (extract_acquisitions.py), or point TOS_INPUT at its output.")
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
    InputFiles=[_INPUT],
    TreeName="siwecaldecoded",
    OutputFile=_OUTPUT,
    PedestalFile=ped_hg,
    MipFile=mip_hg,
    PedestalFileLowGain=ped_lg,
    MipFileLowGain=mip_lg,
    AdcSaturationThreshold=switch,   # on (adc_high - ped_hg)
    GainRatio=k,
    GainIntercept=c,
    # Geometry: without the pad map and the per-slab z the hits carry no
    # position and the viewer comes up empty. Slab 12 is the COB and needs its
    # own fev11 map -- which matters here, since one of the three acquisitions
    # is an all-of-slab-12 event.
    PadMapDefaultFile=os.path.join(_MAPPINGS, "fev10_rotate_chip_channel_x_y_mapping.txt"),
    PadMapSlabOverrides=[f"12:{os.path.join(_MAPPINGS, 'fev11_cob_good_rotate_chip_channel_x_y_mapping.txt')}"],
    SlabZFile=os.path.join(_MAPPINGS, "slab_z_positions.yml"),
    # THE ONE SETTING THIS PIPELINE DELIBERATELY CHANGES.
    # The default MinSlabsHit=10 asks an event to light 10 of the 15 slabs --
    # right for physics, since it is what rejects noise and keeps particles
    # that cross the stack. Here it throws away exactly what we came to look
    # at: two of the three acquisitions are confined to 2 and 3 slabs (a
    # whole-slab-12 coincidence and a chip retrigger in slab 1). With the
    # default, building these three acquisitions yields 7 events and ALL of
    # them come from the one detector-wide acquisition -- the other two vanish
    # silently. Set to 1 so localised events survive.
    # Consequence to keep in mind while looking at the result: this is NOT the
    # physics selection, so the events here are not comparable with those in a
    # normal ecal_*.root, and plenty of them are noise by construction.
    MinSlabsHit=int(os.environ.get("TOS_MIN_SLABS_HIT", "1")),
)

ApplicationMgr(TopAlg=[builder],
               EvtSel="NONE",
               EvtMax=1,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
