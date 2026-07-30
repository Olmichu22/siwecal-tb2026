#
# Gaudi/k4FWCore steering: ecal tree -> EDM4hep, keeping only the events that
# pass a chain of cuts, with the cuts written into the output file name.
#
#   ECAL_FILE=/path/ecal_TB2026CERN_run_000012.root \
#   ECAL_CUTS="moliere:45:,weighte:5000:" \
#       k4run gaudi_source/options/run_event_cut.py
#
#   -> ecal_TB2026CERN_run_000012_moliere_gt45_weighte_gt5000.edm4hep.root
#
# ECAL_CUTS is a comma-separated list of "variable:min:max". Either bound may be
# left empty for a one-sided cut, and the name reflects which:
#
#   moliere:45:          -> moliere_gt45          (>= 45)
#   mip_likeness:0.1:0.8 -> mip_likeness_0.1_0.8  (0.1 <= x <= 0.8)
#   weighte::4000        -> weighte_lt4000        (<= 4000)
#
# Cuts are ANDed: one EcalEventCut instance per entry, each rejecting events on
# its own, so adding a cut is adding one more item to the list and one more
# string to the name. Nothing about the surviving events is rewritten -- they
# keep every collection they had, so the output is a drop-in for anything that
# reads the unfiltered file (siwecal_validation included).
#
# Variables are the base scalars of EcalShowerVars: nhit, zbary, energy,
# mip_likeness, weighte, bar_x, bar_y, bar_r, moliere, transverse_rms,
# is_shower, shower_start, shower_max, shower_end, shower_start_10,
# shower_end_10, shower_length, first_layer, last_layer, n_layers_hit,
# e_over_nhit.
#
import os

import ROOT
from Gaudi.Configuration import INFO
from Configurables import EventDataSvc
from Configurables import EcalToEDM4hep, EcalPidTransformer, EcalEventCut
from k4FWCore import ApplicationMgr, IOSvc

_NO_BOUND = 1e30  # matches EcalEventCut's "no bound" default

ecal_file = os.environ.get("ECAL_FILE", "")
tree_name = os.environ.get("ECAL_TREE", "ecal")
if not ecal_file:
    raise SystemExit("Set ECAL_FILE to the input ecal_<run>.root")

_stem = os.path.splitext(os.path.basename(ecal_file))[0]
_label = _stem[len("ecal_"):] if _stem.startswith("ecal_") else _stem


def _fmt(value):
    """45.0 -> '45', 0.1 -> '0.1'. Keeps file names readable and stable."""
    return f"{value:g}"


def parse_cuts(spec):
    """'var:min:max,...' -> [(name, min, max, name_fragment), ...]."""
    cuts = []
    for item in [s.strip() for s in spec.split(",") if s.strip()]:
        parts = item.split(":")
        if len(parts) != 3:
            raise SystemExit(f"Bad cut '{item}': expected variable:min:max "
                             f"(a bound may be empty, e.g. moliere:45:)")
        var, lo_s, hi_s = (p.strip() for p in parts)
        if not var:
            raise SystemExit(f"Bad cut '{item}': no variable name")
        if not lo_s and not hi_s:
            raise SystemExit(f"Bad cut '{item}': at least one bound is required")
        lo = float(lo_s) if lo_s else -_NO_BOUND
        hi = float(hi_s) if hi_s else _NO_BOUND
        if lo > hi:
            raise SystemExit(f"Bad cut '{item}': min > max, it can never pass")
        if lo_s and hi_s:
            frag = f"{var}_{_fmt(lo)}_{_fmt(hi)}"
        elif lo_s:
            frag = f"{var}_gt{_fmt(lo)}"
        else:
            frag = f"{var}_lt{_fmt(hi)}"
        cuts.append((var, lo, hi, frag))
    return cuts


cuts = parse_cuts(os.environ.get("ECAL_CUTS", ""))
if not cuts:
    raise SystemExit("Set ECAL_CUTS, e.g. ECAL_CUTS='moliere:45:,weighte:5000:'")

suffix = "_".join(frag for *_, frag in cuts)
out_file = os.environ.get(
    "ECAL_CUT_OUT",
    os.path.join(os.path.dirname(ecal_file), f"ecal_{_label}_{suffix}.edm4hep.root"),
)
print("Cuts:      ", ", ".join(f"{v} in [{lo:g}, {hi:g}]" for v, lo, hi, _ in cuts))
print("Output File:", out_file)

_f = ROOT.TFile.Open(ecal_file)
n_events = int(_f.Get(tree_name).GetEntries())
_f.Close()

svc = IOSvc("IOSvc")
svc.Output = out_file

# Physics mode by DEFAULT here (hit cut >= 0.5 MIP, no mip05_/mip1_ variant
# blocks), unlike run_pid.py, whose defaults are the visualizer ones. That is
# what run_pid_batch.py produces for the reconstruction area, so a cut file made
# here sits on the same energy scale as the unfiltered ecal_<run>.edm4hep.root
# it is meant to be compared with.
hit_mip_cut = float(os.environ.get("ECAL_HIT_MIP_CUT", "0.5"))
_raw_mip = os.environ.get("ECAL_MIP_THRESHOLDS", "")
mip_thresholds = [float(t) for t in _raw_mip.split(",") if t.strip()]

source = EcalToEDM4hep("EcalToEDM4hep", InputFile=ecal_file, TreeName=tree_name,
                       HitMipCut=hit_mip_cut)
pid = EcalPidTransformer("EcalPidTransformer",
                         InputCaloHits=["ECalHits"], OutputClusters=["ECalPid"],
                         MipThresholds=mip_thresholds)

cut_algs = []
for idx, (var, lo, hi, _) in enumerate(cuts):
    cut_algs.append(EcalEventCut(f"EcalEventCut_{idx}_{var}",
                                 InputClusters="ECalPid", Variable=var, Min=lo, Max=hi))

ApplicationMgr(TopAlg=[source, pid] + cut_algs,
               EvtSel="NONE",
               EvtMax=n_events,
               ExtSvc=[EventDataSvc("EventDataSvc")],
               OutputLevel=INFO)
