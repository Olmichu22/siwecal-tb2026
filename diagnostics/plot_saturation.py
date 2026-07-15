"""High-gain saturation: adc_high vs adc_low, hit by hit.

The high gain stays linear up to adc_high ~ 1950 and bends over there; the two
horizontal lines show the old (1200) and new (1900) switch-to-low-gain points.

The band is traced with the per-column MODE, not the mean: a second, diffuse
population sits at low adc_high across the whole adc_low range, and it drags a
ProfileX mean downwards at large adc_low, faking a turnover that is not there.

  python3 plot_saturation.py [outdir]   (default: this directory)
"""
import array
import glob
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gErrorIgnoreLevel = ROOT.kFatal
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.environ.get("RUN", "TB2026CERN_run_000013")
TH = os.environ.get("TH", "230")   # which threshold folder this run belongs to
# Chain the run's decoded chunks. There is no merged <RUN>.root any more -- the
# chunks ARE the decoded data (the event builder chains them too).
CHUNKS = f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi/{RUN}/chunks"
OLD_SWITCH, NEW_SWITCH = 1200, 1900
N_ACQ = 8000  # acquisitions read; plenty for the band, keeps this a ~1 min job

t = ROOT.TChain("siwecaldecoded")
for chunk in sorted(glob.glob(os.path.join(CHUNKS, "chunk_*.root")))[:40]:
    t.AddFile(chunk)
if t.GetEntries() == 0:
    raise SystemExit(f"ERROR: no decoded chunks under {CHUNKS}")

c = ROOT.TCanvas("c", "saturation", 1250, 900)
ROOT.gPad.SetMargin(0.11, 0.14, 0.11, 0.07)
ROOT.gPad.SetLogz()

h2 = ROOT.TH2F("h2", f"{RUN}: high gain vs low gain, hit by hit;adc_low (raw ADC);adc_high (raw ADC)",
               210, 150, 680, 260, 0, 2600)
t.Draw("adc_high:adc_low>>h2", "hitbit_high==1", "goff", N_ACQ)
h2.Draw("colz")

# Ridge = the MEDIAN of adc_high in each adc_low column, over ALL the hits.
#
# No cut. An earlier version took the per-column MODE after cutting adc_high > 700
# to "isolate the band", and that cut fabricated its own answer: below adc_low
# ~300 the band has not risen above 700 yet, so the cut kept nothing but each
# column's upper tail and the mode pinned itself to the cut edge -- drawing a flat
# row of points at 700 straight through a region where the data is visibly piled
# at 250. The median of the untouched column cannot do that: where the hits are at
# 250, it says 250.
ridge = ROOT.TGraph()
q = array.array("d", [0.0])
prob = array.array("d", [0.5])
for ix in range(1, h2.GetNbinsX() + 1):
    col = h2.ProjectionY("col", ix, ix)
    if col.GetEntries() < 50:
        continue
    col.GetQuantiles(1, q, prob)
    ridge.SetPoint(ridge.GetN(), h2.GetXaxis().GetBinCenter(ix), q[0])
ridge.SetMarkerStyle(20)
ridge.SetMarkerSize(0.6)
ridge.SetMarkerColor(ROOT.kBlack)
ridge.Draw("P same")

keep = []
leg = ROOT.TLegend(0.15, 0.79, 0.42, 0.90)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.028)
for y, colour, style, label in ((OLD_SWITCH, ROOT.kRed + 1, 2, "previous switch"),
                                (NEW_SWITCH, ROOT.kGreen + 2, 1, "new switch")):
    ln = ROOT.TLine(150, y, 680, y)
    ln.SetLineColor(colour)
    ln.SetLineWidth(3)
    ln.SetLineStyle(style)
    ln.Draw()
    keep.append(ln)
    leg.AddEntry(ln, label, "l")
leg.Draw()
keep.append(leg)

# Per-threshold folder AND a run-stamped name. It used to be a flat, fixed
# "saturation_hg_vs_lg.png", so running this on a second run silently overwrote the
# first one's plot -- which is exactly what happened with run72 and run13.
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)
png = os.path.join(OUT, f"saturation_hg_vs_lg_{RUN.split('_')[-1]}.png")
c.SaveAs(png)
print("saved", png)
