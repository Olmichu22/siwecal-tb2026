"""Energy calibration reproducibility: the FIT run against the TEST run, same
threshold, same per-threshold calibration. If the calibration is sound the two
runs' sum_energy distributions land on top of each other (same beam energy) or
scale linearly (different energy); a shift means the calibration does not
transport between runs.

Reads the reconstructed ``ecal`` tree's hit_energy per event. Overlays the two
runs normalised to unit area, and reports mean / sigma / resolution for each.

Usage:  FIT=<ecal_fit.root> TEST=<ecal_test.root> TH=210 python3 diagnostics/plot_energy_fit_vs_test.py [outdir]
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
TH = os.environ.get("TH", "210")
FIT = os.environ["FIT"]
TEST = os.environ["TEST"]
FIT_LABEL = os.environ.get("FIT_LABEL", "fit")
TEST_LABEL = os.environ.get("TEST_LABEL", "test")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)


def sum_energy(path):
    f = ROOT.TFile.Open(path)
    t = f.Get("ecal") if f else None
    if not t:
        raise SystemExit(f"no 'ecal' tree in {path}")
    vals = []
    for e in t:
        s = 0.0
        for v in e.hit_energy:
            s += v
        vals.append(s)
    return np.array(vals)


ef = sum_energy(FIT)
et = sum_energy(TEST)

c = ROOT.TCanvas("c", "fit vs test", 1000, 780)
c.SetMargin(0.12, 0.05, 0.12, 0.09)
hi = float(np.percentile(np.concatenate([ef, et]), 99.5))
keep = []
BLUE, RED = ROOT.kAzure + 1, ROOT.kRed + 1
stats = {}
for arr, lab, col in ((ef, FIT_LABEL, BLUE), (et, TEST_LABEL, RED)):
    h = ROOT.TH1F(f"h_{lab}", f"th{TH}: energy, fit vs test;sum_energy  [MIP];fraction of events",
                  80, 0, hi)
    h.SetDirectory(0)
    for v in arr:
        h.Fill(v)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(col)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(col, 0.18)
    stats[lab] = (float(np.mean(arr)), float(np.std(arr)), len(arr))
    keep.append((h, lab, col))

ymax = max(h.GetMaximum() for h, _, _ in keep)
leg = ROOT.TLegend(0.52, 0.72, 0.94, 0.89)
leg.SetBorderSize(0)
leg.SetTextSize(0.028)
for i, (h, lab, col) in enumerate(keep):
    h.GetYaxis().SetRangeUser(0, ymax * 1.25)
    h.Draw("hist" if i == 0 else "hist same")
    mu, sig, n = stats[lab]
    leg.AddEntry(h, f"{lab}: #mu={mu:.0f}, #sigma/#mu={sig / mu:.3f}  (N={n:,})", "f")
leg.Draw()
keep.append(leg)

muf = stats[FIT_LABEL][0]
mut = stats[TEST_LABEL][0]
d = ROOT.TLatex()
d.SetNDC(True)
d.SetTextSize(0.028)
d.DrawLatex(0.52, 0.66, f"mean test/fit = {mut / muf:.3f}  ({100 * (mut / muf - 1):+.1f}%)")

out = os.path.join(OUT, f"energy_fit_vs_test_th{TH}.png")
c.SaveAs(out)
print(f"th{TH}: {FIT_LABEL} mu={muf:.1f} sig/mu={stats[FIT_LABEL][1] / muf:.3f} | "
      f"{TEST_LABEL} mu={mut:.1f} sig/mu={stats[TEST_LABEL][1] / mut:.3f} | test/fit={mut / muf:.3f}")
print(f"saved {out}")
