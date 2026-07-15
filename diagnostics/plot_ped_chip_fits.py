"""Per-channel pedestal fits of one full chip (64 channels), with the Gaussian
fit overlaid -- the pedestal analogue of plot_mip_chip_shapes.py. For each
channel the SCA with the most entries is shown (pedestals are fit per (channel,
SCA); one representative SCA per channel keeps the gallery readable). The fit
mirrors PedestalMipCalib.h::fitPedestalSca: locate the peak, fit a Gaussian in a
narrow (+/-8) and a wide (+/-16) window, keep the narrower sigma.

Usage:  python3 diagnostics/plot_ped_chip_fits.py <hist.root> <high|low> <slab> <chip> [outdir]
"""
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kError

_HERE = os.path.dirname(os.path.abspath(__file__))
HIST = sys.argv[1]
GAIN = sys.argv[2] if len(sys.argv) > 2 else "low"
SLAB = int(sys.argv[3]) if len(sys.argv) > 3 else 13
CHIP = int(sys.argv[4]) if len(sys.argv) > 4 else 5
OUT = sys.argv[5] if len(sys.argv) > 5 else os.path.join(_HERE, "th210")
os.makedirs(OUT, exist_ok=True)
TH = os.environ.get("TH", "210")
NSCA = 15


def best_sca_hist(f, chn):
    """The (SCA) pedestal histogram of this channel with the most entries."""
    best, bestn = None, 0
    for sca in range(NSCA):
        h = f.Get(f"ped_{GAIN}_s{SLAB}_c{CHIP}_ch{chn}_sca{sca}")
        if h and h.GetEntries() > bestn:
            best, bestn = h, h.GetEntries()
    return best


def fit_pedestal(h, tag):
    """Gaussian pedestal fit, mirroring fitPedestalSca (narrow vs wide window)."""
    peak = h.GetXaxis().GetBinCenter(h.GetMaximumBin())
    f0 = ROOT.TF1(f"g0_{tag}", "gaus", peak - 8, peak + 8)
    f1 = ROOT.TF1(f"g1_{tag}", "gaus", peak - 16, peak + 16)
    h.Fit(f0, "RQ0")
    h.Fit(f1, "RQ0")
    chosen = f0 if f0.GetParameter(2) < f1.GetParameter(2) else f1
    return chosen, chosen.GetParameter(1), chosen.GetParameter(2)


f = ROOT.TFile.Open(HIST)
if not f or f.IsZombie():
    raise SystemExit(f"cannot open {HIST}")

c = ROOT.TCanvas("c", "ped fits", 2600, 2600)
c.Divide(8, 8, 0.001, 0.001)
keep = []

for chn in range(64):
    pad = c.cd(chn + 1)
    pad.SetMargin(0.10, 0.03, 0.10, 0.08)
    h = best_sca_hist(f, chn)
    if not h or h.GetEntries() < 50:
        continue
    h = h.Clone(f"ped_{chn}"); h.SetDirectory(0)
    peak = h.GetXaxis().GetBinCenter(h.GetMaximumBin())
    h.GetXaxis().SetRangeUser(peak - 25, peak + 25)
    h.SetTitle(f"ch{chn};ADC;")
    h.SetLineColor(ROOT.kGray + 2)
    h.SetFillColorAlpha(ROOT.kGray, 0.35)
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetLabelSize(0.06)
    h.GetXaxis().SetTitleSize(0.06)
    h.Draw("hist")
    g, mean, sig = fit_pedestal(h, f"{chn}")
    g.SetLineColor(ROOT.kAzure + 1); g.SetLineWidth(2); g.SetNpx(200)
    g.Draw("l same")
    txt = ROOT.TLatex(); txt.SetNDC(True); txt.SetTextSize(0.075); txt.SetTextColor(ROOT.kAzure + 2)
    txt.DrawLatex(0.13, 0.84, f"#mu {mean:.1f}")
    txt.DrawLatex(0.13, 0.74, f"#sigma {sig:.2f}")
    keep += [h, g, txt]

c.cd(0)
title = ROOT.TLatex(); title.SetNDC(True); title.SetTextSize(0.009); title.SetTextAlign(23)
title.DrawLatex(0.5, 0.998, f"th{TH} pedestal fits (Gaussian)   {GAIN} gain   slab {SLAB} chip {CHIP}   (best-SCA per channel)")
keep.append(title)

out = os.path.join(OUT, f"ped_fits_{GAIN}_s{SLAB}_c{CHIP}.png")
c.SaveAs(out)
print(f"saved {out}")
