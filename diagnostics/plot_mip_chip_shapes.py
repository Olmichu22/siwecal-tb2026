"""Per-channel combined MIP spectra of one full chip (64 channels), with the
full-range langau fit overlaid -- for eyeballing the SHAPE of the distributions
in either gain. The histogram reconstruction (per-SCA truncated-mean pedestal
subtraction) and the fit setup are replicated verbatim from
PedestalMipCalibrator.cpp / PedestalMipCalib.h::fitLangau, so this is what the
calibrator actually fits.

Usage:  python3 diagnostics/plot_mip_chip_shapes.py <hist.root> <high|low> <slab> <chip> [outdir]
"""
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gErrorIgnoreLevel = ROOT.kError

_HERE = os.path.dirname(os.path.abspath(__file__))
HIST = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    _HERE, "..", "calibration", "MuonCalib_gaudi_test_nsh5", "mips", "th230", "hist_run004_nsh5.root")
GAIN = sys.argv[2] if len(sys.argv) > 2 else "high"
SLAB = int(sys.argv[3]) if len(sys.argv) > 3 else 13
CHIP = int(sys.argv[4]) if len(sys.argv) > 4 else 5
OUT = sys.argv[5] if len(sys.argv) > 5 else os.path.join(_HERE, "th230")
os.makedirs(OUT, exist_ok=True)
HIGH = GAIN == "high"
NSCA = 15
TH = os.environ.get("TH", "230")  # label only

# langau parameter limits, verbatim from PedestalMipCalib.h::fitLangau.
if HIGH:
    PLLO = [0.1, 14.0, 1.0, 0.5]
    PLHI = [20.0, 220.0, 1.0e8, 10.0]
    XVIEW = (0.0, 80.0)
else:
    PLLO = [0.1, 2.0, 1.0, 0.05]
    PLHI = [5.0, 15.0, 1.0e8, 2.0]
    XVIEW = (-5.0, 25.0)

ROOT.gInterpreter.Declare(r'''
double langaufun_cpp(double* x, double* par) {
  const double invsq2pi = 0.3989422804014, mpshift = -0.22278298, np = 100.0, sc = 5.0;
  const double mpc = par[1] - mpshift * par[0];
  const double xlow = x[0] - sc * par[3], xupp = x[0] + sc * par[3];
  const double step = (xupp - xlow) / np;
  double sum = 0.0;
  for (double i = 1.0; i <= np / 2; i += 1.0) {
    double xx = xlow + (i - .5) * step;
    double fl = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fl * TMath::Gaus(x[0], xx, par[3]);
    xx = xupp - (i - .5) * step;
    fl = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fl * TMath::Gaus(x[0], xx, par[3]);
  }
  return par[2] * step * sum * invsq2pi / par[3];
}
''')


def clamp(v, lo, hi):
    return min(max(v, lo), hi)


def combined_hist(f, chn):
    h = ROOT.TH1F(f"h_{chn}", "", 900, -100.5, 799.5)
    h.SetDirectory(0)
    for sca in range(NSCA):
        mh = f.Get(f"mip_{GAIN}_s{SLAB}_c{CHIP}_ch{chn}_sca{sca}")
        ph = f.Get(f"ped_{GAIN}_s{SLAB}_c{CHIP}_ch{chn}_sca{sca}")
        if not mh or not ph or ph.GetEntries() <= 0:
            continue
        ph.GetXaxis().SetRangeUser(ph.GetMean() - 20, ph.GetMean() + 20)
        pm = ph.GetMean()
        if pm <= 0:
            continue
        for k in range(900):
            y = mh.GetBinContent(k)
            if y > 0:
                h.Fill(int(mh.GetXaxis().GetBinCenter(k) - pm), y)
    return h


def fit_fullrange(hist, tag):
    """Full-range langau fit, mirroring fitLangau for the selected gain."""
    if HIGH:
        fr0, fr1 = 0.0, 100.0
        hist.GetXaxis().SetRangeUser(fr0, fr1)
    else:
        fr0 = max(hist.GetMean() - 2.0 * hist.GetRMS(), 2.0)
        hist.GetXaxis().SetRangeUser(fr0, 25.0)
        fr1 = hist.GetMean() * 1.1 + 2.0 * hist.GetRMS()
    lan = ROOT.TF1(f"lan_{tag}", "landau", fr0, fr1)
    hist.Fit(lan, "QRE")
    hist.GetXaxis().SetRangeUser(fr0, fr1)
    sv = [clamp(lan.GetParameter(2), PLLO[0], PLHI[0]),
          clamp(lan.GetParameter(1), PLLO[1], PLHI[1]),
          hist.Integral("width") * 1.2,
          clamp(hist.GetRMS() * 0.1, PLLO[3], PLHI[3])]
    lg = ROOT.TF1(f"langau_{tag}", ROOT.langaufun_cpp, fr0, fr1, 4)
    lg.SetParameters(*sv)
    for i in range(4):
        lg.SetParLimits(i, PLLO[i], PLHI[i])
    hist.Fit(lg, "RBQM")
    return lg, lg.GetParameter(1)


f = ROOT.TFile.Open(HIST)
if not f or f.IsZombie():
    raise SystemExit(f"cannot open {HIST}")

c = ROOT.TCanvas("c", "mip shapes", 2600, 2600)
c.Divide(8, 8, 0.001, 0.001)
keep = []

for chn in range(64):
    pad = c.cd(chn + 1)
    pad.SetMargin(0.10, 0.03, 0.10, 0.08)
    h = combined_hist(f, chn)
    h.GetXaxis().SetRangeUser(*XVIEW)
    h.SetTitle(f"ch{chn};ADC (ped-sub);")
    h.SetLineColor(ROOT.kGray + 2)
    h.SetFillColorAlpha(ROOT.kGray, 0.35)
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetLabelSize(0.06)
    h.GetXaxis().SetTitleSize(0.06)
    h.Draw("hist")
    keep.append(h)
    if h.Integral() < 200:
        continue
    hf = h.Clone(f"hf_{chn}"); hf.SetDirectory(0)
    g, mpv = fit_fullrange(hf, f"{chn}")
    g.SetLineColor(ROOT.kRed + 1); g.SetLineWidth(2); g.SetNpx(300)
    g.Draw("l same")
    ln = ROOT.TLine(mpv, 0, mpv, h.GetMaximum()); ln.SetLineColor(ROOT.kRed + 1); ln.SetLineStyle(2); ln.Draw()
    txt = ROOT.TLatex(); txt.SetNDC(True); txt.SetTextSize(0.085); txt.SetTextColor(ROOT.kRed + 1)
    txt.DrawLatex(0.45, 0.84, f"MPV {mpv:.1f}")
    keep += [hf, g, ln, txt]

c.cd(0)
title = ROOT.TLatex(); title.SetNDC(True); title.SetTextSize(0.009); title.SetTextAlign(23)
title.DrawLatex(0.5, 0.998, f"th{TH} MIP spectra + full-range langau   {GAIN} gain   slab {SLAB} chip {CHIP}")
keep.append(title)

out = os.path.join(OUT, f"mip_shapes_{GAIN}_s{SLAB}_c{CHIP}.png")
c.SaveAs(out)
print(f"saved {out}")
