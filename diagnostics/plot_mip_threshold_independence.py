"""Is the MIP really threshold-independent? Measure it, don't argue it.

The whole calibration now rests on one claim: the MIP MPV in ADC is a property of the
PREAMP, not of the discriminator, so th220's MIP table is valid for a th230 run. The
evidence so far is indirect -- the two spectra overlap on their tails, which they could
not do if the gains differed -- and that is an argument about SHAPE, not a number.

The difficulty is obvious: you cannot just fit th230's MIP spectrum and compare, because
th230's trigger has EATEN the peak. That bias is the thing being corrected, so a fit
that includes the peak region is circular.

So fit where the trigger cannot reach. th230's discriminator cuts at 18-22 ADC; th220's
at 11-13. Above ~35 ADC NEITHER trigger has removed anything -- both spectra are intact
and comparable. A Landau-Gauss MPV is a SHAPE parameter, so fitting the tail alone still
constrains it. Do that, identically, to both thresholds:

    if the tail-only MPVs agree  -> the true MIP does not depend on the threshold, and
                                    using th220's table for a th230 run is correct.
    if they differ               -> the gains really are different, the whole
                                    "th220 MIP for everyone" decision is wrong, and the
                                    energy scale of every th230 run is wrong with it.

For contrast the panel also shows the STANDARD full-range fits, which give 21.69 (th220)
and 32.40 (th230) -- the 1.47x discrepancy this is meant to explain.

Panels:
  1. Both spectra, tail-normalised, with the tail-only fit region shaded.
  2. The tail-only fits drawn on top of the data.
  3. The four MPVs side by side: full-range vs tail-only, each threshold.
  4. Ratio of the two spectra, bin by bin. Flat and at 1 above the fit region = same
     gain; the collapse below it is th230's missing peak.

Usage:  python3 diagnostics/plot_mip_threshold_independence.py
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUT, exist_ok=True)

H = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/calib_fill_scratch/hist"
SRC = {
    "th230": (f"{H}/th230/merged_th230.root", ROOT.kRed + 1),
    "th220": (f"{H}/TB2026CERN_run_000060/merged_run.root", ROOT.kBlue + 1),
}
SLABS = [int(s) for s in os.environ.get("SLABS", "5,8,11,14").split(",")]

# Where NEITHER trigger has removed anything. th230 cuts at 18-22 ADC, th220 at 11-13;
# 35 is clear of both with room to spare. The upper end is where the statistics die.
TAIL_LO = float(os.environ.get("TAIL_LO", "35"))
TAIL_HI = float(os.environ.get("TAIL_HI", "140"))
# The full range the calibrator itself uses -- the biased fit, shown for contrast.
FULL_LO, FULL_HI = 8.0, 140.0


def langaufun(x, par):
    """Landau convolved with a Gaussian. par = [width, MPV, area, sigma]."""
    invsq2pi = 0.3989422804014327
    mpshift = -0.22278298
    npts, sc = 100, 5.0
    mpc = par[1] - mpshift * par[0]
    xlow = x[0] - sc * par[3]
    xupp = x[0] + sc * par[3]
    step = (xupp - xlow) / npts
    s = 0.0
    for i in range(1, npts // 2 + 1):
        xx = xlow + (i - 0.5) * step
        s += ROOT.TMath.Landau(xx, mpc, par[0]) / par[0] * ROOT.TMath.Gaus(x[0], xx, par[3])
        xx = xupp - (i - 0.5) * step
        s += ROOT.TMath.Landau(xx, mpc, par[0]) / par[0] * ROOT.TMath.Gaus(x[0], xx, par[3])
    return par[2] * step * s * invsq2pi / par[3]


def combined(path, tag):
    """Whole-detector MIP spectrum, each channel-SCA pedestal-subtracted first."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"ERROR: cannot open {path}")
    h = ROOT.TH1F(f"h_{tag}", "", 300, -20.5, 279.5)
    h.SetDirectory(0)
    n = 0
    for slab in SLABS:
        for chip in range(16):
            for chn in range(64):
                for sca in range(15):
                    mh = f.Get(f"mip_high_s{slab}_c{chip}_ch{chn}_sca{sca}")
                    ph = f.Get(f"ped_high_s{slab}_c{chip}_ch{chn}_sca{sca}")
                    if not mh or not ph or ph.GetEntries() == 0:
                        continue
                    ph.GetXaxis().SetRangeUser(ph.GetMean() - 20, ph.GetMean() + 20)
                    pm = ph.GetMean()
                    if pm <= 0:
                        continue
                    for k in range(mh.GetNbinsX() + 1):
                        y = mh.GetBinContent(k)
                        if y > 0:
                            h.Fill(mh.GetXaxis().GetBinCenter(k) - pm, y)
                    n += 1
    f.Close()
    print(f"  {tag}: pooled {n:,} channel-SCAs, {h.GetIntegral() and h.Integral():,.0f} entries")
    return h


def fit_mpv(h, lo, hi, tag):
    """Landau-Gauss MPV over [lo, hi]. Returns (mpv, err)."""
    fn = ROOT.TF1(f"lg_{tag}", langaufun, lo, hi, 4)
    # Seed the MPV from the histogram's own peak inside the window, so the tail-only
    # fit is not handed the answer we are testing for.
    sub = h.Clone(f"sub_{tag}")
    sub.SetDirectory(0)
    sub.GetXaxis().SetRangeUser(lo, hi)
    seed = sub.GetXaxis().GetBinCenter(sub.GetMaximumBin())
    fn.SetParameters(3.0, max(seed, 12.0), h.Integral(h.FindBin(lo), h.FindBin(hi)) * 2, 3.0)
    fn.SetParNames("Width", "MPV", "Area", "GSigma")
    fn.SetParLimits(0, 0.3, 15)
    fn.SetParLimits(1, 5, 90)
    fn.SetParLimits(3, 0.3, 20)
    h.Fit(fn, "RQ0", "", lo, hi)
    return fn.GetParameter(1), fn.GetParError(1), fn


print("[build] combined MIP spectra, slabs", SLABS)
hists, tail_fit, full_fit, res = {}, {}, {}, {}
for th, (path, colour) in SRC.items():
    h = combined(path, th)
    h.SetLineColor(colour)
    h.SetLineWidth(2)
    hists[th] = h

for th in SRC:
    m_t, e_t, f_t = fit_mpv(hists[th], TAIL_LO, TAIL_HI, f"tail_{th}")
    m_f, e_f, f_f = fit_mpv(hists[th], FULL_LO, FULL_HI, f"full_{th}")
    tail_fit[th], full_fit[th] = f_t, f_f
    res[th] = dict(tail=(m_t, e_t), full=(m_f, e_f))
    print(f"  {th}: full-range MPV = {m_f:6.2f} +/- {e_f:.2f}   "
          f"TAIL-ONLY MPV = {m_t:6.2f} +/- {e_t:.2f}")

r_full = res["th230"]["full"][0] / res["th220"]["full"][0]
r_tail = res["th230"]["tail"][0] / res["th220"]["tail"][0]
e_tail = r_tail * np.hypot(res["th230"]["tail"][1] / res["th230"]["tail"][0],
                           res["th220"]["tail"][1] / res["th220"]["tail"][0])

c = ROOT.TCanvas("c", "mip threshold independence", 1700, 1200)
c.Divide(2, 2, 0.005, 0.005)
keep = []

# Normalise on the tail so the two are comparable at all.
norm = {}
for th, h in hists.items():
    n = h.Integral(h.FindBin(TAIL_LO), h.FindBin(TAIL_HI))
    hh = h.Clone(f"n_{th}")
    hh.SetDirectory(0)
    hh.Scale(1.0 / n)
    norm[th] = hh
    keep.append(hh)

# ---- 1. the two spectra, with the fit region shaded ----------------------
pad = c.cd(1)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
pad.SetLogy()
ymax = max(h.GetMaximum() for h in norm.values())
for i, (th, h) in enumerate(norm.items()):
    h.SetTitle("Combined MIP spectrum, tail-normalised;ADC - pedestal;entries (norm.)")
    h.GetXaxis().SetRangeUser(-5, 150)
    h.GetYaxis().SetRangeUser(ymax * 1e-4, ymax * 5)
    h.Draw("hist" if i == 0 else "hist same")
band = ROOT.TBox(TAIL_LO, ymax * 1e-4, TAIL_HI, ymax * 5)
band.SetFillColorAlpha(ROOT.kGreen + 2, 0.12)
band.SetLineColor(ROOT.kGreen + 2)
band.SetLineWidth(2)
band.Draw("l")
leg = ROOT.TLegend(0.45, 0.66, 0.96, 0.89)
leg.SetBorderSize(0)
leg.SetFillColorAlpha(ROOT.kWhite, 0.80)
leg.SetTextSize(0.026)
for th, h in norm.items():
    leg.AddEntry(h, f"{th}", "l")
leg.AddEntry(band, f"fit region: ADC {TAIL_LO:.0f}-{TAIL_HI:.0f}", "f")
leg.AddEntry(0, "no trigger reaches here -- both spectra intact", "")
leg.Draw()
keep += [band, leg]

# ---- 2. the tail-only fits ----------------------------------------------
pad = c.cd(2)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
pad.SetLogy()
leg2 = ROOT.TLegend(0.40, 0.62, 0.96, 0.89)
leg2.SetBorderSize(0)
leg2.SetFillColorAlpha(ROOT.kWhite, 0.80)
leg2.SetTextSize(0.026)
for i, (th, h) in enumerate(norm.items()):
    h.SetTitle("TAIL-ONLY fits: the MPV where no trigger interfered;"
               "ADC - pedestal;entries (norm.)")
    h.GetXaxis().SetRangeUser(-5, 150)
    h.Draw("hist" if i == 0 else "hist same")
    scale = 1.0 / hists[th].Integral(hists[th].FindBin(TAIL_LO), hists[th].FindBin(TAIL_HI))
    f = tail_fit[th]
    g = ROOT.TF1(f"draw_{th}", lambda x, p, f=f, s=scale: f.Eval(x[0]) * s,
                 TAIL_LO - 12, TAIL_HI, 0)
    g.SetLineColor(SRC[th][1])
    g.SetLineWidth(4)
    g.SetLineStyle(2)
    g.Draw("same")
    m, e = res[th]["tail"]
    leg2.AddEntry(g, f"{th}: tail-only MPV = {m:.2f} #pm {e:.2f} ADC", "l")
    keep.append(g)
leg2.AddEntry(0, "", "")
leg2.AddEntry(0, f"ratio th230/th220 = {r_tail:.3f} #pm {e_tail:.3f}", "")
verdict = "SAME MIP -- threshold-independent" if abs(r_tail - 1) < 0.05 else "MIPs DIFFER"
leg2.AddEntry(0, f"-> {verdict}", "")
leg2.Draw()
keep.append(leg2)

# ---- 3. the four MPVs ----------------------------------------------------
pad = c.cd(3)
pad.SetMargin(0.13, 0.05, 0.16, 0.10)
fr = ROOT.TH1F("fr", "MPV: what the fit says, and what the trigger did to it;;MPV  [ADC]",
               4, 0, 4)
labels = [("th220 full", res["th220"]["full"], ROOT.kBlue + 1),
          ("th230 full", res["th230"]["full"], ROOT.kRed + 1),
          ("th220 TAIL", res["th220"]["tail"], ROOT.kBlue + 1),
          ("th230 TAIL", res["th230"]["tail"], ROOT.kRed + 1)]
for i, (lab, _, _) in enumerate(labels):
    fr.GetXaxis().SetBinLabel(i + 1, lab)
fr.GetXaxis().SetLabelSize(0.05)
fr.SetMinimum(0)
fr.SetMaximum(max(v[0] for _, v, _ in labels) * 1.45)
fr.Draw()
keep.append(fr)
for i, (lab, (m, e), colour) in enumerate(labels):
    g = ROOT.TGraphErrors(1, np.array([i + 0.5]), np.array([m]),
                          np.array([0.0]), np.array([e]))
    g.SetMarkerStyle(21 if i < 2 else 20)
    g.SetMarkerSize(2.0)
    g.SetMarkerColor(colour)
    g.SetLineColor(colour)
    g.SetLineWidth(3)
    g.Draw("P same")
    t = ROOT.TLatex(i + 0.5, m + fr.GetMaximum() * 0.05, f"{m:.2f}")
    t.SetTextAlign(21)
    t.SetTextSize(0.042)
    t.SetTextColor(colour)
    t.Draw()
    keep += [g, t]
t3 = ROOT.TLatex()
t3.SetTextSize(0.038)
t3.SetNDC()
t3.DrawLatex(0.18, 0.84, f"full-range ratio  th230/th220 = {r_full:.3f}   <- the bias")
t3.DrawLatex(0.18, 0.78, f"TAIL-ONLY  ratio  th230/th220 = {r_tail:.3f} #pm {e_tail:.3f}")
keep.append(t3)

# ---- 4. bin-by-bin ratio -------------------------------------------------
pad = c.cd(4)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
rat = norm["th230"].Clone("rat")
rat.SetDirectory(0)
rat.Divide(norm["th220"])
rat.SetTitle("Ratio th230 / th220, bin by bin;ADC - pedestal;th230 / th220")
rat.GetXaxis().SetRangeUser(0, 120)
rat.GetYaxis().SetRangeUser(0, 2.0)
rat.SetLineColor(ROOT.kBlack)
rat.SetLineWidth(2)
rat.SetMarkerStyle(20)
rat.SetMarkerSize(0.7)
rat.Draw("P")
one = ROOT.TLine(0, 1, 120, 1)
one.SetLineColor(ROOT.kGreen + 2)
one.SetLineWidth(3)
one.SetLineStyle(2)
one.Draw()
b4 = ROOT.TBox(TAIL_LO, 0, min(TAIL_HI, 120), 2.0)
b4.SetFillColorAlpha(ROOT.kGreen + 2, 0.10)
b4.Draw()
rat.Draw("P same")
leg4 = ROOT.TLegend(0.34, 0.68, 0.96, 0.89)
leg4.SetBorderSize(0)
leg4.SetFillColorAlpha(ROOT.kWhite, 0.80)
leg4.SetTextSize(0.026)
leg4.AddEntry(one, "1 = identical spectra", "l")
leg4.AddEntry(b4, "region no trigger touches", "f")
leg4.AddEntry(0, "flat at 1 here = SAME GAIN", "")
leg4.AddEntry(0, "the collapse at low ADC is th230's missing peak", "")
leg4.Draw()
keep += [rat, one, b4, leg4]

out = os.path.join(OUT, "mip_threshold_independence.png")
c.SaveAs(out)

print("\n" + "=" * 78)
print("  IS THE MIP THRESHOLD-INDEPENDENT?")
print("=" * 78)
print(f"\n  {'':10s} {'full-range MPV':>18s} {'TAIL-ONLY MPV':>18s}")
print("  " + "-" * 50)
for th in ("th220", "th230"):
    mf, ef = res[th]["full"]
    mt, et = res[th]["tail"]
    print(f"  {th:10s} {mf:12.2f} +/-{ef:4.2f} {mt:12.2f} +/-{et:4.2f}")
print(f"\n  full-range ratio th230/th220 = {r_full:.3f}   <- the 1.47x bias we are explaining")
print(f"  TAIL-ONLY  ratio th230/th220 = {r_tail:.3f} +/- {e_tail:.3f}")
print(f"\n  The tail is the region NEITHER trigger reaches (ADC {TAIL_LO:.0f}-{TAIL_HI:.0f}).")
if abs(r_tail - 1) < 0.05:
    print(f"  -> The MPVs agree to {100 * abs(r_tail - 1):.1f}%. The true MIP does NOT depend on")
    print(f"     the threshold. Using th220's MIP table for a th230 run is CORRECT.")
else:
    print(f"  -> The MPVs DIFFER by {100 * abs(r_tail - 1):.1f}%. The gains are not the same, and")
    print(f"     'th220 MIP for everyone' is WRONG. Stop and rethink.")
print(f"\nsaved {out}")
