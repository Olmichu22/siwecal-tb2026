"""Is a MIP measurable in the low gain at all? Signal against noise, per channel.

The low-gain MIP calibration has been chased through two explanations already --
the trigger truncates the peak, the Langau fit is biased. Both are true. Neither
is the deepest one, which is simply this:

    a MIP in the low gain is 2.2 ADC, and the low-gain noise is 0.9 ADC.

That is a signal-to-noise of 2.4. In the high gain the same MIP sits at 21.7 ADC
against 1.95 ADC of noise -- 11.1 sigma, a peak you can see with your eyes.

And the low-gain noise does NOT scale down with the gain. If the low gain were
merely the high gain divided by ten, its noise would be 1.95/10 = 0.19 ADC. It is
0.92 -- nearly five times more. The preamp's low-gain branch has its own noise
floor, independent of the signal it is amplifying. So dividing the signal by ten
does not divide the noise by ten: it destroys the signal-to-noise.

Three panels:
  1. Where the MIP sits relative to the noise, HIGH gain: the peak stands clear.
  2. The same for the LOW gain: the MIP is inside the noise distribution.
  3. Signal-to-noise per channel, both gains, on the same axis.

Usage:  TH=th230 python3 diagnostics/plot_mip_snr.py [outdir]
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(_HERE)
TH = os.environ.get("TH", "th230")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, TH)
os.makedirs(OUT, exist_ok=True)

K_ELEC = 0.0999   # measured gain ratio, raw data, 21,101 hits in the linear region
MIP_HG = 21.69    # the true high-gain MIP: th220, where two independent tools agree

TABLES = {
    "th230": ("calibration/MuonCalib_gaudi/pedestals/th230/Pedestal_TB2026CERN_run_000004_highgain.txt",
              "calibration/MuonCalib_gaudi/pedestals/th230/Pedestal_TB2026CERN_run_000004_lowgain.txt"),
    "th220": ("calibration/MuonCalib_gaudi/pedestals/th220/Pedestal_TB2026CERN_run_000th220_highgain.txt",
              "calibration/MuonCalib_gaudi/pedestals/th220/Pedestal_TB2026CERN_run_000th220_lowgain.txt"),
    "th210": ("calibration/MuonCalib_gaudi/pedestals/th210/Pedestal_TB2026CERN_run_000th210_highgain.txt",
              "calibration/MuonCalib_gaudi/pedestals/th210/Pedestal_TB2026CERN_run_000th210_lowgain.txt"),
}


def widths(path):
    """Pedestal width (the electronic noise) of every (channel, SCA)."""
    out = []
    for line in open(os.path.join(REPO, path)):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 4:
            continue
        for w in (float(x) for x in f[5::3]):
            if 0 < w < 50 and not np.isnan(w):
                out.append(w)
    return np.array(out)


wh, wl = widths(TABLES[TH][0]), widths(TABLES[TH][1])
MIP_LG = MIP_HG * K_ELEC          # what a MIP gives in the low gain: 2.17 ADC
n = min(len(wh), len(wl))
snr_h = MIP_HG / wh
snr_l = MIP_LG / wl

c = ROOT.TCanvas("c", "MIP signal to noise", 2400, 850)
c.Divide(3, 1, 0.005, 0.005)
keep = []
BLUE, RED, GREEN = ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2


def noise_panel(pad_id, gain, w, mip, colour):
    pad = c.cd(pad_id)
    pad.SetMargin(0.13, 0.05, 0.13, 0.10)
    pad.SetLogy()
    hi = max(mip * 3.0, float(np.percentile(w, 99)) * 1.2)
    h = ROOT.TH1F(f"n{pad_id}", f"{TH}: {gain} gain -- pedestal vs MIP;"
                                f"ADC;channel-SCAs", 120, 0, hi)
    h.SetDirectory(0)
    for v in w:
        h.Fill(v)
    h.SetLineColor(colour)
    h.SetLineWidth(2)
    h.SetFillColorAlpha(colour, 0.25)
    h.Draw("hist")
    keep.append(h)

    ln = ROOT.TLine(mip, 0, mip, h.GetMaximum())
    ln.SetLineColor(ROOT.kBlack)
    ln.SetLineWidth(3)
    ln.Draw()
    keep.append(ln)

    over = float((w > mip).mean())
    leg = ROOT.TLegend(0.48, 0.70, 0.95, 0.89)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.017)
    leg.AddEntry(h, f"pedestal width, median {np.median(w):.2f} ADC", "f")
    leg.AddEntry(ln, f"1 MIP = {mip:.2f} ADC", "l")
    leg.AddEntry(0, "", "")
    leg.AddEntry(0, f"median S/N = {mip / np.median(w):.1f}", "")
    leg.AddEntry(0, f"channels with pedestal > 1 MIP: {100 * over:.2f}%", "")
    leg.Draw()
    keep.append(leg)


noise_panel(1, "HIGH", wh, MIP_HG, BLUE)
noise_panel(2, "LOW", wl, MIP_LG, RED)

# ---- 3. S/N per channel, both gains on one axis --------------------------
#
# Normalised to unit area and drawn on a LINEAR y axis. The two gains hold very
# different numbers of channels at very different S/N, so a raw log-y overlay is
# unreadable -- one curve is a spike at the left edge and the other a low bump.
# What matters here is the SHAPE of each distribution and where it sits, not how
# many channels each contains, and only a normalised linear plot shows that.
pad = c.cd(3)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
snr_hists = []
for i, (s, gname, colour) in enumerate(((snr_h, "high gain", BLUE), (snr_l, "low gain", RED))):
    h = ROOT.TH1F(f"s{i}", f"{TH}: MIP / pedestal, per channel-SCA;"
                           f"MIP / pedestal width;fraction of channel-SCAs", 60, 0, 20)
    h.SetDirectory(0)
    for v in s:
        h.Fill(v)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(colour)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(colour, 0.20)
    snr_hists.append((h, gname, np.median(s)))
    keep.append(h)

ymax = max(h.GetMaximum() for h, _, _ in snr_hists)
leg = ROOT.TLegend(0.50, 0.72, 0.95, 0.89)
leg.SetBorderSize(0)
leg.SetTextSize(0.017)
for i, (h, gname, med) in enumerate(snr_hists):
    h.GetYaxis().SetRangeUser(0, ymax * 1.30)
    h.Draw("hist" if i == 0 else "hist same")
    leg.AddEntry(h, f"{gname}: median S/N = {med:.1f}", "f")

one = ROOT.TLine(1, 0, 1, ymax * 1.30)
one.SetLineColor(ROOT.kBlack)
one.SetLineWidth(3)
one.SetLineStyle(2)
one.Draw()
leg.AddEntry(one, "S/N = 1", "l")
leg.Draw()
keep += [one, leg]

# The low-gain noise floor does NOT scale with the gain; kept in the text
# summary below (the panel that drew it was redundant with panels 1-2).
scaled = wh[:n] * K_ELEC          # what the LG noise WOULD be if it scaled

out = os.path.join(OUT, f"{TH}_mip_snr.png")
c.SaveAs(out)

print(f"{TH}: is a MIP measurable in each gain?\n")
print(f"  {'':14s} {'1 MIP':>8s} {'noise':>8s} {'S/N':>7s} {'noise > 1 MIP':>15s}")
print("  " + "-" * 58)
print(f"  {'HIGH gain':14s} {MIP_HG:8.2f} {np.median(wh):8.2f} {MIP_HG / np.median(wh):7.1f} "
      f"{100 * (wh > MIP_HG).mean():14.1f}%")
print(f"  {'LOW gain':14s} {MIP_LG:8.2f} {np.median(wl):8.2f} {MIP_LG / np.median(wl):7.1f} "
      f"{100 * (wl > MIP_LG).mean():14.1f}%")
print(f"\n  If the low-gain noise scaled with the gain it would be "
      f"{np.median(scaled):.2f} ADC. It is {np.median(wl):.2f} -- "
      f"{np.median(wl) / np.median(scaled):.1f}x more.")
print(f"\nsaved {out}")
