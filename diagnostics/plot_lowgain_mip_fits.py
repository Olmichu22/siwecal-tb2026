"""Can a MIP even be fitted in the low gain? Look at the spectra and decide.

The high-gain and low-gain MIP calibrations disagree with the electronics. Measured
on run_000013: the raw gain ratio (adc_low-ped)/(adc_high-ped) is 0.0999 -- a clean
10:1, over 21,101 hits in the linear region -- but the MIP tables imply
MIP_lg/MIP_hg = 0.13 (ours), 0.16 and 0.23 (the reference tool's th220/th230). Three
values, same detector, same electronics, spread by a factor 1.7. A well-measured
quantity does not do that.

This draws, per channel: the high-gain MIP spectrum with its Langau fit (left) and
the low-gain one with its own fit (right), on the same footing. Marked on each:
the fitted MPV, and -- on the low-gain panel -- where the MPV *must* be if the
electronics is to be believed (MIP_hg x 0.0999).

The point is not subtle once you see it: in the high gain the MIP peak stands ~32
ADC clear of the pedestal and the Langau has something to grip. In the low gain the
same peak is ~3 ADC out, buried in the pedestal shoulder, and the fit is grabbing
whatever the tail gives it.

Usage:  TH=th230 python3 diagnostics/plot_lowgain_mip_fits.py [outdir]
"""
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROOT.gInterpreter.AddIncludePath(REPO + "/gaudi_source/include")
ROOT.gInterpreter.Declare('#include "k4SiWEcalReco/PedestalMipCalib.h"')

_HERE = os.path.dirname(os.path.abspath(__file__))
TH = os.environ.get("TH", "th230")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, TH)
os.makedirs(OUT, exist_ok=True)
# HIST can be overridden: the th220 threshold-merge does not exist (its DAG failed),
# but a single th220 run's own histograms answer the same question.
HIST = os.environ.get(
    "HIST",
    f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/calib_fill_scratch/hist/{TH}/merged_{TH}.root")
K_ELEC = 0.0999  # measured raw gain ratio, run_000013 linear region

# HIST may be a COMMA-SEPARATED LIST. A single muon run does not have the ~500
# entries per channel this needs (run_000060: ~88), because it is one run, whereas
# th230's merge pools all 1,369 chunks of run_000004. Summing several th220 runs is
# exactly what the threshold-merge stage would have done had its DAG not failed.
FILES = [p for p in HIST.split(",") if p.strip()]
handles = []
for p in FILES:
    handle = ROOT.TFile.Open(p)
    if not handle or handle.IsZombie():
        raise SystemExit(f"ERROR: cannot open {p}")
    handles.append(handle)
print(f"[hist] pooling {len(handles)} file(s)")


def combined(gain, slab, chip, chn):
    """Pedestal-subtracted MIP spectrum, all 15 SCAs summed, pooled over every
    input file -- exactly what the calibrator fits (PedestalMipCalibrator::
    writeMipTable), and exactly what the threshold merge would have handed it."""
    nb, lo, hi = (900, -100.5, 799.5)
    h = ROOT.TH1F(f"h_{gain}_{slab}_{chip}_{chn}", "", nb, lo, hi)
    h.SetDirectory(0)
    for f in handles:
        for sca in range(15):
            mh = f.Get(f"mip_{gain}_s{slab}_c{chip}_ch{chn}_sca{sca}")
            ph = f.Get(f"ped_{gain}_s{slab}_c{chip}_ch{chn}_sca{sca}")
            if not mh or not ph or ph.GetEntries() == 0:
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


# Pick the four channels with the most low-gain statistics: if the fit fails
# THERE, it fails everywhere.
cands = []
for slab in (2, 5, 8, 11, 14):
    for chip in (0, 4, 8):
        for chn in (10, 26, 40, 55):
            hl = combined("low", slab, chip, chn)
            hh = combined("high", slab, chip, chn)
            if hl.Integral() > 500 and hh.Integral() > 500:
                cands.append((hl.Integral(), slab, chip, chn))
cands.sort(reverse=True)
picks = cands[:4]
if not picks:
    raise SystemExit(
        f"ERROR: no channel has >500 entries in BOTH gains across the {len(handles)} file(s) given. "
        f"Pool more runs: HIST=file1.root,file2.root,...")

c = ROOT.TCanvas("c", "low gain MIP fits", 1500, 1700)
c.Divide(2, 4, 0.002, 0.002)

keep = []
print(f"{'channel':16s} {'gain':5s} {'N':>8s} {'fit MPV':>9s} {'chi2/ndf':>9s}  note")
print("-" * 78)
for i, (_, slab, chip, chn) in enumerate(picks):
    mpv_hg = None
    for j, gain in enumerate(("high", "low")):
        pad = c.cd(2 * i + j + 1)
        pad.SetMargin(0.13, 0.05, 0.13, 0.10)
        h = combined(gain, slab, chip, chn)
        r = ROOT.k4siwecal.fitLangau(h, gain == "high")
        if gain == "high":
            mpv_hg = r.mpv
            h.GetXaxis().SetRangeUser(-20, 120)
        else:
            h.GetXaxis().SetRangeUser(-20, 30)
        h.SetTitle(f"slab {slab} chip {chip} ch {chn} -- {gain} gain;"
                   f"ADC - pedestal;hits")
        h.SetLineColor(ROOT.kBlue + 1)
        h.Draw("hist")
        keep.append(h)

        fn = h.GetFunction("langau_tmp")
        if fn:
            fn.SetLineColor(ROOT.kRed + 1)
            fn.SetLineWidth(2)
            fn.Draw("same")

        ymax = h.GetMaximum()
        ln = ROOT.TLine(r.mpv, 0, r.mpv, ymax)
        ln.SetLineColor(ROOT.kRed + 1)
        ln.SetLineWidth(2)
        ln.SetLineStyle(2)
        ln.Draw()
        keep.append(ln)

        leg = ROOT.TLegend(0.50, 0.62, 0.94, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.045)
        leg.AddEntry(ln, f"fitted MPV = {r.mpv:.2f} ADC", "l")
        leg.AddEntry(h, f"N = {h.Integral():.0f}, chi2/ndf = {r.chi2ndf:.2f}", "")

        note = ""
        if gain == "low":
            # Where the MPV has to be if the 10:1 electronics is right.
            demanded = mpv_hg * K_ELEC
            dl = ROOT.TLine(demanded, 0, demanded, ymax)
            dl.SetLineColor(ROOT.kGreen + 2)
            dl.SetLineWidth(3)
            dl.Draw()
            keep.append(dl)
            leg.AddEntry(dl, f"electronics demands {demanded:.2f} ADC", "l")
            note = f"fit is {100 * (r.mpv / demanded - 1):+.0f}% off what 10:1 demands"
        leg.Draw()
        keep.append(leg)
        label = f"s{slab} c{chip} ch{chn}"
        print(f"{label:16s} {gain:5s} {h.Integral():8.0f} {r.mpv:9.2f} {r.chi2ndf:9.2f}  {note}")

out = os.path.join(OUT, f"{TH}_lowgain_mip_fits.png")
c.SaveAs(out)
print(f"\nsaved {out}")
