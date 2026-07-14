"""th220 vs th230: same detector, same MIP, two trigger thresholds.

The fitted MIP MPV comes out at 21.7 ADC for th220 and 32.4 for th230 -- 49%
higher, for the SAME electronics (FeedbackCap 15, HoldDelay 110, FSPeakTime 2 in
both runs; the only thing that differs is ThresholdDAC, 220 -> 230).

The trigger explains it. The MIP histogram is only filled when the channel FIRES,
so the discriminator threshold cuts the spectrum from the left. In th230 that cut
lands at 18-22 ADC -- right on top of the MIP peak itself (~22 ADC). Half the
distribution does not exist in the data, and the peak that survives is pushed up.

This plot shows it directly. Both spectra are normalised on the TAIL (ADC > 40),
which sits far above either threshold and which neither one touches -- so if the
gains are the same, they must overlap there. They do. What does not overlap is the
peak region: th230 is MISSING the entries below ~22 ADC, exactly where its trigger
cut them away.

Usage:  PEAK_HI=50 python3 diagnostics/plot_th220_vs_th230.py [outdir] [slab]
"""
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# This one compares thresholds, so it belongs in neither folder: compare/.
_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUT, exist_ok=True)
SLAB = int(sys.argv[2]) if len(sys.argv) > 2 else 11

H = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/calib_fill_scratch/hist"
SRC = {
    "th230  (run_000004, DAC=230)": (H + "/th230/merged_th230.root", ROOT.kRed + 1),
    "th220  (run_000060, DAC=220)": (H + "/TB2026CERN_run_000060/merged_run.root", ROOT.kBlue + 1),
}
TAIL_LO, TAIL_HI = 40, 120   # region NO trigger cuts into: normalise here
PEAK_HI = int(os.environ.get("PEAK_HI", "30"))   # boundary: peak region vs above it


def combined(path, tag):
    """Whole-layer MIP spectrum, pedestal-subtracted, all 15 SCAs summed."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"ERROR: cannot open {path}")
    h = ROOT.TH1F(f"h_{tag}", "", 900, -100.5, 799.5)
    h.SetDirectory(0)
    for chip in range(16):
        for chn in range(64):
            for sca in range(15):
                mh = f.Get(f"mip_high_s{SLAB}_c{chip}_ch{chn}_sca{sca}")
                ph = f.Get(f"ped_high_s{SLAB}_c{chip}_ch{chn}_sca{sca}")
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
    f.Close()
    return h


def integral(h, lo, hi):
    return h.Integral(h.FindBin(lo), h.FindBin(hi))


hists, stats = {}, {}
for i, (label, (path, colour)) in enumerate(SRC.items()):
    h = combined(path, f"t{i}")
    # Normalise on the tail: if the gains are equal, the two must overlap there.
    tail = integral(h, TAIL_LO, TAIL_HI)
    if tail <= 0:
        raise SystemExit(f"ERROR: {label} has no tail")
    h.Scale(1.0 / tail)
    h.SetLineColor(colour)
    h.SetLineWidth(2)
    hists[label] = h
    stats[label] = {
        "peak": h.GetXaxis().GetBinCenter(h.GetMaximumBin()),
        "lo": integral(h, 0, PEAK_HI),      # peak region
        "hi": integral(h, PEAK_HI, 200),    # above it
    }

c = ROOT.TCanvas("c", "th220 vs th230", 1250, 900)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.09)

# The Y range is set by the TALLER of the two, or th220 -- the one that still has
# its whole peak -- runs off the top of the frame and you see nothing.
ymax = max(h.GetMaximum() for h in hists.values())
first = True
for label, h in hists.items():
    h.SetTitle(f"slab {SLAB}: combined MIP spectrum, normalised on the tail "
               f"(ADC {TAIL_LO}-{TAIL_HI});ADC - pedestal;entries (norm.)")
    h.GetXaxis().SetRangeUser(-10, 100)
    h.GetYaxis().SetRangeUser(0, ymax * 1.35)
    h.Draw("hist" if first else "hist same")
    first = False

keep = []
for x, colour, style, lab in ((PEAK_HI, ROOT.kGray + 2, 2, None),):
    ln = ROOT.TLine(x, 0, x, ymax * 1.05)
    ln.SetLineColor(colour)
    ln.SetLineWidth(2)
    ln.SetLineStyle(style)
    ln.Draw()
    keep.append(ln)

labels = list(hists)
a, b = stats[labels[0]], stats[labels[1]]      # a = th230, b = th220
r_lo = a["lo"] / b["lo"] if b["lo"] else float("nan")
r_hi = a["hi"] / b["hi"] if b["hi"] else float("nan")

leg = ROOT.TLegend(0.42, 0.50, 0.95, 0.90)
leg.SetBorderSize(1)
leg.SetFillColorAlpha(ROOT.kWhite, 0.85)
leg.SetTextSize(0.026)
for label, h in hists.items():
    leg.AddEntry(h, f"{label}   peak = {stats[label]['peak']:.0f} ADC", "l")
leg.AddEntry(0, "", "")
leg.AddEntry(0, "ratio th230 / th220 (tail-normalised):", "")
leg.AddEntry(0, f"   peak region, ADC < {PEAK_HI} .... {r_lo:.3f}", "")
leg.AddEntry(0, f"   above,       ADC > {PEAK_HI} .... {r_hi:.3f}", "")
leg.AddEntry(0, "", "")
leg.AddEntry(0, "Same preamp (FeedbackCap 15) in both runs.", "")
leg.AddEntry(0, "Only ThresholdDAC differs: 220 -> 230.", "")
leg.Draw()

out = os.path.join(OUT, f"th220_vs_th230_slab{SLAB}_split{PEAK_HI}.png")
c.SaveAs(out)

print(f"slab {SLAB}, normalised on the tail, ADC {TAIL_LO}-{TAIL_HI}\n")
for label in labels:
    s = stats[label]
    print(f"  {label}:  peak = {s['peak']:.0f} ADC   "
          f"integral(ADC<{PEAK_HI}) = {s['lo']:.3f}   integral(ADC>{PEAK_HI}) = {s['hi']:.3f}")
print(f"\n  th230/th220 ratio in the PEAK REGION (ADC < {PEAK_HI}): {r_lo:.3f}")
print(f"  th230/th220 ratio ABOVE it            (ADC > {PEAK_HI}): {r_hi:.3f}")
print(f"\n  Reading: normalised on the tail, th230 is MISSING "
      f"{100 * (1 - r_lo):.0f}% of the peak-region entries.")
print(f"  Those are the hits its higher trigger threshold never recorded.")
print(f"\nsaved {out}")
