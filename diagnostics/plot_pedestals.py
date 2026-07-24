"""Pedestal diagnostics: look at the spectra, don't just trust the table.

The pedestal has been argued to be sound on the strength of two number
comparisons: our table agrees with the reference tool's to 0.035 ADC (high gain)
and 0.029 ADC (low gain), and the trigger threshold sits ~11 sigma away from it so
it cannot be truncating anything. Both true. Neither is a substitute for looking.

And there is a specific reason to look: the pedestal file's own header says
"remove channels/sca with two pedestals peaks from the analysis" -- a DOUBLE-PEAKED
pedestal is a known failure mode of this detector. A Gaussian fitted to a
double-peaked distribution lands between the peaks, on a value no sample ever took.

Six panels:
  1-2. Example pedestal spectra, high and low gain, several channels overlaid --
       is each one a clean single Gaussian?
  3.   Pedestal mean per SCA. SCA0 is known to behave differently; if its pedestal
       is offset, everything filled from SCA0 inherits that offset -- and in th230,
       67% of the MIP entries ARE SCA0.
  4.   Distribution of pedestal means (should be a tight band per gain).
  5.   Distribution of pedestal widths (the noise; a channel with a huge width is
       a channel whose pedestal is not a pedestal).
  6.   Width vs mean -- outliers pop out of the blob.

Usage:  TH=th230 python3 diagnostics/plot_pedestals.py [outdir]
        TH=th220 HIST=<run>/merged_run.root python3 diagnostics/plot_pedestals.py
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
TH = os.environ.get("TH", "th230")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, TH)
os.makedirs(OUT, exist_ok=True)
HIST = os.environ.get(
    "HIST",
    f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/calib_fill_scratch/hist/{TH}/merged_{TH}.root")

FILES = [p for p in HIST.split(",") if p.strip()]
handles = []
for p in FILES:
    h = ROOT.TFile.Open(p)
    if not h or h.IsZombie():
        raise SystemExit(f"ERROR: cannot open {p}")
    handles.append(h)
print(f"[hist] {TH}: pooling {len(handles)} file(s)")


def ped_spectrum(gain, slab, chip, chn, sca):
    """The raw per-SCA pedestal histogram, summed over the input files."""
    out = None
    for f in handles:
        h = f.Get(f"ped_{gain}_s{slab}_c{chip}_ch{chn}_sca{sca}")
        if not h or h.GetEntries() == 0:
            continue
        if out is None:
            out = h.Clone(f"p_{gain}_{slab}_{chip}_{chn}_{sca}")
            out.SetDirectory(0)
        else:
            out.Add(h)
    return out


c = ROOT.TCanvas("c", "pedestals", 1800, 1200)
c.Divide(3, 2, 0.004, 0.004)
keep = []
COLS = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7]

# ---- 1-2. Example spectra, both gains ------------------------------------
CHANS = [(5, 4, 26), (8, 4, 26), (11, 4, 26), (14, 4, 26)]
for j, gain in enumerate(("high", "low")):
    pad = c.cd(j + 1)
    pad.SetMargin(0.13, 0.05, 0.13, 0.10)
    pad.SetLogy()
    leg = ROOT.TLegend(0.60, 0.62, 0.95, 0.88)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.032)
    first = True
    for i, (slab, chip, chn) in enumerate(CHANS):
        h = ped_spectrum(gain, slab, chip, chn, 0)
        if not h or h.GetEntries() < 100:
            continue
        h.SetLineColor(COLS[i % len(COLS)])
        h.SetLineWidth(2)
        h.SetTitle(f"{TH}: pedestal spectra, {gain} gain, SCA0;ADC;entries")
        h.GetXaxis().SetRangeUser(h.GetMean() - 25, h.GetMean() + 25)
        h.Draw("hist" if first else "hist same")
        first = False
        leg.AddEntry(h, f"s{slab} c{chip} ch{chn}: mean {h.GetMean():.1f}, RMS {h.GetRMS():.2f}", "l")
        keep.append(h)
    leg.Draw()
    keep.append(leg)

# ---- 3. Pedestal mean per SCA -------------------------------------------
pad = c.cd(3)
pad.SetMargin(0.14, 0.05, 0.13, 0.10)
gr = {}
for gi, gain in enumerate(("high", "low")):
    g = ROOT.TGraph()
    for sca in range(15):
        vals = []
        for slab in (2, 5, 8, 11, 14):
            for chip in (0, 4, 8, 12):
                for chn in (10, 26, 42, 58):
                    h = ped_spectrum(gain, slab, chip, chn, sca)
                    if h and h.GetEntries() > 50:
                        vals.append(h.GetMean())
        if vals:
            g.SetPoint(g.GetN(), sca, float(np.median(vals)) - (0 if gi == 0 else 0))
    g.SetMarkerStyle(20 + gi * 4)
    g.SetMarkerSize(1.3)
    g.SetMarkerColor(COLS[gi])
    g.SetLineColor(COLS[gi])
    gr[gain] = g
mg = ROOT.TMultiGraph()
mg.SetTitle(f"{TH}: median pedestal mean per SCA;SCA;pedestal mean (ADC)")
for g in gr.values():
    mg.Add(g, "LP")
mg.Draw("A")
leg3 = ROOT.TLegend(0.55, 0.75, 0.93, 0.90)
leg3.SetBorderSize(0)
leg3.SetTextSize(0.035)
for gain, g in gr.items():
    leg3.AddEntry(g, f"{gain} gain", "lp")
leg3.Draw()
keep += [mg, leg3] + list(gr.values())

# ---- 4-6. Table-level distributions --------------------------------------
def read_table(path):
    means, widths = [], []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 4:
            continue
        for m, w in zip((float(x) for x in f[3::3]), (float(x) for x in f[5::3])):
            if 0 < m < 500 and 0 < w < 50 and not (np.isnan(m) or np.isnan(w)):
                means.append(m)
                widths.append(w)
    return np.array(means), np.array(widths)


TABLES = {
    "th230": ("calibration/MuonCalib_gaudi/pedestals/th230/Pedestal_TB2026CERN_run_000004_highgain.txt",
              "calibration/MuonCalib_gaudi/pedestals/th230/Pedestal_TB2026CERN_run_000004_lowgain.txt"),
    "th220": ("calibration/MuonCalib_gaudi/pedestals/th220/Pedestal_TB2026CERN_run_000th220_highgain.txt",
              "calibration/MuonCalib_gaudi/pedestals/th220/Pedestal_TB2026CERN_run_000th220_lowgain.txt"),
    "th210": ("calibration/MuonCalib_gaudi/pedestals/th210/Pedestal_TB2026CERN_run_000th210_highgain.txt",
              "calibration/MuonCalib_gaudi/pedestals/th210/Pedestal_TB2026CERN_run_000th210_lowgain.txt"),
}
repo = os.path.dirname(_HERE)
hp, lp = (os.path.join(repo, p) for p in TABLES[TH])
mh, wh = read_table(hp)
ml, wl = read_table(lp)

for idx, (title, arrs, nb, lo, hi, xlab) in enumerate([
        ("pedestal MEAN", (mh, ml), 120, 200, 320, "pedestal mean (ADC)"),
        ("pedestal WIDTH (noise)", (wh, wl), 100, 0, 12, "pedestal width (ADC)")]):
    pad = c.cd(4 + idx)
    pad.SetMargin(0.13, 0.05, 0.13, 0.10)
    pad.SetLogy()
    leg = ROOT.TLegend(0.58, 0.70, 0.95, 0.88)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.032)
    for gi, (arr, gname) in enumerate(zip(arrs, ("high", "low"))):
        h = ROOT.TH1F(f"d{idx}{gi}", f"{TH}: {title};{xlab};SCA channels", nb, lo, hi)
        h.SetDirectory(0)
        for v in arr:
            h.Fill(v)
        h.SetLineColor(COLS[gi])
        h.SetLineWidth(2)
        h.Draw("hist" if gi == 0 else "hist same")
        leg.AddEntry(h, f"{gname}: median {np.median(arr):.2f}", "l")
        keep.append(h)
    leg.Draw()
    keep.append(leg)

pad = c.cd(6)
pad.SetMargin(0.13, 0.14, 0.13, 0.10)
pad.SetLogz()
h2 = ROOT.TH2F("h2", f"{TH}: pedestal width vs mean, HIGH gain;pedestal mean (ADC);width (ADC)",
               120, 200, 320, 100, 0, 12)
for m, w in zip(mh, wh):
    h2.Fill(m, w)
h2.Draw("colz")
keep.append(h2)

out = os.path.join(OUT, f"{TH}_pedestals.png")
c.SaveAs(out)

print(f"\n{TH} pedestals, from the tables:")
print(f"  HIGH gain: mean {np.median(mh):7.2f} ADC   width {np.median(wh):5.2f} ADC   ({len(mh):,} SCA channels)")
print(f"  LOW  gain: mean {np.median(ml):7.2f} ADC   width {np.median(wl):5.2f} ADC   ({len(ml):,} SCA channels)")
print(f"  HIGH gain width outliers (> 5 ADC): {int((wh > 5).sum()):,}  ({100 * (wh > 5).mean():.2f}%)")
print(f"  LOW  gain width outliers (> 5 ADC): {int((wl > 5).sum()):,}  ({100 * (wl > 5).mean():.2f}%)")
print(f"\nsaved {out}")
