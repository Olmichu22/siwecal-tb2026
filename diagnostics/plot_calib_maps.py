"""Full-detector maps of the deployed calibration: pedestal mean and MIP MPV for
every channel, both gains. Reads the calibration TABLES directly (not the
histograms), so every (slab, chip, channel) that has a value is shown.

Four 2D maps on one canvas -- rows are (slab x 16 + chip), columns are channel:
  * pedestal mean, high gain     * pedestal mean, low gain
  * MIP MPV, high gain           * MIP MPV, low gain
Pedestal value per channel is the mean over its valid SCAs. Empty cells (masked /
no value) are left blank. The colour scale is clipped to the 2-98 percentile so a
few outliers do not wash out the map. (Slab 12 is the COB slab: its MIP MPV is
~1.6x the others by design, so its row reads bright -- expected, not a fault.)

Usage:  TH=220 python3 diagnostics/plot_calib_maps.py [outdir]
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(_HERE)
TH = os.environ.get("TH", "220")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
NROW = NSL * NCHIP  # 240

G = f"calibration/MuonCalib_gaudi/pedestals/th{TH}"
M = f"calibration/MuonCalib_gaudi/mips/th{TH}"
tag = "TB2026CERN_run_000004" if TH == "230" else f"TB2026CERN_run_000th{TH}"
PED_HG = f"{G}/Pedestal_{tag}_highgain.txt"
PED_LG = f"{G}/Pedestal_{tag}_lowgain.txt"
MIP_HG = f"{M}/MIP_pedestalsubmode1_{tag}_highgain.txt"
MIP_LG = f"{M}/MIP_pedestalsubmode1_{tag}_lowgain.txt"


def read_ped_mean(path):
    """Per-channel pedestal mean, averaged over its valid SCAs. NaN where none."""
    acc = np.full((NROW, NCHN), np.nan)
    if not os.path.exists(os.path.join(REPO, path)):
        return acc
    for line in open(os.path.join(REPO, path)):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 4:
            continue
        slab, chip, chn = int(f[0]), int(f[1]), int(f[2])
        if not (0 <= slab < NSL and 0 <= chip < NCHIP and 0 <= chn < NCHN):
            continue
        vals = [v for v in (float(x) for x in f[3::3]) if 0 < v < 600 and not np.isnan(v)]
        if vals:
            acc[slab * NCHIP + chip, chn] = float(np.mean(vals))
    return acc


def read_mip_mpv(path):
    """Per-channel MIP MPV (col 3). NaN where <=0 / missing."""
    acc = np.full((NROW, NCHN), np.nan)
    if not os.path.exists(os.path.join(REPO, path)):
        return acc
    for line in open(os.path.join(REPO, path)):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 5:
            continue
        slab, chip, chn = int(f[0]), int(f[1]), int(f[2])
        mpv = float(f[3])
        if 0 <= slab < NSL and 0 <= chip < NCHIP and 0 <= chn < NCHN and mpv > 0:
            acc[slab * NCHIP + chip, chn] = mpv
    return acc


PANELS = [
    ("pedestal mean, HIGH gain", read_ped_mean(PED_HG), "ADC"),
    ("pedestal mean, LOW gain", read_ped_mean(PED_LG), "ADC"),
    ("MIP MPV, HIGH gain", read_mip_mpv(MIP_HG), "ADC"),
    ("MIP MPV, LOW gain", read_mip_mpv(MIP_LG), "ADC"),
]

c = ROOT.TCanvas("c", "calib maps", 2000, 1600)
c.Divide(2, 2, 0.008, 0.012)
keep = []

for i, (title, arr, unit) in enumerate(PANELS, 1):
    pad = c.cd(i)
    pad.SetMargin(0.11, 0.14, 0.11, 0.10)
    h = ROOT.TH2F(f"m{i}", f"th{TH}: {title};channel;slab #times 16 + chip",
                  NCHN, 0, NCHN, NROW, 0, NROW)
    h.SetDirectory(0)
    finite = arr[np.isfinite(arr)]
    if finite.size:
        lo, hi = np.percentile(finite, [2, 98])
        if hi <= lo:
            hi = lo + 1
        h.GetZaxis().SetRangeUser(float(lo), float(hi))
    for r in range(NROW):
        for ch in range(NCHN):
            v = arr[r, ch]
            if np.isfinite(v):
                h.SetBinContent(ch + 1, r + 1, float(v))
    h.GetZaxis().SetTitle(unit)
    h.Draw("colz")
    n = int(np.isfinite(arr).sum())
    med = float(np.nanmedian(arr)) if n else 0.0
    lab = ROOT.TLatex()
    lab.SetNDC(True)
    lab.SetTextSize(0.030)
    lab.DrawLatex(0.12, 0.915, f"filled {n}/{NROW * NCHN}   median {med:.2f} {unit}")
    keep += [h, lab]

out = os.path.join(OUT, f"calib_maps_th{TH}.png")
c.SaveAs(out)
print(f"saved {out}")
for title, arr, unit in PANELS:
    n = int(np.isfinite(arr).sum())
    if n:
        print(f"  {title:28s}: filled {n:5d}, median {np.nanmedian(arr):7.2f} {unit}")
    else:
        print(f"  {title:28s}: NO DATA")
