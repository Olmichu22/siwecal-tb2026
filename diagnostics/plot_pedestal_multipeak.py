"""Study of MULTI-PEAKED pedestals (double, triple) across the whole detector.

A single (channel, SCA) pedestal should be one clean Gaussian. Some are not:
they carry a second (or third) baseline population a few ADC away. This scans
every (slab, chip, channel, SCA) pedestal histogram with enough statistics,
counts its peaks with TSpectrum, and asks WHERE and WHY:

  1. 2D map: fraction of that (slab, chip)'s SCA-cells that are multi-peaked --
     is the effect clustered on particular chips/slabs, or uniform?
  2. Multi-peak fraction vs SCA NUMBER -- is it tied to the SCA's place in the
     acquisition pipeline (e.g. only the early SCAs)?
  3. Peak-to-peak offset of the second population (ADC) -- one characteristic
     spacing points to one mechanism (common-mode step), a broad range to many.
  4. Peaks-per-cell histogram (1 / 2 / 3+), high vs low gain side by side.

Usage:  HIST=<grid.root> TH=210 python3 diagnostics/plot_pedestal_multipeak.py [outdir]
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

_HERE = os.path.dirname(os.path.abspath(__file__))
TH = os.environ.get("TH", "210")
HIST = os.environ["HIST"]
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)

NSL, NCHIP, NCHN, NSCA = 15, 16, 64, 15
MIN_ENTRIES = 200          # below this a peak count is not trustworthy
SIGMA, THRESH = 2.0, 0.30  # TSpectrum: peak width in bins, and height fraction


def scan(f, gain):
    """Return per-cell peak counts and, for multi-peak cells, the peak spread.
    Aggregates: nmap/mmap (cells, multi-cells per slab*chip), per-SCA counts,
    offsets, and the 1/2/3+ tally."""
    nmap = np.zeros((NSL * NCHIP,))          # cells scanned per (slab,chip)
    mmap = np.zeros((NSL * NCHIP,))          # multi-peak cells per (slab,chip)
    per_sca = np.zeros((NSCA, 2))            # [scanned, multi] per SCA
    offsets = []
    tally = {1: 0, 2: 0, 3: 0}               # 3 = "3 or more"
    sp = ROOT.TSpectrum(10)
    for slab in range(NSL):
        for chip in range(NCHIP):
            row = slab * NCHIP + chip
            for chn in range(NCHN):
                for sca in range(NSCA):
                    h = f.Get(f"ped_{gain}_s{slab}_c{chip}_ch{chn}_sca{sca}")
                    if not h or h.GetEntries() < MIN_ENTRIES:
                        continue
                    npk = sp.Search(h, SIGMA, "nobackground", THRESH)
                    nmap[row] += 1
                    per_sca[sca, 0] += 1
                    tally[min(npk, 3) if npk >= 1 else 1] += 1
                    if npk >= 2:
                        mmap[row] += 1
                        per_sca[sca, 1] += 1
                        xs = sorted(sp.GetPositionX()[i] for i in range(npk))
                        offsets.append(xs[-1] - xs[0])
    return nmap, mmap, per_sca, np.array(offsets), tally


f = ROOT.TFile.Open(HIST)
if not f or f.IsZombie():
    raise SystemExit(f"cannot open {HIST}")

# The map + per-SCA + offsets come from the low gain (where the user sees it);
# the peaks-per-cell tally is shown for both gains.
nmap, mmap, per_sca, offsets, tally_lo = scan(f, "low")
_, _, _, _, tally_hi = scan(f, "high")

tot = nmap.sum()
multi = mmap.sum()
print(f"th{TH} LOW gain: {int(tot):,} (channel,SCA) pedestals with >={MIN_ENTRIES} entries")
print(f"  multi-peak (>=2): {int(multi):,} ({100 * multi / max(tot, 1):.1f}%)   "
      f"[1: {tally_lo[1]:,}  2: {tally_lo[2]:,}  3+: {tally_lo[3]:,}]")
if len(offsets):
    print(f"  peak-to-peak offset: median {np.median(offsets):.1f} ADC, "
          f"IQR [{np.percentile(offsets, 25):.1f}, {np.percentile(offsets, 75):.1f}]")

c = ROOT.TCanvas("c", "pedestal multipeak", 2000, 1600)
c.Divide(2, 2, 0.008, 0.012)
keep = []

# 1. map of multi-peak fraction per (slab,chip)
c.cd(1).SetMargin(0.11, 0.14, 0.11, 0.10)
c.cd(1).SetRightMargin(0.16)
hm = ROOT.TH2F("hm", f"th{TH}: multi-peak fraction of pedestals (low gain);"
                     f"chip;slab", NCHIP, 0, NCHIP, NSL, 0, NSL)
hm.SetDirectory(0)
for slab in range(NSL):
    for chip in range(NCHIP):
        row = slab * NCHIP + chip
        if nmap[row] > 0:
            hm.SetBinContent(chip + 1, slab + 1, mmap[row] / nmap[row])
hm.GetZaxis().SetTitle("fraction")
hm.Draw("colz")
keep.append(hm)

# 2. multi-peak fraction vs SCA number
c.cd(2).SetMargin(0.12, 0.05, 0.12, 0.10)
hs = ROOT.TH1F("hs", f"th{TH}: multi-peak fraction vs SCA number;SCA;fraction multi-peak",
               NSCA, 0, NSCA)
hs.SetDirectory(0)
for sca in range(NSCA):
    if per_sca[sca, 0] > 0:
        hs.SetBinContent(sca + 1, per_sca[sca, 1] / per_sca[sca, 0])
hs.SetLineColor(ROOT.kAzure + 1)
hs.SetLineWidth(3)
hs.SetFillColorAlpha(ROOT.kAzure + 1, 0.3)
hs.SetMinimum(0)
hs.Draw("hist")
keep.append(hs)

# 3. peak-to-peak offset distribution
c.cd(3).SetMargin(0.12, 0.05, 0.12, 0.10)
if len(offsets):
    ho = ROOT.TH1F("ho", f"th{TH}: 2nd-population offset;peak-to-peak spacing  [ADC];multi-peak cells",
                   60, 0, float(np.percentile(offsets, 99)) + 2)
    ho.SetDirectory(0)
    for v in offsets:
        ho.Fill(v)
    ho.SetLineColor(ROOT.kRed + 1)
    ho.SetLineWidth(3)
    ho.SetFillColorAlpha(ROOT.kRed + 1, 0.25)
    ho.Draw("hist")
    keep.append(ho)

# 4. peaks-per-cell tally, high vs low gain
c.cd(4).SetMargin(0.12, 0.05, 0.12, 0.10)
ht = ROOT.TH1F("ht", f"th{TH}: peaks per pedestal cell;peaks;fraction of cells", 3, 0.5, 3.5)
hh = ROOT.TH1F("hh", "", 3, 0.5, 3.5)
for h, tally, col in ((ht, tally_lo, ROOT.kAzure + 1), (hh, tally_hi, ROOT.kOrange + 7)):
    s = sum(tally.values()) or 1
    for k in (1, 2, 3):
        h.SetBinContent(k, tally[k] / s)
    h.SetLineColor(col)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(col, 0.25)
ht.GetXaxis().SetBinLabel(1, "1")
ht.GetXaxis().SetBinLabel(2, "2")
ht.GetXaxis().SetBinLabel(3, "3+")
ht.SetMaximum(1.1)
ht.Draw("hist")
hh.Draw("hist same")
leg = ROOT.TLegend(0.5, 0.75, 0.93, 0.88)
leg.SetBorderSize(0)
leg.SetTextSize(0.03)
leg.AddEntry(ht, "low gain", "f")
leg.AddEntry(hh, "high gain", "f")
leg.Draw()
keep += [ht, hh, leg]

out = os.path.join(OUT, f"pedestal_multipeak_th{TH}.png")
c.SaveAs(out)
print(f"saved {out}")
