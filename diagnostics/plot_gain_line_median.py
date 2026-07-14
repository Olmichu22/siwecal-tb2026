"""The low-gain / high-gain line, traced by MEDIANS instead of a least-squares fit.

A least-squares fit of (adc_low - ped_lg) against (adc_high - ped_hg) reports an
intercept of +11.6 ADC. That fit is worthless and you can see it is worthless: the
line it draws does not pass through the densest part of the band. Least squares
minimises the SQUARE of the residual, so a handful of hits far off the band pull it
harder than the hundreds of thousands of hits sitting on it.

So do not fit the hits. Trace the band:

  1. slice the plane in adc_high,
  2. take the MEDIAN of adc_low in each slice -- a single outlier cannot move it,
  3. fit a straight line THROUGH THE MEDIANS in the region where the high gain is
     linear, and read the intercept.

The two lines are drawn on top of each other here. The median line is the one that
follows the band. Its intercept is the number worth quoting.

Panels:
  1. The plane, with the medians and both lines. Where does each line actually go?
  2. Zoom on the origin: the pedestal point (0,0) is where the band MUST cross if
     the pedestals are right. The median line misses it by the intercept.
  3. The residual of each median slice about the median line -- is a straight line
     even the right model?
  4. The intercept chip by chip, each one fitted through ITS OWN medians.

Usage:  RUN=TB2026CERN_run_000012 TH=230 python3 diagnostics/plot_gain_line_median.py
"""
import glob
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(_HERE)
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUT, exist_ok=True)

RUN = os.environ.get("RUN", "TB2026CERN_run_000012")
TH = os.environ.get("TH", "230")
N_ACQ = int(os.environ.get("N_ACQ", "6000"))
LIN_LO, LIN_HI = 150, 1500     # where the high gain is unambiguously linear
MIP_LG = 2.08

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
CELLS = NSL * NCHIP * NSCA * NCHN


def read_ped(path):
    a = np.zeros((NSL, NCHIP, NSCA, NCHN))
    for line in open(os.path.join(REPO, path)):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 4:
            continue
        slab, chip, chn = int(f[0]), int(f[1]), int(f[2])
        if not (0 <= slab < NSL and 0 <= chip < NCHIP and 0 <= chn < NCHN):
            continue
        for sca, m in enumerate(float(x) for x in f[3::3]):
            if sca < NSCA and 0 < m < 500 and not np.isnan(m):
                a[slab, chip, sca, chn] = m
    return a


G = f"calibration/MuonCalib_gaudi/pedestals/th{TH}"
tag = "TB2026CERN_run_000004" if TH == "230" else "TB2026CERN_run_000th220"
PH = read_ped(f"{G}/Pedestal_{tag}_highgain.txt")
PL = read_ped(f"{G}/Pedestal_{tag}_lowgain.txt")
HAVE = (PH > 0) & (PL > 0)

chain = ROOT.TChain("siwecaldecoded")
for c in sorted(glob.glob(f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/"
                          f"rundata_converted_gaudi/{RUN}/chunks/chunk_*.root"))[:40]:
    chain.AddFile(c)
if chain.GetEntries() == 0:
    raise SystemExit(f"ERROR: no decoded chunks for {RUN}")

buf = {n: np.zeros(CELLS, dtype=np.int32) for n in ("adc_low", "adc_high", "hitbit_high")}
for n, b in buf.items():
    chain.SetBranchAddress(n, b)
ncol = np.zeros(NSL * NCHIP, dtype=np.int32)
slbid = np.zeros(NSL, dtype=np.int32)
chain.SetBranchAddress("nColumns", ncol)
chain.SetBranchAddress("slboard_id", slbid)

# 2-D histogram of the plane, plus the same thing per chip. Histograms, not point
# lists: 90 M cells will not fit in memory, and a histogram is all a median needs.
XE = np.arange(-50.0, 1800.01, 12.5)          # adc_high - ped
YE = np.arange(-25.0, 200.01, 0.5)            # adc_low  - ped
NX, NY = len(XE) - 1, len(YE) - 1
H = np.zeros((NX, NY), dtype=np.int64)
HC = np.zeros((NSL * NCHIP, NX, NY), dtype=np.int32)

sca_ax = np.arange(NSCA).reshape(1, NSCA, 1)
n_ent = min(N_ACQ, chain.GetEntries())
print(f"[read] {RUN} (th{TH}): {n_ent:,} acquisitions")

_chip_ix, _, _ = np.meshgrid(np.arange(NCHIP), np.arange(NSCA), np.arange(NCHN),
                             indexing="ij")

for e in range(n_ent):
    chain.GetEntry(e)
    ah = buf["adc_high"].reshape(NSL, NCHIP, NSCA, NCHN)
    al = buf["adc_low"].reshape(NSL, NCHIP, NSCA, NCHN)
    hb = buf["hitbit_high"].reshape(NSL, NCHIP, NSCA, NCHN)
    cols = ncol.reshape(NSL, NCHIP)
    for slot in range(NSL):
        slab = int(slbid[slot])
        if not (0 <= slab < NSL):
            continue
        live = sca_ax < cols[slot].reshape(NCHIP, 1, 1)
        m = live & HAVE[slab] & (hb[slot] == 1) & (ah[slot] > 0) & (al[slot] > 0)
        if not m.any():
            continue
        x = ah[slot][m] - PH[slab][m]
        y = al[slot][m] - PL[slab][m]
        ix = np.searchsorted(XE, x, side="right") - 1
        iy = np.searchsorted(YE, y, side="right") - 1
        good = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        if not good.any():
            continue
        np.add.at(H, (ix[good], iy[good]), 1)
        g = slab * NCHIP + _chip_ix[m][good]
        np.add.at(HC, (g, ix[good], iy[good]), 1)
    if e and e % 1000 == 0:
        print(f"  ... {e:,}/{n_ent:,}   {int(H.sum()):,} hits")

NHIT = int(H.sum())
print(f"[done] {NHIT:,} triggering hits with a pedestal in both gains")

XC = 0.5 * (XE[:-1] + XE[1:])
YC = 0.5 * (YE[:-1] + YE[1:])


def col_median(col):
    tot = col.sum()
    if tot < 50:
        return np.nan
    return float(YC[np.searchsorted(np.cumsum(col), 0.5 * tot)])


def median_line(h2, lo=LIN_LO, hi=LIN_HI, nmin=50):
    """Slice, take the median of each slice, fit a line THROUGH THE MEDIANS."""
    xs, ys = [], []
    for i in range(NX):
        if not (lo <= XC[i] <= hi) or h2[i].sum() < nmin:
            continue
        m = col_median(h2[i])
        if not np.isnan(m):
            xs.append(XC[i])
            ys.append(m)
    if len(xs) < 5:
        return None
    xs, ys = np.array(xs), np.array(ys)
    k, c0 = np.polyfit(xs, ys, 1)
    return k, c0, xs, ys


fit = median_line(H)
if fit is None:
    raise SystemExit("ERROR: not enough populated slices to trace the band")
K, C0, MX, MY = fit

# The bad fit, for contrast: least squares on the HITS, weighted by the counts. This
# is the one that reports +11.6 and draws a line through empty space.
wi, wj = np.nonzero(H)
sel = (XC[wi] >= LIN_LO) & (XC[wi] <= LIN_HI)
w = H[wi, wj][sel].astype(float)
K_LS, C_LS = np.polyfit(XC[wi][sel], YC[wj][sel], 1, w=np.sqrt(w))

print(f"\n  median line (through the slice medians): slope {K:.4f}, intercept {C0:+.2f} ADC")
print(f"  least squares on the hits              : slope {K_LS:.4f}, intercept {C_LS:+.2f} ADC")

c = ROOT.TCanvas("c", "gain line", 1700, 1200)
c.Divide(2, 2, 0.005, 0.005)
keep = []
RED, GREEN, BLACK = ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kBlack


def plane(pad_id, xlo, xhi, ylo, yhi, title, with_mip):
    global keep
    pad = c.cd(pad_id)
    pad.SetMargin(0.13, 0.14, 0.13, 0.10)
    pad.SetLogz()
    nbx = max(1, int((xhi - xlo) / (XE[1] - XE[0])))
    nby = max(1, int((yhi - ylo) / (YE[1] - YE[0])))
    h = ROOT.TH2F(f"pl{pad_id}", f"{title};adc_high - pedestal;adc_low - pedestal",
                  nbx, xlo, xhi, nby, ylo, yhi)
    h.SetDirectory(0)
    for i in range(NX):
        if not (xlo <= XC[i] < xhi):
            continue
        for j in np.nonzero(H[i])[0]:
            if ylo <= YC[j] < yhi:
                h.Fill(XC[i], YC[j], float(H[i, j]))
    h.Draw("colz")
    keep.append(h)

    med = ROOT.TGraph()
    for x, y in zip(*[np.array(a) for a in ([XC[i] for i in range(NX) if H[i].sum() >= 50],
                                            [col_median(H[i]) for i in range(NX)
                                             if H[i].sum() >= 50])]):
        if xlo <= x < xhi and not np.isnan(y):
            med.SetPoint(med.GetN(), x, y)
    med.SetMarkerStyle(20)
    med.SetMarkerSize(0.9)
    med.SetMarkerColor(BLACK)
    med.Draw("P same")
    keep.append(med)

    f1 = ROOT.TF1(f"f{pad_id}", "[0]*x+[1]", xlo, xhi)
    f1.SetParameters(K, C0)
    f1.SetLineColor(GREEN)
    f1.SetLineWidth(3)
    f1.Draw("same")
    f2 = ROOT.TF1(f"g{pad_id}", "[0]*x+[1]", xlo, xhi)
    f2.SetParameters(K_LS, C_LS)
    f2.SetLineColor(RED)
    f2.SetLineWidth(3)
    f2.SetLineStyle(2)
    f2.Draw("same")
    z = ROOT.TMarker(0, 0, 29)
    z.SetMarkerColor(ROOT.kCyan + 1)
    z.SetMarkerSize(3.0)
    z.Draw()
    keep += [f1, f2, z]

    leg = ROOT.TLegend(0.16, 0.62, 0.66, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.78)
    leg.SetTextSize(0.026)
    leg.AddEntry(med, "median of adc_low in each slice", "p")
    leg.AddEntry(f1, f"line through the MEDIANS: {K:.4f}x {C0:+.2f}", "l")
    leg.AddEntry(f2, f"least squares on the hits: {K_LS:.4f}x {C_LS:+.2f}", "l")
    leg.AddEntry(z, "(0,0): where the band must cross", "p")
    if with_mip:
        ml = ROOT.TLine(xlo, MIP_LG, xhi, MIP_LG)
        ml.SetLineColor(ROOT.kMagenta + 1)
        ml.SetLineStyle(3)
        ml.SetLineWidth(2)
        ml.Draw()
        leg.AddEntry(ml, f"1 low-gain MIP = {MIP_LG:.2f} ADC", "l")
        keep.append(ml)
    leg.Draw()
    keep.append(leg)


plane(1, -50, 1700, -25, 175,
      f"{RUN} (th{TH}): the band, and two ways to fit it", False)
plane(2, -50, 300, -8, 32,
      f"{RUN} (th{TH}): zoom on the pedestal point", True)

# ---- 3. Residual of each median about the median line --------------------
pad = c.cd(3)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
res = ROOT.TGraph()
for x in XC:
    i = int(np.searchsorted(XC, x))
    if i >= NX or H[i].sum() < 50:
        continue
    m = col_median(H[i])
    if not np.isnan(m):
        res.SetPoint(res.GetN(), x, m - (K * x + C0))
res.SetMarkerStyle(20)
res.SetMarkerSize(0.9)
res.SetMarkerColor(GREEN)
res.SetTitle(f"{RUN}: residual of each slice median about the median line;"
             f"adc_high - pedestal;median(adc_low) - line  [ADC]")
res.Draw("AP")
res.GetYaxis().SetRangeUser(-6, 6)
z3 = ROOT.TLine(res.GetXaxis().GetXmin(), 0, res.GetXaxis().GetXmax(), 0)
z3.SetLineColor(BLACK)
z3.SetLineStyle(2)
z3.SetLineWidth(2)
z3.Draw()
b1 = ROOT.TLine(LIN_LO, -6, LIN_LO, 6)
b2 = ROOT.TLine(LIN_HI, -6, LIN_HI, 6)
for b in (b1, b2):
    b.SetLineColor(ROOT.kGray + 2)
    b.SetLineStyle(2)
    b.Draw()
leg3 = ROOT.TLegend(0.16, 0.74, 0.70, 0.89)
leg3.SetBorderSize(0)
leg3.SetTextSize(0.028)
leg3.AddEntry(res, "slice medians about the line", "p")
leg3.AddEntry(b1, f"fit region: {LIN_LO}-{LIN_HI} ADC", "l")
leg3.AddEntry(0, "flat = a straight line is the right model", "")
leg3.Draw()
keep += [res, z3, b1, b2, leg3]

# ---- 4. The intercept, chip by chip, each from its own medians -----------
pad = c.cd(4)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
ints, slopes = [], []
for g in range(NSL * NCHIP):
    if HC[g].sum() < 2000:
        continue
    saved = H.copy()
    r = None
    xs, ys = [], []
    for i in range(NX):
        if not (LIN_LO <= XC[i] <= LIN_HI) or HC[g, i].sum() < 30:
            continue
        tot = HC[g, i].sum()
        m = float(YC[np.searchsorted(np.cumsum(HC[g, i]), 0.5 * tot)])
        xs.append(XC[i])
        ys.append(m)
    if len(xs) >= 5:
        kk, cc = np.polyfit(np.array(xs), np.array(ys), 1)
        slopes.append(kk)
        ints.append(cc)
ints, slopes = np.array(ints), np.array(slopes)
h4 = ROOT.TH1F("h4", f"{RUN} (th{TH}): intercept per chip, each fitted through ITS OWN medians;"
                     f"intercept at adc_high = pedestal  [ADC];chips", 50, -6, 12)
h4.SetDirectory(0)
for v in ints:
    h4.Fill(v)
h4.SetLineColor(GREEN)
h4.SetLineWidth(3)
h4.SetFillColorAlpha(GREEN, 0.22)
h4.Draw("hist")
ym = h4.GetMaximum() * 1.4
h4.GetYaxis().SetRangeUser(0, ym)
z4 = ROOT.TLine(0, 0, 0, ym)
z4.SetLineColor(ROOT.kCyan + 1)
z4.SetLineWidth(3)
z4.Draw()
m4 = ROOT.TLine(MIP_LG, 0, MIP_LG, ym)
m4.SetLineColor(ROOT.kMagenta + 1)
m4.SetLineStyle(3)
m4.SetLineWidth(3)
m4.Draw()
leg4 = ROOT.TLegend(0.46, 0.64, 0.96, 0.89)
leg4.SetBorderSize(0)
leg4.SetTextSize(0.026)
leg4.AddEntry(h4, f"median {np.median(ints):+.2f} ADC  ({len(ints):,} chips)", "f")
leg4.AddEntry(z4, "0 = pedestals consistent between gains", "l")
leg4.AddEntry(m4, f"1 low-gain MIP = {MIP_LG:.2f} ADC", "l")
leg4.AddEntry(0, f"median slope k = {np.median(slopes):.4f}  (1/k = {1 / np.median(slopes):.1f})", "")
leg4.Draw()
keep += [h4, z4, m4, leg4]

out = os.path.join(OUT, f"gain_line_median_{RUN.split('_')[-1]}_th{TH}.png")
c.SaveAs(out)

print(f"\n{RUN} (th{TH}): tracing the band with medians, {NHIT:,} hits\n")
print(f"  {'estimator':40s} {'slope':>9s} {'1/slope':>9s} {'intercept':>11s}")
print("  " + "-" * 74)
print(f"  {'line through the SLICE MEDIANS':40s} {K:9.4f} {1 / K:9.2f} {C0:+11.2f}")
print(f"  {'least squares on the hits (biased)':40s} {K_LS:9.4f} {1 / K_LS:9.2f} {C_LS:+11.2f}")
print(f"\n  per chip, each from its own medians: intercept {np.median(ints):+.2f} ADC, "
      f"slope {np.median(slopes):.4f}  ({len(ints):,} chips)")
print(f"  the intercept is {abs(np.median(ints)) / MIP_LG:.2f} low-gain MIPs")
print(f"\nsaved {out}")
