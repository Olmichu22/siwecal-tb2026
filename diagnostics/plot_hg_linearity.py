"""Where does the HIGH gain actually stop being linear? That is where the switch goes.

The switch sits at raw adc_high = 1900 on the strength of an eyeball ("the band bends
around 1950"). Looking at the anchor-fit plot it seems to bend earlier -- but that plot
draws its fitted line only out to 1500, because 1500 is where the FIT REGION ends, not
where the data bends. The end of a green line is not a measurement. This is.

Use the LOW gain as the ruler. It is linear over this whole range (it saturates far
above, at ~4x the deposits we see) while the high gain is the one that saturates. So:

  * slice in adc_low, which never lies,
  * take the MEDIAN adc_high in each slice,
  * compare it against what the anchor line says it should be:
        adc_high_expected = (adc_low - ped_lg - c) / k
    with k, c fitted in a region low enough to be unquestionably linear (150-1000),
  * the point where the measured adc_high falls below that by more than a tolerance is
    where the high gain has gone non-linear, and no hit past it should be trusted.

The residual is quoted in PERCENT of adc_high, because "10 ADC low" means something
different at 400 than at 2000.

Usage:  RUN=TB2026CERN_run_000012 TH=230 python3 diagnostics/plot_hg_linearity.py
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

RUN = os.environ.get("RUN", "TB2026CERN_run_000012")
TH = os.environ.get("TH", "230")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)

N_ACQ = int(os.environ.get("N_ACQ", "8000"))
# Fit the line where the high gain CANNOT be doubted. Deliberately well below any
# suspected turn-over, so the reference line is not itself contaminated by the bend
# it is meant to detect.
SAFE_LO, SAFE_HI = 150, 1000
TOL = [1.0, 2.0, 5.0]          # % non-linearity, the levels we quote
CUR_SWITCH = 1900

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
CELLS = NSL * NCHIP * NSCA * NCHN
# Slice in adc_low (the ruler); histogram adc_high in each slice.
LE = np.arange(-10.0, 320.01, 1.0)
HE = np.arange(-50.0, 2600.01, 5.0)
NL, NH = len(LE) - 1, len(HE) - 1
LC = 0.5 * (LE[:-1] + LE[1:])
HC = 0.5 * (HE[:-1] + HE[1:])


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
PED_HG_MED = float(np.median(PH[HAVE]))

chain = ROOT.TChain("siwecaldecoded")
for c in sorted(glob.glob(f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/"
                          f"rundata_converted_gaudi/{RUN}/chunks/chunk_*.root"))[:60]:
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

M = np.zeros((NL, NH), dtype=np.int64)     # adc_high spectrum in each adc_low slice
sca_ax = np.arange(NSCA).reshape(1, NSCA, 1)
n_ent = min(N_ACQ, chain.GetEntries())
print(f"[read] {RUN} (th{TH}): {n_ent:,} acquisitions")

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
        x = al[slot][m] - PL[slab][m]      # adc_low  - ped   (the ruler)
        y = ah[slot][m] - PH[slab][m]      # adc_high - ped   (the suspect)
        ix = np.searchsorted(LE, x, side="right") - 1
        iy = np.searchsorted(HE, y, side="right") - 1
        g = (ix >= 0) & (ix < NL) & (iy >= 0) & (iy < NH)
        if g.any():
            np.add.at(M, (ix[g], iy[g]), 1)
    if e and e % 2000 == 0:
        print(f"  ... {e:,}/{n_ent:,}   {int(M.sum()):,} hits")

print(f"[done] {int(M.sum()):,} triggering hits")


def slice_median(i, nmin=100):
    tot = M[i].sum()
    if tot < nmin:
        return np.nan, 0
    return float(HC[np.searchsorted(np.cumsum(M[i]), 0.5 * tot)]), int(tot)


# Median adc_high per adc_low slice.
lx, hy, hn = [], [], []
for i in range(NL):
    m, n = slice_median(i)
    if not np.isnan(m) and LC[i] > 0:
        lx.append(LC[i])
        hy.append(m)
        hn.append(n)
lx, hy, hn = np.array(lx), np.array(hy), np.array(hn)

# The reference line, fitted where the high gain cannot be doubted.
safe = (hy >= SAFE_LO) & (hy <= SAFE_HI)
SL, SC = np.polyfit(lx[safe], hy[safe], 1)        # adc_high = SL*adc_low + SC
print(f"\n  reference line, fitted on adc_high in [{SAFE_LO}, {SAFE_HI}] "
      f"({int(safe.sum())} slices):")
print(f"     adc_high - ped = {SL:.3f} * (adc_low - ped) {SC:+.2f}   (slope 1/k = {SL:.2f})")

expected = SL * lx + SC
resid_pct = 100.0 * (hy - expected) / np.where(expected > 0, expected, np.nan)

# Where does the high gain first fall below the line by more than each tolerance, and
# STAY there? A single noisy slice is not a turn-over; require it to hold from there on.
onset = {}
for tol in TOL:
    bad = resid_pct < -tol
    hit = None
    for j in range(len(lx)):
        if expected[j] < SAFE_HI:
            continue
        if bad[j] and bad[j:j + 8].mean() > 0.75:   # it stays bad, not one bin
            hit = expected[j]
            break
    onset[tol] = hit

print(f"\n  high-gain non-linearity: where the measured adc_high first falls BELOW the line")
print(f"  {'tolerance':>12s}   {'onset (adc_high - ped)':>24s}   {'onset (RAW adc_high)':>21s}")
print("  " + "-" * 64)
for tol in TOL:
    v = onset[tol]
    if v is None:
        print(f"  {tol:>10.0f} %   {'never reached':>24s}   {'--':>21s}")
    else:
        print(f"  {tol:>10.0f} %   {v:>24.0f}   {v + PED_HG_MED:>21.0f}")
print(f"\n  the switch is currently at RAW adc_high = {CUR_SWITCH} "
      f"(= {CUR_SWITCH - PED_HG_MED:.0f} pedestal-subtracted)")
cur_res = float(np.interp(CUR_SWITCH - PED_HG_MED, expected, resid_pct))
print(f"  at that point the high gain is already {abs(cur_res):.1f}% non-linear")
for probe in (1200, 1500, 1700):
    r = float(np.interp(probe - PED_HG_MED, expected, resid_pct))
    print(f"  at RAW adc_high = {probe}: {abs(r):.1f}% non-linear")

# ------------------------------------------------------------------ plot ----
c = ROOT.TCanvas("c", "hg linearity", 1700, 780)
c.Divide(2, 1, 0.005, 0.005)
keep = []
RED, GREEN, ORANGE, BLUE = ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kOrange + 7, ROOT.kBlue + 1

# 1. adc_high against adc_low, with the line that the LOW gain says it should follow
pad = c.cd(1)
pad.SetMargin(0.13, 0.14, 0.13, 0.09)
pad.SetLogz()
h2 = ROOT.TH2F("h2", f"{RUN} (th{TH}): the high gain against the low gain (the ruler);"
                     f"adc_low - pedestal  [ADC];adc_high - pedestal  [ADC]",
               150, 0, 300, 130, 0, 2600)
for i in range(NL):
    if not (0 <= LC[i] < 300):
        continue
    for j in np.nonzero(M[i])[0]:
        if 0 <= HC[j] < 2600:
            h2.Fill(LC[i], HC[j], float(M[i, j]))
h2.Draw("colz")
keep.append(h2)

gmed = ROOT.TGraph()
for x, y in zip(lx, hy):
    gmed.SetPoint(gmed.GetN(), x, y)
gmed.SetMarkerStyle(20)
gmed.SetMarkerSize(0.8)
gmed.SetMarkerColor(ROOT.kBlack)
gmed.Draw("P same")
fl = ROOT.TF1("fl", "[0]*x+[1]", 0, 300)
fl.SetParameters(SL, SC)
fl.SetLineColor(GREEN)
fl.SetLineWidth(4)
fl.SetLineStyle(2)
fl.Draw("same")
sw = ROOT.TLine(0, CUR_SWITCH - PED_HG_MED, 300, CUR_SWITCH - PED_HG_MED)
sw.SetLineColor(ORANGE)
sw.SetLineWidth(4)
sw.Draw()
keep += [gmed, fl, sw]
leg = ROOT.TLegend(0.16, 0.66, 0.62, 0.89)
leg.SetBorderSize(0)
leg.SetFillColorAlpha(ROOT.kWhite, 0.80)
leg.SetTextSize(0.028)
leg.AddEntry(gmed, "median adc_high per adc_low slice", "p")
leg.AddEntry(fl, f"linear, fitted on adc_high < {SAFE_HI}", "l")
leg.AddEntry(sw, f"current switch (raw {CUR_SWITCH})", "l")
leg.AddEntry(0, "the band falls BELOW the line = saturation", "")
leg.Draw()
keep.append(leg)

# 2. THE MEASUREMENT: non-linearity in %, against raw adc_high
pad = c.cd(2)
pad.SetMargin(0.13, 0.05, 0.13, 0.09)
g = ROOT.TGraph()
for xe, r in zip(expected + PED_HG_MED, resid_pct):
    if 200 < xe < 2400 and not np.isnan(r):
        g.SetPoint(g.GetN(), xe, r)
g.SetMarkerStyle(20)
g.SetMarkerSize(1.0)
g.SetMarkerColor(BLUE)
g.SetLineColor(BLUE)
g.SetLineWidth(2)
g.SetTitle(f"{RUN} (th{TH}): high-gain non-linearity;"
           f"raw adc_high  [ADC];(measured - linear) / linear   [%]")
g.Draw("ALP")
g.GetYaxis().SetRangeUser(-25, 8)
g.GetXaxis().SetLimits(200, 2400)
keep.append(g)

z = ROOT.TLine(200, 0, 2400, 0)
z.SetLineColor(ROOT.kBlack)
z.SetLineWidth(2)
z.Draw()
keep.append(z)
leg2 = ROOT.TLegend(0.16, 0.14, 0.72, 0.42)
leg2.SetBorderSize(0)
leg2.SetFillColorAlpha(ROOT.kWhite, 0.80)
leg2.SetTextSize(0.027)
leg2.AddEntry(g, "how far the high gain falls short", "lp")
for tol, colour in zip(TOL, (ROOT.kGray + 2, ORANGE, RED)):
    tl = ROOT.TLine(200, -tol, 2400, -tol)
    tl.SetLineColor(colour)
    tl.SetLineStyle(2)
    tl.SetLineWidth(2)
    tl.Draw()
    keep.append(tl)
    v = onset[tol]
    leg2.AddEntry(tl, f"-{tol:.0f}%  reached at raw adc_high = "
                      f"{('%.0f' % (v + PED_HG_MED)) if v else 'never'}", "l")
swv = ROOT.TLine(CUR_SWITCH, -25, CUR_SWITCH, 8)
swv.SetLineColor(ORANGE)
swv.SetLineWidth(4)
swv.Draw()
keep.append(swv)
leg2.AddEntry(swv, f"current switch ({CUR_SWITCH}): {abs(cur_res):.1f}% non-linear there", "l")
leg2.Draw()
keep.append(leg2)

out = os.path.join(OUT, f"hg_linearity_{RUN.split('_')[-1]}_th{TH}.png")
c.SaveAs(out)
print(f"\nsaved {out}")
