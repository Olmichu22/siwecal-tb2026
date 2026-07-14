"""The fit that anchors the low gain to the high gain -- the one the event builder uses.

Above the switch (raw adc_high >= 1500) the high-gain preamp has gone non-linear -- it
is already 1.2% low there and 3.9% low by 1900 -- so the energy has to come from the
low gain. The old way divided the low-gain ADC by MIP_lg, and MIP_lg is unmeasurable: a
low-gain MIP is 2.08 ADC against ~1 ADC of noise, so the Landau fit lands at 3.2-4.5
instead. Energy read that way is under-read by that factor.

Instead, anchor the low gain to the high gain, which IS well calibrated. In the
FIDUCIAL REGION -- high-gain ADC, pedestal-subtracted, between 150 and 1200, which sits
below BOTH the switch and the onset of non-linearity -- the two branches trace out a
straight line:

    adc_low - ped_lg  =  k * (adc_high - ped_hg)  +  c

Fit it there, then INVERT it to turn a low-gain reading into the high-gain ADC the
unsaturated preamp would have given, and divide by MIP_hg:

    adc_high_equiv = (adc_low - ped_lg - c) / k          energy = adc_high_equiv / MIP_hg

MIP_lg never appears. This is what the plot shows: the line is measured only in the
shaded fiducial band, and used only beyond the switch -- an extrapolation of a
straight line into the region where the other branch has gone blind.

The line is traced by the MEDIAN of adc_low in each slice of adc_high, NOT by a
least-squares fit to the hits. Least squares on this band is pulled off the data by its
tails: it reports c = +3.13 where the medians say +1.60, and you can see it miss the
dense part of the band. Both lines are drawn, so you can.

Usage:  RUN=TB2026CERN_run_000072 TH=220 python3 diagnostics/plot_gain_anchor_fit.py
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

RUN = os.environ.get("RUN", "TB2026CERN_run_000072")
TH = os.environ.get("TH", "220")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)

N_ACQ = int(os.environ.get("N_ACQ", "8000"))
# The fiducial region, pedestal-subtracted. It must sit BELOW the switch (there is no
# point fitting where we no longer use the high gain) and below the onset of high-gain
# non-linearity, which plot_hg_linearity measures at raw adc_high ~ 1270 (1% level).
# 1200 pedestal-subtracted is raw ~1445 -- inside both bounds.
FID_LO, FID_HI = 150, 1200
SWITCH_RAW = 1500                   # AdcSaturationThreshold, on the RAW high-gain ADC
MIP_LG = 2.08
MIP_HG = 21.69

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
CELLS = NSL * NCHIP * NSCA * NCHN
XE = np.arange(-50.0, 2400.01, 12.5)
YE = np.arange(-25.0, 320.01, 0.5)
NX, NY = len(XE) - 1, len(YE) - 1
XC = 0.5 * (XE[:-1] + XE[1:])
YC = 0.5 * (YE[:-1] + YE[1:])


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
SWITCH = SWITCH_RAW - PED_HG_MED     # the switch, on the pedestal-subtracted axis

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

H = np.zeros((NX, NY), dtype=np.int64)
sca_ax = np.arange(NSCA).reshape(1, NSCA, 1)
n_ent = min(N_ACQ, chain.GetEntries())
print(f"[read] {RUN} (th{TH}): {n_ent:,} acquisitions, pedestal-subtracted switch at "
      f"{SWITCH:.0f} ADC (raw {SWITCH_RAW})")

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
        g = (ix >= 0) & (ix < NX) & (iy >= 0) & (iy < NY)
        if g.any():
            np.add.at(H, (ix[g], iy[g]), 1)
    if e and e % 2000 == 0:
        print(f"  ... {e:,}/{n_ent:,}   {int(H.sum()):,} hits")

NHIT = int(H.sum())
print(f"[done] {NHIT:,} triggering hits")


def col_median(i):
    tot = H[i].sum()
    if tot < 50:
        return np.nan
    return float(YC[np.searchsorted(np.cumsum(H[i]), 0.5 * tot)])


# The fit: a straight line through the slice MEDIANS, inside the fiducial band only.
fx, fy = [], []
for i in range(NX):
    if FID_LO <= XC[i] <= FID_HI:
        m = col_median(i)
        if not np.isnan(m):
            fx.append(XC[i])
            fy.append(m)
K, C0 = np.polyfit(np.array(fx), np.array(fy), 1)

# The bad fit, for contrast: least squares on the hits over the same region.
wi, wj = np.nonzero(H)
sel = (XC[wi] >= FID_LO) & (XC[wi] <= FID_HI)
K_LS, C_LS = np.polyfit(XC[wi][sel], YC[wj][sel], 1,
                        w=np.sqrt(H[wi, wj][sel].astype(float)))

print(f"\n  fiducial region: adc_high - ped in [{FID_LO}, {FID_HI}]  ({len(fx)} slices)")
print(f"  MEDIAN line   : adc_low - ped = {K:.4f} * (adc_high - ped) + {C0:+.2f}   (1/k = {1 / K:.2f})")
print(f"  least squares : adc_low - ped = {K_LS:.4f} * (adc_high - ped) + {C_LS:+.2f}")

c = ROOT.TCanvas("c", "anchor fit", 1700, 780)
c.Divide(2, 1, 0.005, 0.005)
keep = []
RED, GREEN, ORANGE = ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kOrange + 7


def draw(pad_id, xlo, xhi, ylo, yhi, title, zoom):
    global keep
    pad = c.cd(pad_id)
    pad.SetMargin(0.12, 0.14, 0.13, 0.09)
    pad.SetLogz()
    nbx = max(1, int((xhi - xlo) / (XE[1] - XE[0])))
    nby = max(1, int((yhi - ylo) / (YE[1] - YE[0])))
    h = ROOT.TH2F(f"p{pad_id}", f"{title};adc_high - pedestal  [ADC];adc_low - pedestal  [ADC]",
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

    # The fiducial band: the ONLY place the line is measured.
    band = ROOT.TBox(max(xlo, FID_LO), ylo, min(xhi, FID_HI), yhi)
    band.SetFillColorAlpha(GREEN, 0.10)
    band.SetLineColor(GREEN)
    band.SetLineWidth(2)
    band.Draw("l")
    keep.append(band)

    med = ROOT.TGraph()
    for i in range(NX):
        if not (xlo <= XC[i] < xhi):
            continue
        m = col_median(i)
        if not np.isnan(m) and ylo <= m < yhi:
            med.SetPoint(med.GetN(), XC[i], m)
    med.SetMarkerStyle(20)
    med.SetMarkerSize(0.8)
    med.SetMarkerColor(ROOT.kBlack)
    med.Draw("P same")
    keep.append(med)

    # The fit: SOLID where it was measured, DASHED where it is extrapolated -- and it
    # is only ever USED in the dashed part, past the switch.
    f_in = ROOT.TF1(f"fin{pad_id}", "[0]*x+[1]", max(xlo, FID_LO), min(xhi, FID_HI))
    f_in.SetParameters(K, C0)
    f_in.SetLineColor(GREEN)
    f_in.SetLineWidth(4)
    f_in.Draw("same")
    f_ex = ROOT.TF1(f"fex{pad_id}", "[0]*x+[1]", xlo, xhi)
    f_ex.SetParameters(K, C0)
    f_ex.SetLineColor(GREEN)
    f_ex.SetLineWidth(3)
    f_ex.SetLineStyle(2)
    f_ex.Draw("same")
    keep += [f_in, f_ex]

    f_ls = ROOT.TF1(f"fls{pad_id}", "[0]*x+[1]", xlo, xhi)
    f_ls.SetParameters(K_LS, C_LS)
    f_ls.SetLineColor(RED)
    f_ls.SetLineWidth(3)
    f_ls.SetLineStyle(7)
    f_ls.Draw("same")
    keep.append(f_ls)

    sw = None
    if xlo <= SWITCH <= xhi:
        sw = ROOT.TLine(SWITCH, ylo, SWITCH, yhi)
        sw.SetLineColor(ORANGE)
        sw.SetLineWidth(4)
        sw.Draw()
        keep.append(sw)

    leg = ROOT.TLegend(0.14, 0.60, 0.68, 0.90) if not zoom else ROOT.TLegend(0.14, 0.62, 0.72, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.80)
    leg.SetTextSize(0.028)
    leg.AddEntry(med, "median of adc_low per slice", "p")
    leg.AddEntry(f_in, f"FIT (medians): k = {K:.4f}, c = {C0:+.2f}   1/k = {1 / K:.2f}", "l")
    leg.AddEntry(f_ls, f"least squares on hits: k = {K_LS:.4f}, c = {C_LS:+.2f}", "l")
    if sw:
        leg.AddEntry(sw, f"switch to low gain (raw adc_high = {SWITCH_RAW})", "l")
    leg.AddEntry(band, f"fiducial region: the line is fitted HERE ONLY", "f")
    leg.Draw()
    keep.append(leg)


draw(1, -50, 2350, -20, 300,
     f"{RUN} (th{TH}): the line that anchors the low gain to the high gain", False)
draw(2, -30, 260, -8, 30,
     f"{RUN} (th{TH}): zoom on the pedestal point -- the intercept c", True)

pad = c.cd(2)
zero = ROOT.TMarker(0, 0, 29)
zero.SetMarkerColor(ROOT.kCyan + 1)
zero.SetMarkerSize(3.2)
zero.Draw()
mip = ROOT.TLine(-30, MIP_LG, 260, MIP_LG)
mip.SetLineColor(ROOT.kMagenta + 1)
mip.SetLineStyle(3)
mip.SetLineWidth(3)
mip.Draw()
lab = ROOT.TLatex()
lab.SetTextSize(0.030)
lab.SetTextColor(ROOT.kMagenta + 1)
lab.DrawLatex(150, MIP_LG + 1.0, f"1 low-gain MIP = {MIP_LG:.2f} ADC")
lab.SetTextColor(ROOT.kCyan + 1)
lab.DrawLatex(10, -5.5, "(0,0): both pedestals")
keep += [zero, mip, lab]

out = os.path.join(OUT, f"gain_anchor_fit_{RUN.split('_')[-1]}_th{TH}.png")
c.SaveAs(out)

print(f"\n  the event builder inverts this line:")
print(f"     adc_high_equiv = (adc_low - ped_lg - {C0:.2f}) / {K:.4f}")
print(f"     energy         = adc_high_equiv / MIP_hg")
print(f"  applied only past the switch (raw adc_high >= {SWITCH_RAW}). MIP_lg is never used.")
print(f"\n  the intercept c is {C0 / MIP_LG:.2f} low-gain MIPs, and shifts the energy by "
      f"c/k = {C0 / K:.1f} high-gain ADC = {C0 / K / MIP_HG:.2f} MIP per saturated hit")
print(f"\nsaved {out}")
