"""Does the low-gain baseline MOVE when the channel fires? The decisive test.

A robust median-in-bins of (adc_low - ped_lg) against (adc_high - ped_hg) leaves an
intercept of +2.3 to +2.6 ADC at zero signal. Both axes are already
pedestal-subtracted, so that intercept must be zero -- a hit that deposits nothing
reads out at the pedestal in BOTH gains, by definition. It is not zero, and 2.3 ADC
is more than a whole low-gain MIP (2.08 ADC).

There are exactly two ways that can happen, and they need opposite fixes:

  (a) THE PEDESTAL IS WRONG.  A bug: the low-gain pedestal table is 2.3 ADC below
      where the baseline actually sits. Fix the calibration.

  (b) THE BASELINE MOVES.  Real physics/electronics: the pedestal is measured from
      samples that did NOT trigger, and the low-gain baseline shifts when the preamp
      fires. Then the table is RIGHT and the model E = (adc - ped) / MIP is wrong;
      the transfer function is adc_low = k * adc_high + c, with c to be calibrated.

They are told apart by asking the question the pedestal is measured on. Take one
channel and one SCA, and look at adc_low in two populations:

   N  -- the channel did NOT trigger (hitbit_high == 0). This is, by construction,
         exactly the sample the pedestal table is built from.
   T0 -- the channel DID trigger, but the high gain read out at its pedestal
         (|adc_high - ped_hg| < WIN), i.e. no signal was sampled.

Same channel, same SCA, same pedestal. If adc_low sits +2.3 ADC higher in T0 than in
N, the baseline genuinely moves when the channel fires -- (b). If the two agree, the
intercept is an artefact and the pedestal is the suspect -- (a).

AND THE CONTROL, which is what makes this a test rather than a story: the window cut
|adc_high - ped_hg| < WIN is applied to N as well. If selecting on the high gain
drags the low gain up by itself -- because their noise is correlated -- then the
control will show the same +2.3 and this whole measurement means nothing. The
control MUST come out at zero for panel 1 to be worth reading.

Panels:
  1. adc_low - ped_lg in N and in T0, pooled. The measurement.
  2. THE CONTROL: the same window cut on non-triggering samples. Must be ~0.
  3. The shift channel by channel: one number per (channel, SCA), against the MIP.
  4. The shift per SCA and per slab: is it uniform (electronics) or structured
     (a pedestal artefact -- SCA0 is known to misbehave)?

Usage:  RUN=TB2026CERN_run_000012 TH=230 python3 diagnostics/plot_pedestal_shift.py
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
N_ACQ = int(os.environ.get("N_ACQ", "3000"))
WIN = float(os.environ.get("WIN", "5"))   # "the high gain read out at its pedestal"
MIP_LG = 2.08                              # one low-gain MIP, from k * MIP_hg

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
CELLS = NSL * NCHIP * NSCA * NCHN


def read_ped(path):
    """Pedestal table -> array[slab][chip][sca][chn], 0 where the table has nothing."""
    a = np.zeros((NSL, NCHIP, NSCA, NCHN), dtype=np.float64)
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
HAVE_PED = (PH > 0) & (PL > 0)
print(f"[ped] th{TH}: {int(HAVE_PED.sum()):,} channel-SCAs with a pedestal in both gains")

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

# A HISTOGRAM per group, not a sum: every number this script quotes, at every level
# of aggregation, is a MEDIAN. A mean over a few hundred samples of a distribution
# with tails is exactly the outlier-sensitive estimator that made the earlier
# least-squares version of this measurement report +11.6 ADC where the median of the
# same data said +2.3. Sums cannot give a median, so we keep the shape.
#
# Groups are indexed (slab, chip, sca) and flattened; the 64 channels are pooled,
# since a single channel-SCA sees only a couple of T0 samples in a whole run.
EDGES = np.arange(-20.0, 25.01, 0.5)
NBIN = len(EDGES) - 1
NGRP = NSL * NCHIP * NSCA
hgrp_N = np.zeros((NGRP, NBIN), dtype=np.int64)
hgrp_T = np.zeros((NGRP, NBIN), dtype=np.int64)
hgrp_C = np.zeros((NGRP, NBIN), dtype=np.int64)

# group index for every cell of one slab's (chip, sca, chn) block
_chip_ix, _sca_ix, _ = np.meshgrid(np.arange(NCHIP), np.arange(NSCA), np.arange(NCHN),
                                   indexing="ij")
GRP_IN_SLAB = (_chip_ix * NSCA + _sca_ix)          # (chip, sca, chn) -> chip*15+sca

sca_ax = np.arange(NSCA).reshape(1, NSCA, 1)     # broadcast over (chip, sca, chn)
n_ent = min(N_ACQ, chain.GetEntries())
print(f"[read] {RUN}: {n_ent:,} acquisitions, window |adc_high - ped_hg| < {WIN:g} ADC")

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
        # Only SCAs the chip actually read out this acquisition; the rest is garbage.
        live = sca_ax < cols[slot].reshape(NCHIP, 1, 1)
        ok = live & HAVE_PED[slab] & (ah[slot] > 0) & (al[slot] > 0)
        if not ok.any():
            continue

        dl = al[slot] - PL[slab]                 # adc_low - pedestal
        dh = ah[slot] - PH[slab]                 # adc_high - pedestal
        at_ped = np.abs(dh) < WIN                # "the high gain saw nothing"
        fired = hb[slot] == 1

        m_N = ok & ~fired                        # what the pedestal is measured from
        m_T = ok & fired & at_ped                # fired, but no high-gain signal
        m_C = ok & ~fired & at_ped               # CONTROL: same cut, no trigger

        base = slab * NCHIP * NSCA
        for m, H in ((m_N, hgrp_N), (m_T, hgrp_T), (m_C, hgrp_C)):
            if not m.any():
                continue
            b = np.clip(np.searchsorted(EDGES, dl[m], side="right") - 1, 0, NBIN - 1)
            np.add.at(H, (base + GRP_IN_SLAB[m], b), 1)

    if e and e % 1000 == 0:
        print(f"  ... {e:,}/{n_ent:,}   N={int(hgrp_N.sum()):,}  "
              f"T0={int(hgrp_T.sum()):,}  control={int(hgrp_C.sum()):,}")

CTR = 0.5 * (EDGES[:-1] + EDGES[1:])


def median_of(hst):
    """Median of a binned distribution (1-D histogram of counts)."""
    tot = hst.sum()
    if tot == 0:
        return float("nan")
    return float(CTR[np.searchsorted(np.cumsum(hst), 0.5 * tot)])


def mean_of(hst):
    tot = hst.sum()
    return float((CTR * hst).sum() / tot) if tot else float("nan")


hist_N, hist_T, hist_C = hgrp_N.sum(0), hgrp_T.sum(0), hgrp_C.sum(0)
med_N, mean_N, n_N = median_of(hist_N), mean_of(hist_N), int(hist_N.sum())
med_T, mean_T, n_T = median_of(hist_T), mean_of(hist_T), int(hist_T.sum())
med_C, mean_C, n_C = median_of(hist_C), mean_of(hist_C), int(hist_C.sum())
SHIFT = med_T - med_N
CTRL = med_C - med_N

c = ROOT.TCanvas("c", "pedestal shift", 1700, 1200)
c.Divide(2, 2, 0.005, 0.005)
keep = []
BLUE, RED, GREEN, GRAY = ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kGray + 2


def to_th1(hst, name, colour):
    h = ROOT.TH1F(name, "", len(EDGES) - 1, EDGES[0], EDGES[-1])
    h.SetDirectory(0)
    for i, v in enumerate(hst):
        h.SetBinContent(i + 1, v)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(colour)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(colour, 0.20)
    keep.append(h)
    return h


def overlay(pad_id, title, pairs, note):
    pad = c.cd(pad_id)
    pad.SetMargin(0.13, 0.05, 0.13, 0.10)
    hs = [(to_th1(hst, f"p{pad_id}_{i}", col), lab, med, n)
          for i, (hst, lab, col, med, n) in enumerate(pairs)]
    ymax = max(h.GetMaximum() for h, _, _, _ in hs)
    leg = ROOT.TLegend(0.50, 0.62, 0.96, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.75)
    leg.SetTextSize(0.026)
    for i, (h, lab, med, n) in enumerate(hs):
        h.SetTitle(f"{title};adc_low - pedestal (ADC);fraction of samples")
        h.GetYaxis().SetRangeUser(0, ymax * 1.45)
        h.Draw("hist" if i == 0 else "hist same")
        leg.AddEntry(h, f"{lab}: median {med:+.2f} ADC  (N = {n:,})", "f")
    zero = ROOT.TLine(0, 0, 0, ymax * 1.45)
    zero.SetLineColor(ROOT.kBlack)
    zero.SetLineWidth(2)
    zero.SetLineStyle(2)
    zero.Draw()
    leg.AddEntry(zero, "0 = the tabulated pedestal", "l")
    leg.AddEntry(0, "", "")
    for line in note:
        leg.AddEntry(0, line, "")
    leg.Draw()
    keep.extend([zero, leg])


overlay(1, f"{RUN}: low-gain baseline, same channel and SCA",
        [(hist_N, "did NOT trigger  (= the pedestal sample)", BLUE, med_N, n_N),
         (hist_T, f"TRIGGERED, |adc_high-ped| < {WIN:g}", RED, med_T, n_T)],
        [f"shift when the channel fires: {SHIFT:+.2f} ADC",
         f"= {abs(SHIFT) / MIP_LG:.2f} low-gain MIPs"])

overlay(2, f"{RUN}: CONTROL -- same window cut, no trigger",
        [(hist_N, "did NOT trigger, all samples", BLUE, med_N, n_N),
         (hist_C, f"did NOT trigger, |adc_high-ped| < {WIN:g}", GRAY, med_C, n_C)],
        [f"bias from the window cut alone: {CTRL:+.2f} ADC",
         "this MUST be ~0, or panel 1 means nothing"])

# ---- 3. The shift, chip by chip ------------------------------------------
#
# Aggregated per CHIP, not per channel-SCA. A channel-SCA sees only a couple of T0
# samples in a whole run -- a channel that fires but samples no signal is rare by
# construction -- and a mean over 2 samples of a distribution 1.6 ADC wide is noise.
# A chip pools 64 x 15 of them, which is enough to say something.
pad = c.cd(3)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)


HT = hgrp_T.reshape(NSL, NCHIP, NSCA, NBIN)
HN = hgrp_N.reshape(NSL, NCHIP, NSCA, NBIN)


def collapse(axes):
    """MEDIAN shift per group, pooling the histograms over `axes`. Never a mean."""
    sT = HT.sum(axis=axes).reshape(-1, NBIN)
    sN = HN.sum(axis=axes).reshape(-1, NBIN)
    ok = (sT.sum(1) >= 20) & (sN.sum(1) >= 500)
    d = np.array([median_of(sT[i]) - median_of(sN[i]) for i in range(len(sT))])
    return d[ok], ok, sT.sum(1)


per_ch, _, _ = collapse((2,))         # pool the SCAs -> one median per (slab, chip)
h3 = ROOT.TH1F("h3", f"{RUN}: the shift, chip by chip;"
                     f"mean(adc_low | fired, no HG signal) - mean(adc_low | no trigger)  [ADC];"
                     f"chips", 60, -6, 14)
h3.SetDirectory(0)
for v in per_ch:
    h3.Fill(v)
h3.SetLineColor(RED)
h3.SetLineWidth(3)
h3.SetFillColorAlpha(RED, 0.22)
h3.Draw("hist")
ymax3 = h3.GetMaximum() * 1.35
h3.GetYaxis().SetRangeUser(0, ymax3)
z3 = ROOT.TLine(0, 0, 0, ymax3)
z3.SetLineColor(ROOT.kBlack)
z3.SetLineWidth(2)
z3.SetLineStyle(2)
z3.Draw()
m3 = ROOT.TLine(MIP_LG, 0, MIP_LG, ymax3)
m3.SetLineColor(GREEN)
m3.SetLineWidth(3)
m3.Draw()
leg3 = ROOT.TLegend(0.48, 0.66, 0.96, 0.89)
leg3.SetBorderSize(0)
leg3.SetTextSize(0.026)
leg3.AddEntry(h3, f"median {np.median(per_ch):+.2f} ADC  ({len(per_ch):,} chips)", "f")
leg3.AddEntry(z3, "0 = no shift (pedestal is the baseline)", "l")
leg3.AddEntry(m3, f"1 low-gain MIP = {MIP_LG:.2f} ADC", "l")
leg3.AddEntry(0, f"above zero: {100 * float((per_ch > 0).mean()):.0f}% of chips", "")
leg3.Draw()
keep += [h3, z3, m3, leg3]

# ---- 4. Is the shift structured? per SCA and per slab --------------------
pad = c.cd(4)
pad.SetMargin(0.13, 0.13, 0.13, 0.10)
g_sca, g_slab = ROOT.TGraph(), ROOT.TGraph()
sca_med, slab_med = [], []

v_sca, ok_sca, n_sca = collapse((0, 1))     # -> one median per SCA
for j, sca in enumerate(np.flatnonzero(ok_sca)):
    g_sca.SetPoint(g_sca.GetN(), int(sca), float(v_sca[j]))
    sca_med.append((int(sca), float(v_sca[j]), int(n_sca[sca])))

v_slab, ok_slab, n_slab = collapse((1, 2))  # -> one median per slab
for j, slab in enumerate(np.flatnonzero(ok_slab)):
    g_slab.SetPoint(g_slab.GetN(), int(slab), float(v_slab[j]))
    slab_med.append((int(slab), float(v_slab[j]), int(n_slab[slab])))

g_sca.SetMarkerStyle(20)
g_sca.SetMarkerSize(1.5)
g_sca.SetMarkerColor(RED)
g_sca.SetLineColor(RED)
g_slab.SetMarkerStyle(24)
g_slab.SetMarkerSize(1.5)
g_slab.SetMarkerColor(BLUE)
g_slab.SetLineColor(BLUE)
mg = ROOT.TMultiGraph()
mg.SetTitle(f"{RUN}: is the shift uniform, or structured?;SCA  /  slab;"
            f"median shift (ADC)")
mg.Add(g_sca, "LP")
mg.Add(g_slab, "LP")
mg.Draw("A")
mg.GetYaxis().SetRangeUser(min(-1.0, mg.GetYaxis().GetXmin()), max(6.0, mg.GetYaxis().GetXmax()))
z4 = ROOT.TLine(mg.GetXaxis().GetXmin(), 0, mg.GetXaxis().GetXmax(), 0)
z4.SetLineColor(ROOT.kBlack)
z4.SetLineStyle(2)
z4.Draw()
leg4 = ROOT.TLegend(0.42, 0.72, 0.96, 0.89)
leg4.SetBorderSize(0)
leg4.SetTextSize(0.026)
leg4.AddEntry(g_sca, "median shift per SCA", "lp")
leg4.AddEntry(g_slab, "median shift per slab", "lp")
leg4.AddEntry(0, "flat = an electronics effect; structured = a pedestal artefact", "")
leg4.Draw()
keep += [mg, g_sca, g_slab, z4, leg4]

out = os.path.join(OUT, f"pedestal_shift_{RUN.split('_')[-1]}_th{TH}.png")
c.SaveAs(out)

print(f"\n{RUN}, th{TH}: does the low-gain baseline move when the channel fires?\n")
print(f"  {'population':52s} {'N':>12s} {'median':>9s} {'mean':>9s}")
print("  " + "-" * 86)
print(f"  {'did NOT trigger  (what the pedestal is built from)':52s} {n_N:12,} "
      f"{med_N:9.2f} {mean_N:9.2f}")
print(f"  {f'TRIGGERED, high gain at pedestal (|dh| < {WIN:g})':52s} {n_T:12,} "
      f"{med_T:9.2f} {mean_T:9.2f}")
print(f"  {f'CONTROL: no trigger, same window cut':52s} {n_C:12,} "
      f"{med_C:9.2f} {mean_C:9.2f}")
print(f"\n  shift when the channel FIRES ...... {SHIFT:+.2f} ADC "
      f"({abs(SHIFT) / MIP_LG:.2f} low-gain MIPs)")
print(f"  bias from the window cut ALONE .... {CTRL:+.2f} ADC   "
      f"<-- must be ~0 or the measurement above is meaningless")
print(f"\n  per chip: median {np.median(per_ch):+.2f} ADC over {len(per_ch):,} chips")
print("\n  shift per SCA:")
for sca, v, n in sca_med:
    print(f"     SCA {sca:2d}   {v:+6.2f} ADC   ({n:,} T0 samples)")
print("\n  shift per slab:")
for slab, v, n in slab_med:
    print(f"     slab {slab:2d}  {v:+6.2f} ADC   ({n:,} T0 samples)")
print(f"\nsaved {out}")
