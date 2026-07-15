"""What is the ERROR on the low-gain intercept? Is run12 vs run72 even a difference?

The median-traced band gives an intercept of +1.33 ADC in run12 (th230) and +0.95 in
run72 (th220), and it is tempting to read that gap as physics -- the shift depends on
the run. But no error was ever put on either number, and without one the comparison
is empty. A 0.4 ADC gap between two numbers each uncertain by 0.5 ADC is nothing.

Three independent estimates of the error, because a single one can lie:

  1. FORMAL. Each slice median carries its own uncertainty, 1.2533 * sigma / sqrt(n)
     (the standard error of a median, not of a mean). Fit the line weighted by those
     and read the covariance. This is the smallest of the three and the least
     trustworthy: it assumes the only thing wobbling is counting noise.

  2. SUBSAMPLES. Cut the run into 4 disjoint blocks of acquisitions, trace the band in
     each, and look at the scatter of the 4 intercepts. This is empirical -- it picks
     up drift, beam structure, anything that varies through the run, not just Poisson.

  3. CHIP SPREAD. The ~60 chips each give their own intercept. The error on their
     median is 1.2533 * sigma / sqrt(N_chips). This one is inflated by any REAL
     chip-to-chip variation, so it is an upper bound.

If the two runs agree inside these, "the shift depends on the run" was noise.

Usage:  python3 diagnostics/plot_gain_line_errors.py
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

RUNS = [("TB2026CERN_run_000012", "230", ROOT.kRed + 1),
        ("TB2026CERN_run_000072", "220", ROOT.kBlue + 1)]
N_ACQ = int(os.environ.get("N_ACQ", "8000"))
NSUB = 4                       # disjoint blocks of acquisitions
LIN_LO, LIN_HI = 150, 1500
MIP_LG = 2.08

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
CELLS = NSL * NCHIP * NSCA * NCHN
XE = np.arange(-50.0, 1800.01, 12.5)
YE = np.arange(-25.0, 200.01, 0.5)
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


def col_median_and_err(col):
    """Median of a binned column, and the standard error OF THE MEDIAN."""
    n = col.sum()
    if n < 50:
        return np.nan, np.nan
    cum = np.cumsum(col)
    med = float(YC[np.searchsorted(cum, 0.5 * n)])
    # sigma from the interquartile range: robust, unlike the RMS which the tails own
    q1 = float(YC[np.searchsorted(cum, 0.25 * n)])
    q3 = float(YC[np.searchsorted(cum, 0.75 * n)])
    sigma = (q3 - q1) / 1.349
    return med, 1.2533 * sigma / np.sqrt(n)


def trace(h2, nmin=50):
    """Fit a line through the slice medians. Returns (k, c, dc_formal, npts)."""
    xs, ys, es = [], [], []
    for i in range(NX):
        if not (LIN_LO <= XC[i] <= LIN_HI) or h2[i].sum() < nmin:
            continue
        m, e = col_median_and_err(h2[i])
        if not np.isnan(m) and e > 0:
            xs.append(XC[i])
            ys.append(m)
            es.append(e)
    if len(xs) < 5:
        return None
    xs, ys, es = np.array(xs), np.array(ys), np.array(es)
    p, cov = np.polyfit(xs, ys, 1, w=1.0 / es, cov=True)
    return p[0], p[1], float(np.sqrt(cov[1, 1])), len(xs)


results = {}
for RUN, TH, colour in RUNS:
    G = f"calibration/MuonCalib_gaudi/pedestals/th{TH}"
    tag = "TB2026CERN_run_000004" if TH == "230" else "TB2026CERN_run_000th220"
    PH = read_ped(f"{G}/Pedestal_{tag}_highgain.txt")
    PL = read_ped(f"{G}/Pedestal_{tag}_lowgain.txt")
    HAVE = (PH > 0) & (PL > 0)

    chain = ROOT.TChain("siwecaldecoded")
    for c in sorted(glob.glob(f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/"
                              f"rundata_converted_gaudi/{RUN}/chunks/chunk_*.root"))[:60]:
        chain.AddFile(c)

    buf = {n: np.zeros(CELLS, dtype=np.int32)
           for n in ("adc_low", "adc_high", "hitbit_high")}
    for n, b in buf.items():
        chain.SetBranchAddress(n, b)
    ncol = np.zeros(NSL * NCHIP, dtype=np.int32)
    slbid = np.zeros(NSL, dtype=np.int32)
    chain.SetBranchAddress("nColumns", ncol)
    chain.SetBranchAddress("slboard_id", slbid)

    HS = np.zeros((NSUB, NX, NY), dtype=np.int64)     # one plane per subsample
    HC = np.zeros((NSL * NCHIP, NX, NY), dtype=np.int32)
    sca_ax = np.arange(NSCA).reshape(1, NSCA, 1)
    _chip_ix, _, _ = np.meshgrid(np.arange(NCHIP), np.arange(NSCA), np.arange(NCHN),
                                 indexing="ij")

    n_ent = min(N_ACQ, chain.GetEntries())
    print(f"[read] {RUN} (th{TH}): {n_ent:,} acquisitions, {NSUB} disjoint blocks")
    # BLOCKS of consecutive acquisitions, not interleaved: interleaving would average
    # away exactly the slow drift through the run that we are trying to detect.
    edges = np.linspace(0, n_ent, NSUB + 1).astype(int)

    for e in range(n_ent):
        chain.GetEntry(e)
        sub = int(np.searchsorted(edges, e, side="right") - 1)
        sub = min(max(sub, 0), NSUB - 1)
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
            if not g.any():
                continue
            np.add.at(HS, (sub, ix[g], iy[g]), 1)
            np.add.at(HC, (slab * NCHIP + _chip_ix[m][g], ix[g], iy[g]), 1)
        if e and e % 2000 == 0:
            print(f"  ... {e:,}/{n_ent:,}")

    H = HS.sum(0)
    k, c0, dc_formal, npts = trace(H)

    # (2) the four subsamples
    subs = []
    for s in range(NSUB):
        r = trace(HS[s], nmin=30)
        if r:
            subs.append((r[0], r[1]))
    sub_c = np.array([b for _, b in subs])
    sub_k = np.array([a for a, _ in subs])
    # spread of the MEAN of NSUB blocks: std/sqrt(n), with the (n-1) correction
    dc_sub = float(sub_c.std(ddof=1) / np.sqrt(len(sub_c))) if len(sub_c) > 1 else np.nan

    # (3) the chips
    chip_c, chip_k = [], []
    for g in range(NSL * NCHIP):
        if HC[g].sum() < 3000:
            continue
        r = trace(HC[g].astype(np.int64), nmin=30)
        if r:
            chip_k.append(r[0])
            chip_c.append(r[1])
    chip_c, chip_k = np.array(chip_c), np.array(chip_k)
    dc_chip = (1.2533 * chip_c.std(ddof=1) / np.sqrt(len(chip_c))
               if len(chip_c) > 1 else np.nan)

    results[RUN] = dict(th=TH, colour=colour, k=k, c0=c0, npts=npts,
                        dc_formal=dc_formal, sub_c=sub_c, sub_k=sub_k, dc_sub=dc_sub,
                        chip_c=chip_c, chip_k=chip_k, dc_chip=dc_chip,
                        nhit=int(H.sum()))
    print(f"  -> slope {k:.4f}, intercept {c0:+.2f} ADC   "
          f"(formal +/- {dc_formal:.2f}, subsamples +/- {dc_sub:.2f}, "
          f"chips +/- {dc_chip:.2f})")

# ---------------------------------------------------------------- plots ----
c = ROOT.TCanvas("c", "errors", 1700, 1200)
c.Divide(2, 2, 0.005, 0.005)
keep = []

# 1. The intercept of each run with all three error bars
pad = c.cd(1)
pad.SetMargin(0.15, 0.05, 0.13, 0.10)
fr = ROOT.TH2F("fr", "Low-gain intercept: is the run-to-run gap significant?;"
                     "intercept at adc_high = pedestal  [ADC];", 10, -1.0, 3.5, 4, 0, 4)
fr.GetYaxis().SetLabelSize(0)
fr.GetYaxis().SetTickLength(0)
fr.Draw()
keep.append(fr)
ytxt = []
row = 3.2
for RUN, r in results.items():
    for label, err, style in ((f"th{r['th']}  formal", r["dc_formal"], 1),
                              (f"th{r['th']}  subsamples", r["dc_sub"], 1)):
        g = ROOT.TGraphErrors(1, np.array([r["c0"]]), np.array([row]),
                              np.array([err]), np.array([0.0]))
        g.SetMarkerStyle(20)
        g.SetMarkerSize(1.6)
        g.SetMarkerColor(r["colour"])
        g.SetLineColor(r["colour"])
        g.SetLineWidth(3)
        g.Draw("P same")
        keep.append(g)
        t = ROOT.TLatex(-0.9, row + 0.13, f"{label}: {r['c0']:+.2f} #pm {err:.2f}")
        t.SetTextSize(0.033)
        t.SetTextColor(r["colour"])
        t.Draw()
        keep.append(t)
        row -= 0.8
z = ROOT.TLine(0, 0, 0, 4)
z.SetLineColor(ROOT.kCyan + 1)
z.SetLineWidth(3)
z.Draw()
mp = ROOT.TLine(MIP_LG, 0, MIP_LG, 4)
mp.SetLineColor(ROOT.kMagenta + 1)
mp.SetLineStyle(3)
mp.SetLineWidth(3)
mp.Draw()
lg = ROOT.TLegend(0.55, 0.15, 0.95, 0.30)
lg.SetBorderSize(0)
lg.SetTextSize(0.030)
lg.AddEntry(z, "0 = no offset", "l")
lg.AddEntry(mp, f"1 low-gain MIP = {MIP_LG:.2f} ADC", "l")
lg.Draw()
keep += [z, mp, lg]

# 2. The subsample intercepts, block by block
pad = c.cd(2)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
mg = ROOT.TMultiGraph()
mg.SetTitle("Intercept in each disjoint block of the run;"
            "block of acquisitions;intercept  [ADC]")
for RUN, r in results.items():
    g = ROOT.TGraph()
    for i, v in enumerate(r["sub_c"]):
        g.SetPoint(g.GetN(), i + 1, v)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.6)
    g.SetMarkerColor(r["colour"])
    g.SetLineColor(r["colour"])
    g.SetLineWidth(2)
    mg.Add(g, "LP")
    keep.append(g)
mg.Draw("A")
mg.GetYaxis().SetRangeUser(-1.5, 4.0)
lg2 = ROOT.TLegend(0.45, 0.72, 0.95, 0.89)
lg2.SetBorderSize(0)
lg2.SetTextSize(0.030)
for i, (RUN, r) in enumerate(results.items()):
    lg2.AddEntry(mg.GetListOfGraphs().At(i),
                 f"th{r['th']}: spread {r['sub_c'].std(ddof=1):.2f} ADC "
                 f"between blocks", "lp")
lg2.Draw()
z2 = ROOT.TLine(0.5, 0, len(results[RUNS[0][0]]["sub_c"]) + 0.5, 0)
z2.SetLineColor(ROOT.kCyan + 1)
z2.SetLineWidth(2)
z2.SetLineStyle(2)
z2.Draw()
keep += [mg, lg2, z2]

# 3. The chip-by-chip intercepts, both runs
pad = c.cd(3)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
lg3 = ROOT.TLegend(0.52, 0.68, 0.96, 0.89)
lg3.SetBorderSize(0)
lg3.SetTextSize(0.028)
first = True
ymax3 = 0
hs3 = []
for RUN, r in results.items():
    h = ROOT.TH1F(f"c_{r['th']}", "Intercept per chip, each traced through its own medians;"
                                  "intercept  [ADC];chips (norm.)", 40, -4, 8)
    h.SetDirectory(0)
    for v in r["chip_c"]:
        h.Fill(v)
    if h.Integral():
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(r["colour"])
    h.SetLineWidth(3)
    h.SetFillColorAlpha(r["colour"], 0.20)
    ymax3 = max(ymax3, h.GetMaximum())
    hs3.append((h, r))
    keep.append(h)
for i, (h, r) in enumerate(hs3):
    h.GetYaxis().SetRangeUser(0, ymax3 * 1.4)
    h.Draw("hist" if i == 0 else "hist same")
    lg3.AddEntry(h, f"th{r['th']}: median {np.median(r['chip_c']):+.2f} "
                    f"#pm {r['dc_chip']:.2f} ADC  ({len(r['chip_c'])} chips)", "f")
z3 = ROOT.TLine(0, 0, 0, ymax3 * 1.4)
z3.SetLineColor(ROOT.kCyan + 1)
z3.SetLineWidth(3)
z3.Draw()
lg3.AddEntry(z3, "0 = no offset", "l")
lg3.Draw()
keep += [z3, lg3]

# 4. The slope, the number that IS solid
pad = c.cd(4)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
lg4 = ROOT.TLegend(0.45, 0.68, 0.96, 0.89)
lg4.SetBorderSize(0)
lg4.SetTextSize(0.028)
ymax4 = 0
hs4 = []
for RUN, r in results.items():
    h = ROOT.TH1F(f"k_{r['th']}", "Gain ratio per chip;"
                                  "slope adc_low / adc_high;chips (norm.)", 40, 0.07, 0.13)
    h.SetDirectory(0)
    for v in r["chip_k"]:
        h.Fill(v)
    if h.Integral():
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(r["colour"])
    h.SetLineWidth(3)
    h.SetFillColorAlpha(r["colour"], 0.20)
    ymax4 = max(ymax4, h.GetMaximum())
    hs4.append((h, r))
    keep.append(h)
for i, (h, r) in enumerate(hs4):
    h.GetYaxis().SetRangeUser(0, ymax4 * 1.4)
    h.Draw("hist" if i == 0 else "hist same")
    med = np.median(r["chip_k"])
    lg4.AddEntry(h, f"th{r['th']}: k = {med:.4f}  (1/k = {1 / med:.2f})", "f")
lg4.Draw()
keep.append(lg4)

out = os.path.join(OUT, "gain_line_errors.png")
c.SaveAs(out)

print("\n" + "=" * 84)
print("  THE INTERCEPT, WITH ERRORS")
print("=" * 84)
print(f"\n  {'run':28s} {'intercept':>10s} {'formal':>9s} {'subsamp':>9s} {'chips':>9s}")
print("  " + "-" * 70)
for RUN, r in results.items():
    name = f"{RUN}  (th{r['th']})"
    print(f"  {name:28s} {r['c0']:+10.2f} {r['dc_formal']:9.2f} "
          f"{r['dc_sub']:9.2f} {r['dc_chip']:9.2f}")

a, b = [results[R[0]] for R in RUNS]
diff = a["c0"] - b["c0"]
# combine on the LARGEST of the three errors: the honest one
ea = max(a["dc_formal"], a["dc_sub"], a["dc_chip"])
eb = max(b["dc_formal"], b["dc_sub"], b["dc_chip"])
sig = diff / np.hypot(ea, eb)
print(f"\n  run12 - run72 = {diff:+.2f} ADC")
print(f"  combined error (using the LARGEST estimate for each) = {np.hypot(ea, eb):.2f} ADC")
print(f"  -> the gap is {abs(sig):.1f} sigma")
if abs(sig) < 2:
    print("\n  NOT SIGNIFICANT. The two runs are consistent with the SAME intercept.")
    w = 1 / ea**2 + 1 / eb**2
    comb = (a["c0"] / ea**2 + b["c0"] / eb**2) / w
    print(f"  Combined intercept: {comb:+.2f} +/- {1 / np.sqrt(w):.2f} ADC "
          f"({comb / MIP_LG:.2f} low-gain MIPs)")
else:
    print("\n  The gap survives. The intercept really does differ between the runs.")

print(f"\n  slope: th230 k = {np.median(a['chip_k']):.4f} (1/k = {1 / np.median(a['chip_k']):.2f}), "
      f"th220 k = {np.median(b['chip_k']):.4f} (1/k = {1 / np.median(b['chip_k']):.2f})")
print(f"\nsaved {out}")
