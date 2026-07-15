"""Is run12's extra 19% a GAIN difference, or PILE-UP from a hotter beam?

run_000012 (th230) and run_000072 (th220) deposit very different amounts of raw ADC:

    nhit_chan   434.5  vs  418.9    ratio 1.037    <- almost the same number of hits
    sum_hg    67446    vs 56457     ratio 1.195    <- but 19% more ADC

Same FeedbackCap (15), HoldDelay (110), FSPeakTime (2). So the electronics are set up
identically, and yet each hit in run12 carries ~15% more charge. Two stories fit:

  GAIN.     The detector response really drifted between the two periods (18 June vs
            20-24 June). If so, one MIP table cannot serve both, and the whole
            "th220 MIP for everyone" decision is under-corrected.

  PILE-UP.  The beam spill rate differed. Two particles landing in the SAME CELL in the
            SAME BCID have their charges SUMMED into one hit. That inflates the ADC per
            hit WITHOUT creating new hits -- which is exactly the signature above.

They make OPPOSITE predictions about the per-hit ADC spectrum, and that is the test:

  GAIN    -> the whole spectrum is STRETCHED by 1.19. THE PEAK MOVES.
  PILE-UP -> the peak stays exactly where it is (a single particle still deposits what
             a single particle deposits) and the HIGH TAIL grows.

A peak that moves cannot be pile-up. A peak that does not move cannot be gain.

Panels:
  1. Per-hit ADC spectrum, both runs, normalised. Where is the peak?
  2. The same, zoomed on the peak, with each peak marked.
  3. Occupancy and rate: hits per acquisition, SCAs used, events per acquisition.
     Pile-up requires a hotter beam -- if run12 is not busier, pile-up is out.
  4. Ratio of the two spectra, bin by bin. Flat = same shape (rate/normalisation);
     rising = run12's hits are systematically bigger (gain).

Usage:  python3 diagnostics/plot_rate_vs_gain.py
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

CONV = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"
RUNS = [("TB2026CERN_run_000012", "230", ROOT.kRed + 1),
        ("TB2026CERN_run_000072", "220", ROOT.kBlue + 1)]
N_ACQ = int(os.environ.get("N_ACQ", "8000"))

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
CELLS = NSL * NCHIP * NSCA * NCHN


def read_ped(th):
    tag = "TB2026CERN_run_000004" if th == "230" else "TB2026CERN_run_000th220"
    p = f"calibration/MuonCalib_gaudi/pedestals/th{th}/Pedestal_{tag}_highgain.txt"
    a = np.zeros((NSL, NCHIP, NSCA, NCHN))
    for line in open(os.path.join(REPO, p)):
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


EDGES = np.arange(0.0, 400.01, 2.0)
CTR = 0.5 * (EDGES[:-1] + EDGES[1:])
data = {}

for run, th, colour in RUNS:
    PH = read_ped(th)
    HAVE = PH > 0
    chain = ROOT.TChain("siwecaldecoded")
    for c in sorted(glob.glob(f"{CONV}/{run}/chunks/chunk_*.root"))[:60]:
        chain.AddFile(c)
    if chain.GetEntries() == 0:
        raise SystemExit(f"ERROR: no chunks for {run}")

    ah = np.zeros(CELLS, dtype=np.int32)
    hb = np.zeros(CELLS, dtype=np.int32)
    ncol = np.zeros(NSL * NCHIP, dtype=np.int32)
    slbid = np.zeros(NSL, dtype=np.int32)
    bcid = np.zeros(NSL * NCHIP * NSCA, dtype=np.int32)
    chain.SetBranchAddress("adc_high", ah)
    chain.SetBranchAddress("hitbit_high", hb)
    chain.SetBranchAddress("nColumns", ncol)
    chain.SetBranchAddress("slboard_id", slbid)

    hist = np.zeros(len(EDGES) - 1)
    hits_per_acq, scas_used, acqwin, delay = [], [], [], []
    sca_ax = np.arange(NSCA).reshape(1, NSCA, 1)
    n_ent = min(N_ACQ, chain.GetEntries())
    print(f"[read] {run} (th{th}): {n_ent:,} acquisitions")

    for e in range(n_ent):
        chain.GetEntry(e)
        acqwin.append(float(chain.acqWindowMs))
        delay.append(float(chain.delayBetweenCycleMs))
        A = ah.reshape(NSL, NCHIP, NSCA, NCHN)
        B = hb.reshape(NSL, NCHIP, NSCA, NCHN)
        cols = ncol.reshape(NSL, NCHIP)
        nh = 0
        for slot in range(NSL):
            slab = int(slbid[slot])
            if not (0 <= slab < NSL):
                continue
            live = sca_ax < cols[slot].reshape(NCHIP, 1, 1)
            m = live & HAVE[slab] & (B[slot] == 1) & (A[slot] > 0)
            if not m.any():
                continue
            v = A[slot][m] - PH[slab][m]
            hist += np.histogram(v, bins=EDGES)[0]
            nh += int(m.sum())
        hits_per_acq.append(nh)
        # SCAs actually read out, averaged over the chips that read anything
        cc = cols[cols > 0]
        scas_used.append(float(cc.mean()) if cc.size else 0.0)
        if e and e % 2000 == 0:
            print(f"  ... {e:,}/{n_ent:,}")

    data[run] = dict(th=th, colour=colour, hist=hist,
                     hits=np.array(hits_per_acq), scas=np.array(scas_used),
                     acqwin=float(np.median(acqwin)), delay=float(np.median(delay)),
                     n_acq=n_ent)
    d = data[run]
    print(f"  hits/acquisition: {d['hits'].mean():.1f}   SCAs used: {d['scas'].mean():.2f}"
          f"   acqWindow {d['acqwin']:.1f} ms   delay {d['delay']:.1f} ms")

a, b = data[RUNS[0][0]], data[RUNS[1][0]]     # a = run12 (th230), b = run72 (th220)


def peak_of(h, lo=10.0):
    """Robust peak: smooth, then take the maximum above `lo` (skip the trigger edge)."""
    k = np.ones(5) / 5.0
    s = np.convolve(h, k, mode="same")
    m = CTR > lo
    return float(CTR[m][np.argmax(s[m])])


def median_of(h, lo=0.0):
    m = CTR > lo
    hh, cc = h[m], CTR[m]
    tot = hh.sum()
    return float(cc[np.searchsorted(np.cumsum(hh), 0.5 * tot)]) if tot else np.nan


# Normalise on the number of hits, so we compare SHAPES.
na, nb = a["hist"].sum(), b["hist"].sum()
ha, hb_ = a["hist"] / na, b["hist"] / nb
pk_a, pk_b = peak_of(a["hist"]), peak_of(b["hist"])
md_a, md_b = median_of(a["hist"]), median_of(b["hist"])
mean_a = float((CTR * a["hist"]).sum() / na)
mean_b = float((CTR * b["hist"]).sum() / nb)

print("\n" + "=" * 74)
print("  GAIN or PILE-UP? What the per-hit ADC spectrum says")
print("=" * 74)
print(f"\n  {'':22s} {'run12 (th230)':>14s} {'run72 (th220)':>14s} {'ratio':>9s}")
print("  " + "-" * 64)
print(f"  {'PEAK of the spectrum':22s} {pk_a:14.1f} {pk_b:14.1f} {pk_a / pk_b:9.3f}")
print(f"  {'median hit ADC':22s} {md_a:14.1f} {md_b:14.1f} {md_a / md_b:9.3f}")
print(f"  {'MEAN hit ADC':22s} {mean_a:14.1f} {mean_b:14.1f} {mean_a / mean_b:9.3f}")
print(f"  {'hits / acquisition':22s} {a['hits'].mean():14.1f} {b['hits'].mean():14.1f} "
      f"{a['hits'].mean() / b['hits'].mean():9.3f}")
print(f"  {'SCAs read out':22s} {a['scas'].mean():14.2f} {b['scas'].mean():14.2f} "
      f"{a['scas'].mean() / b['scas'].mean():9.3f}")
print(f"  {'acqWindow [ms]':22s} {a['acqwin']:14.1f} {b['acqwin']:14.1f} "
      f"{a['acqwin'] / b['acqwin']:9.3f}")
print(f"  {'delay between cyc [ms]':22s} {a['delay']:14.1f} {b['delay']:14.1f} "
      f"{a['delay'] / b['delay']:9.3f}")

print(f"\n  GAIN predicts the PEAK moves by the same factor as the mean ({mean_a / mean_b:.3f}).")
print(f"  PILE-UP predicts the peak does NOT move, and only the mean/tail grows.")
print(f"\n  The peak moved by {pk_a / pk_b:.3f}.  The mean moved by {mean_a / mean_b:.3f}.")

c = ROOT.TCanvas("c", "rate vs gain", 1700, 1200)
c.Divide(2, 2, 0.005, 0.005)
keep = []


def to_th1(h, name, colour):
    t = ROOT.TH1F(name, "", len(EDGES) - 1, EDGES[0], EDGES[-1])
    t.SetDirectory(0)
    for i, v in enumerate(h):
        t.SetBinContent(i + 1, v)
    t.SetLineColor(colour)
    t.SetLineWidth(3)
    t.SetFillColorAlpha(colour, 0.15)
    keep.append(t)
    return t


def spectrum(pad_id, xlo, xhi, logy, title):
    global keep
    pad = c.cd(pad_id)
    pad.SetMargin(0.13, 0.05, 0.13, 0.10)
    if logy:
        pad.SetLogy()
    ts = []
    for h, (run, th, colour), pk in ((ha, RUNS[0], pk_a), (hb_, RUNS[1], pk_b)):
        t = to_th1(h, f"s{pad_id}_{th}", colour)
        t.SetTitle(f"{title};adc_high - pedestal  [ADC];fraction of hits")
        t.GetXaxis().SetRangeUser(xlo, xhi)
        ts.append((t, run, th, colour, pk))
    ymax = max(t.GetMaximum() for t, _, _, _, _ in ts)
    leg = ROOT.TLegend(0.42, 0.60, 0.96, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.80)
    leg.SetTextSize(0.026)
    for i, (t, run, th, colour, pk) in enumerate(ts):
        if not logy:
            t.GetYaxis().SetRangeUser(0, ymax * 1.45)
        t.Draw("hist" if i == 0 else "hist same")
        leg.AddEntry(t, f"{run.split('_')[-1]} (th{th}): peak {pk:.0f} ADC", "f")
    for t, run, th, colour, pk in ts:
        ln = ROOT.TLine(pk, 0, pk, ymax * (1.15 if not logy else 1.0))
        ln.SetLineColor(colour)
        ln.SetLineWidth(3)
        ln.SetLineStyle(2)
        ln.Draw()
        keep.append(ln)
    leg.AddEntry(0, "", "")
    leg.AddEntry(0, f"peak ratio = {pk_a / pk_b:.3f}   mean ratio = {mean_a / mean_b:.3f}", "")
    leg.AddEntry(0, "GAIN: peak moves with the mean", "")
    leg.AddEntry(0, "PILE-UP: peak stays, only the mean moves", "")
    leg.Draw()
    keep.append(leg)


spectrum(1, 0, 350, True, "Per-hit ADC spectrum")
spectrum(2, 0, 90, False, "Zoom on the peak: DOES IT MOVE?")

# ---- 3. occupancy and rate ----------------------------------------------
pad = c.cd(3)
pad.SetMargin(0.13, 0.05, 0.16, 0.10)
fr = ROOT.TH1F("occ", "Is run12 actually a BUSIER beam? (pile-up needs it to be);;"
                      "run12 / run72", 4, 0, 4)
vals = [("hits per\nacquisition", a["hits"].mean() / b["hits"].mean()),
        ("SCAs read\nout", a["scas"].mean() / b["scas"].mean()),
        ("MEAN hit\nADC", mean_a / mean_b),
        ("PEAK hit\nADC", pk_a / pk_b)]
for i, (lab, _) in enumerate(vals):
    fr.GetXaxis().SetBinLabel(i + 1, lab.replace("\n", " "))
fr.GetXaxis().SetLabelSize(0.038)
fr.SetMinimum(0.8)
fr.SetMaximum(1.35)
fr.Draw()
keep.append(fr)
for i, (lab, v) in enumerate(vals):
    bx = ROOT.TBox(i + 0.15, 1.0, i + 0.85, v)
    bx.SetFillColorAlpha(ROOT.kRed + 1 if v > 1.02 else ROOT.kGray + 1, 0.55)
    bx.Draw()
    t = ROOT.TLatex(i + 0.5, v + 0.012, f"{v:.3f}")
    t.SetTextAlign(21)
    t.SetTextSize(0.040)
    t.Draw()
    keep += [bx, t]
one = ROOT.TLine(0, 1, 4, 1)
one.SetLineColor(ROOT.kBlack)
one.SetLineWidth(2)
one.SetLineStyle(2)
one.Draw()
keep.append(one)

# ---- 4. bin-by-bin ratio ------------------------------------------------
pad = c.cd(4)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
g = ROOT.TGraph()
for x, u, v in zip(CTR, ha, hb_):
    if v > 1e-6 and 4 < x < 300:
        g.SetPoint(g.GetN(), x, u / v)
g.SetMarkerStyle(20)
g.SetMarkerSize(0.8)
g.SetMarkerColor(ROOT.kBlack)
g.SetTitle("Ratio of the two spectra, bin by bin;adc_high - pedestal  [ADC];run12 / run72")
g.Draw("AP")
g.GetYaxis().SetRangeUser(0, 3.0)
z = ROOT.TLine(0, 1, 300, 1)
z.SetLineColor(ROOT.kGreen + 2)
z.SetLineWidth(3)
z.SetLineStyle(2)
z.Draw()
leg4 = ROOT.TLegend(0.32, 0.70, 0.96, 0.89)
leg4.SetBorderSize(0)
leg4.SetTextSize(0.026)
leg4.AddEntry(z, "1 = identical shape", "l")
leg4.AddEntry(0, "RISING with ADC  -> run12's hits are bigger (PILE-UP tail)", "")
leg4.AddEntry(0, "SHIFTED as a whole -> stretched spectrum (GAIN)", "")
leg4.Draw()
keep += [g, z, leg4]

out = os.path.join(OUT, "rate_vs_gain.png")
c.SaveAs(out)
print(f"\nsaved {out}")
