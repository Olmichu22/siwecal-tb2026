"""Beam spill structure as the DAQ sees it: is SpillSelect's cut really cutting
on spills, and what shape does the beam have?

PedestalMipCalibrator's SpillSelect (see plot_pedestal_multipeak_spill.py) calls
an acquisition "in spill" when its total accepted-hit occupancy clears a cut.
That is only legitimate if beam-on and beam-off acquisitions actually separate,
and if the ones above the cut group into blocks the way a spill does rather
than scattering at random. This draws the evidence.

Two figures:

  beam_spill_time.png -- the time domain
    1. occupancy vs acquisition index, with the cut drawn: the spill trains
    2. occupancy distribution: the bimodality the cut lives in, and how wide
       the valley is (any cut inside it gives the same split)
    3. spill and inter-spill block lengths, in acquisitions and in seconds
    4. averaged spill profile: occupancy vs position within a spill

  beam_spill_bcid.png -- inside one acquisition
    1. hit time (corrected BCID) distribution, in spill vs out
    2. hits vs SCA slot: which of the 15 SCA memory cells actually fill
    3. filled SCA slots per chip, in spill vs out -- this is WHY the
       out-of-spill pedestal sample is ~all SCA0: with no beam a chip is
       rarely triggered, so it never gets past its first memory cell
    4. mean hit time per SCA slot: how fast the memory fills, which is the
       trigger rate written in the chip's own clock

Acquisitions split across a raw chunk boundary are written twice, as two
COMPLEMENTARY halves with the same acqNumber; everything here groups by
acqNumber and SUMS, which splices them. Never dedup by position.

TIME AXIS: use `corrected_bcid`, never `bcid`. The raw BCID is a 12-bit counter
that WRAPS every 4096 ticks, and an acquisition runs far longer than that
(observed corrected_bcid up to ~38,600 on run_000060, i.e. ~9 wraps). Filled
with raw bcid every SCA slot returns a mean of ~2048 = 4096/2 -- a perfectly
flat, perfectly meaningless plot that looks like "no time structure". With
corrected_bcid the mean rises monotonically with slot number, as a switched-
capacitor array filling in time order must. `corrected_bcid` uses -999 as its
invalid sentinel; guard against it.

Usage:  RUN=<dir with chunks/> [OCC_CUT=200] [MAXCHUNKS=0] \
        python3 diagnostics/plot_beam_spill_structure.py [outdir]
"""
import array
import glob
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUT, exist_ok=True)

RUN = os.environ["RUN"]
OCC_CUT = float(os.environ.get("OCC_CUT", "200"))
MAXCHUNKS = int(os.environ.get("MAXCHUNKS", "0"))
NSL, NCHIP, NSCA = 15, 16, 15
# corrected_bcid is unwrapped, so it runs well past the 12-bit raw range.
NBCID = 65536

paths = sorted(glob.glob(os.path.join(RUN, "chunks", "chunk_*.root")))
if MAXCHUNKS:
    paths = paths[:MAXCHUNKS]
if not paths:
    raise SystemExit(f"no chunks under {RUN}")
print(f"{len(paths)} chunks from {RUN}")

bcid = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
cbcid = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
bad = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
nh = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
acqn = np.zeros(1, dtype=np.int32)
nslb = np.zeros(1, dtype=np.int32)
acqwin = np.zeros(1, dtype=np.float32)
delaycyc = np.zeros(1, dtype=np.float32)

occ_by_acq = {}                               # acqNumber -> summed occupancy
sca_by_acq = {}                               # acqNumber -> filled SCA slots
bcid_hits = {"in": np.zeros(NBCID), "out": np.zeros(NBCID)}
sca_hits = {"in": np.zeros(NSCA), "out": np.zeros(NSCA)}
bcid_sca = {"in": np.zeros((NSCA, NBCID)),    # hits vs (SCA slot, BCID)
            "out": np.zeros((NSCA, NBCID))}
slots_per_chip = {"in": [], "out": []}
period_ms = None

# Pass 1: occupancy per acquisition (needed before anything can be labelled
# in/out, since the label is a property of the whole acquisition).
raw = []
for path in paths:
    f = ROOT.TFile.Open(path)
    t = f.Get("siwecaldecoded")
    t.SetBranchStatus("*", 0)
    for b in ("acqNumber", "n_slboards", "bcid", "corrected_bcid", "badbcid", "nhits",
              "acqWindowMs", "delayBetweenCycleMs"):
        t.SetBranchStatus(b, 1)
    t.SetBranchAddress("acqNumber", acqn)
    t.SetBranchAddress("n_slboards", nslb)
    t.SetBranchAddress("bcid", bcid)
    t.SetBranchAddress("corrected_bcid", cbcid)
    t.SetBranchAddress("badbcid", bad)
    t.SetBranchAddress("nhits", nh)
    t.SetBranchAddress("acqWindowMs", acqwin)
    t.SetBranchAddress("delayBetweenCycleMs", delaycyc)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        nsl = int(nslb[0])
        good = (bad[:nsl] == 0) & (bcid[:nsl] >= 0)
        a = int(acqn[0])
        occ_by_acq[a] = occ_by_acq.get(a, 0) + int(nh[:nsl][good].sum())
        sca_by_acq[a] = sca_by_acq.get(a, 0) + int(good.sum())
        raw.append((path, i, a))
        if period_ms is None and float(acqwin[0]) > 0:
            period_ms = float(acqwin[0]) + max(float(delaycyc[0]), 0.0)
    f.Close()

acqs = np.array(sorted(occ_by_acq))
occ = np.array([occ_by_acq[a] for a in acqs], dtype=float)
inspill = occ > OCC_CUT
print(f"acquisitions (spliced by acqNumber): {len(acqs):,}   "
      f"in spill: {inspill.sum():,} ({100 * inspill.mean():.1f}%)")
print(f"acquisition period from the run config: {period_ms} ms"
      if period_ms else "acqWindowMs not available")

# Pass 2: the per-acquisition detail, now that each acquisition has a label.
label_of = {a: ("in" if occ_by_acq[a] > OCC_CUT else "out") for a in acqs}
by_file = {}
for path, i, a in raw:
    by_file.setdefault(path, []).append((i, a))
for path, items in by_file.items():
    f = ROOT.TFile.Open(path)
    t = f.Get("siwecaldecoded")
    t.SetBranchStatus("*", 0)
    for b in ("n_slboards", "bcid", "corrected_bcid", "badbcid", "nhits"):
        t.SetBranchStatus(b, 1)
    t.SetBranchAddress("n_slboards", nslb)
    t.SetBranchAddress("bcid", bcid)
    t.SetBranchAddress("corrected_bcid", cbcid)
    t.SetBranchAddress("badbcid", bad)
    t.SetBranchAddress("nhits", nh)
    for i, a in items:
        t.GetEntry(i)
        nsl = int(nslb[0])
        good = (bad[:nsl] == 0) & (bcid[:nsl] >= 0)
        lab = label_of[a]
        b_ = cbcid[:nsl][good]
        n_ = nh[:nsl][good]
        m = (b_ >= 0) & (b_ < NBCID)   # -999 = invalid sentinel
        np.add.at(bcid_hits[lab], b_[m], n_[m])
        sl_i, ch_i, sc_i = np.nonzero(good)
        np.add.at(sca_hits[lab], sc_i, nh[:nsl][good])
        bsel = cbcid[:nsl][good]
        ok = (bsel >= 0) & (bsel < NBCID)
        np.add.at(bcid_sca[lab], (sc_i[ok], bsel[ok]), nh[:nsl][good][ok])
        slots = good.reshape(nsl * NCHIP, NSCA).sum(axis=1)
        slots_per_chip[lab].append(slots[slots > 0])
    f.Close()

for lab in ("in", "out"):
    slots_per_chip[lab] = (np.concatenate(slots_per_chip[lab])
                           if slots_per_chip[lab] else np.zeros(0))

# --- spill / inter-spill block structure -------------------------------------
blocks = []
cur, start = inspill[0], 0
for i, v in enumerate(inspill[1:], start=1):
    if v != cur:
        blocks.append((cur, i - start))
        cur, start = v, i
blocks.append((cur, len(inspill) - start))
on = np.array([n for v, n in blocks if v], dtype=float)
off = np.array([n for v, n in blocks if not v], dtype=float)
print(f"spill blocks: {len(on)}  median {np.median(on):.0f} acquisitions"
      if len(on) else "no spill blocks")
if len(off):
    print(f"gaps        : {len(off)}  median {np.median(off):.0f} acquisitions")
if len(on) and len(off):
    # Duty cycle is the robust number: it is a ratio of acquisition COUNTS and
    # so survives whatever readout dead time acqWindowMs+delayBetweenCycleMs
    # leaves out. The seconds below are a LOWER bound on the real durations for
    # exactly that reason.
    print(f"beam duty cycle: {100 * inspill.mean():.1f}% of acquisitions")
if period_ms and len(on):
    print(f"=> spill >~{np.median(on) * period_ms / 1000:.1f} s, "
          f"gap >~{np.median(off) * period_ms / 1000:.1f} s, "
          f"period >~{(np.median(on) + np.median(off)) * period_ms / 1000:.1f} s "
          f"(config period {period_ms:.0f} ms/acquisition; readout dead time not "
          f"included, so these are lower bounds)")

# averaged spill profile: occupancy vs fractional position inside a block
NPROF = 20
prof = np.zeros(NPROF)
profn = np.zeros(NPROF)
idx = 0
for v, n in blocks:
    if v and n >= 5:
        seg = occ[idx:idx + n]
        pos = (np.arange(n) + 0.5) / n * NPROF
        np.add.at(prof, pos.astype(int).clip(0, NPROF - 1), seg)
        np.add.at(profn, pos.astype(int).clip(0, NPROF - 1), 1)
    idx += n
prof = np.divide(prof, np.maximum(profn, 1))

COL_IN, COL_OUT = ROOT.kOrange + 7, ROOT.kAzure + 2
keep = []

# =========================================================== figure 1: time ===
c = ROOT.TCanvas("c", "beam spill structure", 2000, 1600)
c.Divide(2, 2, 0.008, 0.012)

c.cd(1)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
ROOT.gPad.SetLogy()
# The WHOLE run, not a window. This used to draw only the first 1500
# acquisitions, which showed one and a half spill trains out of 29 and read as
# if that were all the beam there was -- with nothing in the title to say
# otherwise. Every acquisition is plotted now; markers only, since a connecting
# line across ~29k points is unreadable.
nshow = len(acqs)
g = ROOT.TGraph(nshow)
for i in range(nshow):
    g.SetPoint(i, i, max(occ[i], 0.5))
ntrain = int((on >= 10).sum()) if len(on) else 0
g.SetTitle(f"occupancy per acquisition -- whole run, {nshow:,} acquisitions, "
           f"{ntrain} spill trains;"
           "acquisition index;accepted hits in the acquisition")
g.SetMarkerColor(ROOT.kBlack)
g.SetMarkerStyle(20)
g.SetMarkerSize(0.25)
g.Draw("AP")
cut = ROOT.TLine(0, OCC_CUT, nshow, OCC_CUT)
cut.SetLineColor(ROOT.kRed + 1)
cut.SetLineStyle(2)
cut.SetLineWidth(3)
cut.Draw()
lat = ROOT.TLatex()
lat.SetNDC()
lat.SetTextSize(0.035)
lat.SetTextColor(ROOT.kRed + 1)
lat.DrawLatex(0.62, 0.855, f"SpillOccupancyCut = {OCC_CUT:.0f}")
keep += [g, cut, lat]

c.cd(2)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
ROOT.gPad.SetLogy()
ROOT.gPad.SetLogx()
edges = array.array("d", np.logspace(np.log10(0.5),
                                     np.log10(max(occ.max(), 10) * 1.5), 61))
hocc = ROOT.TH1F("hocc", "occupancy distribution: the valley the cut sits in;"
                         "accepted hits in the acquisition;acquisitions",
                 len(edges) - 1, edges)
for v in occ:
    hocc.Fill(max(v, 0.5))
hocc.SetLineColor(ROOT.kBlack)
hocc.SetLineWidth(3)
hocc.SetFillColorAlpha(ROOT.kGray + 1, 0.35)
hocc.Draw("hist")
cut2 = ROOT.TLine(OCC_CUT, 0.5, OCC_CUT, hocc.GetMaximum())
cut2.SetLineColor(ROOT.kRed + 1)
cut2.SetLineStyle(2)
cut2.SetLineWidth(3)
cut2.Draw()
keep += [hocc, cut2]

c.cd(3)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
ROOT.gPad.SetLogx()
top = max(on.max() if len(on) else 1, off.max() if len(off) else 1)
# log-x: the spill blocks (~tens) and the gaps (~hundreds) are an order of
# magnitude apart and a linear axis crushes the spills against zero.
bedges = array.array("d", np.logspace(0, np.log10(top * 1.3), 31))
hon = ROOT.TH1F("hon", "block lengths: spill vs gap;"
                       "consecutive acquisitions in the block;blocks",
                len(bedges) - 1, bedges)
hoff = ROOT.TH1F("hoff", "", len(bedges) - 1, bedges)
for v in on:
    hon.Fill(v)
for v in off:
    hoff.Fill(v)
for h, col in ((hon, COL_IN), (hoff, COL_OUT)):
    h.SetLineColor(col)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(col, 0.30)
hon.SetMaximum(1.6 * max(hon.GetMaximum(), hoff.GetMaximum()))
hon.Draw("hist")
hoff.Draw("hist same")
leg3 = ROOT.TLegend(0.45, 0.70, 0.93, 0.88)
leg3.SetBorderSize(0)
sfx = (f" ~{np.median(on) * period_ms / 1000:.1f} s" if period_ms and len(on) else "")
gfx = (f" ~{np.median(off) * period_ms / 1000:.1f} s" if period_ms and len(off) else "")
leg3.AddEntry(hon, f"spill (median {np.median(on):.0f} acq{sfx})", "f")
leg3.AddEntry(hoff, f"gap (median {np.median(off):.0f} acq{gfx})", "f")
leg3.Draw()
keep += [hon, hoff, leg3]

c.cd(4)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
hp = ROOT.TH1F("hp", "averaged spill profile;position within the spill  [%];"
                     "mean accepted hits per acquisition", NPROF, 0, 100)
for i, v in enumerate(prof, start=1):
    hp.SetBinContent(i, v)
hp.SetLineColor(COL_IN)
hp.SetLineWidth(3)
hp.SetFillColorAlpha(COL_IN, 0.30)
hp.SetMinimum(0)
hp.Draw("hist")
keep.append(hp)

p1 = os.path.join(OUT, "beam_spill_time.png")
c.SaveAs(p1)
print(f"saved {p1}")

# ================================================= figure 1b: zoom on 3 spills ===
# The whole-run panel above answers "are the trains regular?" but at ~29k points
# each train is a few pixels wide and its SHAPE is invisible. One pad, three
# consecutive trains, so the rise/plateau/fall (or absence of one) can be read.
# MIN_TRAIN separates a real spill from an isolated acquisition that happens to
# clear the cut. On run_000060, 3 of the 29 above-cut blocks are single
# acquisitions -- pile-up or a cosmic shower, not beam -- and picking one as a
# "spill" to zoom on would frame the window around nothing.
MIN_TRAIN = 10
starts, idx = [], 0
for v, n in blocks:
    if v and n >= MIN_TRAIN:
        starts.append((idx, n))
    idx += n
print(f"real spill trains (>={MIN_TRAIN} acquisitions): {len(starts)} of "
      f"{int(len(on))} above-cut blocks")
if starts:
    take = starts[:3]
    gap = int(np.median(off)) if len(off) else 200
    lo = max(take[0][0] - gap // 2, 0)
    hi = min(take[-1][0] + take[-1][1] + gap // 2, len(occ))
    cz = ROOT.TCanvas("cz", "spill zoom", 1800, 900)
    ROOT.gPad.SetMargin(0.10, 0.04, 0.12, 0.10)
    ROOT.gPad.SetLogy()
    sec = (period_ms or 0) / 1000.0
    xtitle = ("time since the first of these spills  [s]" if sec
              else "acquisition index")
    gz = ROOT.TGraph(hi - lo)
    for i in range(lo, hi):
        x = (i - lo) * sec if sec else i
        gz.SetPoint(i - lo, x, max(occ[i], 0.5))
    gz.SetTitle(f"{len(take)} consecutive spills, "
                f"{hi - lo:,} acquisitions;{xtitle};accepted hits in the acquisition")
    gz.SetLineColor(ROOT.kGray + 2)
    gz.SetMarkerColor(ROOT.kBlack)
    gz.SetMarkerStyle(20)
    gz.SetMarkerSize(0.6)
    gz.Draw("APL")
    xmax = (hi - lo - 1) * sec if sec else hi - 1
    cutz = ROOT.TLine(0 if sec else lo, OCC_CUT, xmax, OCC_CUT)
    cutz.SetLineColor(ROOT.kRed + 1)
    cutz.SetLineStyle(2)
    cutz.SetLineWidth(3)
    cutz.Draw()
    latz = ROOT.TLatex()
    latz.SetNDC()
    latz.SetTextSize(0.035)
    latz.SetTextColor(ROOT.kRed + 1)
    latz.DrawLatex(0.68, 0.855, f"SpillOccupancyCut = {OCC_CUT:.0f}")
    pz = os.path.join(OUT, "beam_spill_zoom.png")
    cz.SaveAs(pz)
    print(f"saved {pz}  (spills of {[n for _, n in take]} acquisitions)")

# ========================================================== figure 2: bcid ===
c2 = ROOT.TCanvas("c2", "inside the acquisition", 2000, 1600)
c2.Divide(2, 2, 0.008, 0.012)
keep2 = []

c2.cd(1)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
ROOT.gPad.SetLogy()
last = max(int(np.nonzero(bcid_hits["in"])[0].max() if bcid_hits["in"].any() else 100),
           int(np.nonzero(bcid_hits["out"])[0].max() if bcid_hits["out"].any() else 100))
nb = 200
hb = {}
for lab, col in (("in", COL_IN), ("out", COL_OUT)):
    h = ROOT.TH1F(f"hb_{lab}", "hit time within the acquisition;"
                               "corrected BCID (unwrapped);hits (area-normalised)", nb, 0, last * 1.02)
    for b, v in enumerate(bcid_hits[lab][:last + 1]):
        if v:
            h.Fill(b, v)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(col)
    h.SetLineWidth(3)
    hb[lab] = h
    keep2.append(h)
hb["in"].SetMaximum(5 * max(hb["in"].GetMaximum(), hb["out"].GetMaximum()))
hb["in"].Draw("hist")
hb["out"].Draw("hist same")
legb = ROOT.TLegend(0.55, 0.74, 0.93, 0.88)
legb.SetBorderSize(0)
legb.AddEntry(hb["in"], "in spill", "l")
legb.AddEntry(hb["out"], "out of spill", "l")
legb.Draw()
keep2.append(legb)

c2.cd(2)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
ROOT.gPad.SetLogy()
hsc = {}
for lab, col in (("in", COL_IN), ("out", COL_OUT)):
    h = ROOT.TH1F(f"hsc_{lab}", "hits per SCA slot;SCA;hits (area-normalised)",
                  NSCA, 0, NSCA)
    for s in range(NSCA):
        h.SetBinContent(s + 1, sca_hits[lab][s])
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(col)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(col, 0.25)
    hsc[lab] = h
    keep2.append(h)
hsc["in"].SetMaximum(8 * max(hsc["in"].GetMaximum(), hsc["out"].GetMaximum()))
hsc["in"].Draw("hist")
hsc["out"].Draw("hist same")
legs = ROOT.TLegend(0.55, 0.74, 0.93, 0.88)
legs.SetBorderSize(0)
legs.AddEntry(hsc["in"], "in spill", "f")
legs.AddEntry(hsc["out"], "out of spill", "f")
legs.Draw()
keep2.append(legs)

c2.cd(3)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
ROOT.gPad.SetLogy()
hsl = {}
for lab, col in (("in", COL_IN), ("out", COL_OUT)):
    h = ROOT.TH1F(f"hsl_{lab}", "filled SCA slots per triggered chip "
                                "(why the quiet sample is all SCA0);"
                                "SCA slots filled;chip-acquisitions (normalised)",
                  NSCA, 0.5, NSCA + 0.5)
    for v in slots_per_chip[lab]:
        h.Fill(v)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(col)
    h.SetLineWidth(3)
    h.SetFillColorAlpha(col, 0.25)
    hsl[lab] = h
    keep2.append(h)
hsl["in"].SetMaximum(8 * max(hsl["in"].GetMaximum(), hsl["out"].GetMaximum()))
hsl["in"].Draw("hist")
hsl["out"].Draw("hist same")
legl = ROOT.TLegend(0.5, 0.74, 0.93, 0.88)
legl.SetBorderSize(0)
med_in = np.median(slots_per_chip["in"]) if len(slots_per_chip["in"]) else 0
med_out = np.median(slots_per_chip["out"]) if len(slots_per_chip["out"]) else 0
legl.AddEntry(hsl["in"], f"in spill (median {med_in:.0f})", "f")
legl.AddEntry(hsl["out"], f"out of spill (median {med_out:.0f})", "f")
legl.Draw()
keep2.append(legl)

c2.cd(4)
ROOT.gPad.SetMargin(0.12, 0.12, 0.12, 0.10)
ROOT.gPad.SetLogz()
# How fast the chip's 15 memory cells fill. An SCA array fills in time order,
# so the mean BCID of slot n rises with n, and the SLOPE is the inverse trigger
# rate: a steep line means the chip waits a long time between triggers, a flat
# one means it fills its whole memory early in the window. Comparing in/out of
# spill turns that into a beam-intensity measurement in the chip's own clock.
mb = {}
for lab, col in (("in", COL_IN), ("out", COL_OUT)):
    g = ROOT.TGraph()
    n = 0
    for s in range(NSCA):
        w = bcid_sca[lab][s]
        tot = w.sum()
        if tot <= 0:
            continue
        g.SetPoint(n, s, float((np.arange(NBCID) * w).sum() / tot))
        n += 1
    g.SetLineColor(col)
    g.SetMarkerColor(col)
    g.SetLineWidth(3)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.4)
    mb[lab] = g
    keep2.append(g)
fr = ROOT.TH1F("fr", "how fast the SCA memory fills "
                     "(slope = inverse trigger rate);SCA slot;mean corrected BCID",
               NSCA, -0.5, NSCA - 0.5)
allm = [mb[l].GetY()[i] for l in mb for i in range(mb[l].GetN())]
fr.SetMinimum(0)
fr.SetMaximum(max(allm) * 1.45 if allm else 1)
fr.Draw("axis")
for lab in ("in", "out"):
    mb[lab].Draw("PL same")
legm = ROOT.TLegend(0.16, 0.74, 0.60, 0.88)
legm.SetBorderSize(0)
legm.AddEntry(mb["in"], "in spill", "lp")
legm.AddEntry(mb["out"], "out of spill", "lp")
legm.Draw()
keep2 += [fr, legm]

p2 = os.path.join(OUT, "beam_spill_bcid.png")
c2.SaveAs(p2)
print(f"saved {p2}")
