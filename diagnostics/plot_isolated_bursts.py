"""The isolated high-occupancy acquisitions BETWEEN spill trains -- what are they?

plot_beam_spill_structure.py finds that a few acquisitions clear
SpillOccupancyCut while sitting alone, not inside a spill train (3 of 29
above-cut blocks on run_000060). They were initially guessed to be cosmic
showers. They are not, and they are not even one phenomenon: this draws each
one along the two axes that separate the possibilities.

  WHERE (hits per slab)      a through-going particle or a shower lights many
                             slabs; a detector pathology lights one.
  WHEN (hits per corrected   one event deposits at a single corrected_bcid;
        BCID)                a retriggering chip smears over the window.

On run_000060 the three fall into three different classes, which is the point
of the figure:

  acq 2631   15 slabs, top chip 16% of the hits, 27 distinct BCIDs
             -> detector-wide activity at several times. Beam-like: halo or
                leakage outside the train, not a single event.
  acq 16038  all 236 hits in slab 12, spread over many of its chips, the top
             chip filling ONE SCA slot, everything at a single BCID (2998)
             -> a whole slab firing coincidentally in one instant. Not a
                retrigger (that fills many SCAs of one chip) and not a track.
                Note slab 12 is the COB slab, which is worth keeping in mind
                before generalising this to the other layers.
  acq 26053  all in slab 1, 3 chips, top chip filling 7 SCA slots spanning
             corrected_bcid 3207 -> 13429
             -> a chip retriggering across the acquisition, the classic
                SKIROC pathology.

CAVEAT, and it is a large one: this is THREE acquisitions in ONE run. Nothing
here is a rate, a classification, or a population. It is three case studies
that refute "they are all cosmic showers" and show the categories are mixed.
Extending to the other eight th220 runs would give ~27 and is the obvious next
step if these matter.

Usage:  RUN=<dir with chunks/> [OCC_CUT=200] [MIN_TRAIN=10] \
        python3 diagnostics/plot_isolated_bursts.py [outdir]
"""
import glob
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUT, exist_ok=True)

RUN = os.environ["RUN"]
OCC_CUT = float(os.environ.get("OCC_CUT", "200"))
MIN_TRAIN = int(os.environ.get("MIN_TRAIN", "10"))
NSL, NCHIP, NSCA, NBCID = 15, 16, 15, 65536

bcid = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
cbcid = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
bad = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
nh = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
acqn = np.zeros(1, dtype=np.int32)
nslb = np.zeros(1, dtype=np.int32)


def bind(t):
    t.SetBranchStatus("*", 0)
    for b in ("acqNumber", "n_slboards", "bcid", "corrected_bcid", "badbcid", "nhits"):
        t.SetBranchStatus(b, 1)
    t.SetBranchAddress("acqNumber", acqn)
    t.SetBranchAddress("n_slboards", nslb)
    t.SetBranchAddress("bcid", bcid)
    t.SetBranchAddress("corrected_bcid", cbcid)
    t.SetBranchAddress("badbcid", bad)
    t.SetBranchAddress("nhits", nh)


paths = sorted(glob.glob(os.path.join(RUN, "chunks", "chunk_*.root")))
if not paths:
    raise SystemExit(f"no chunks under {RUN}")

# Pass 1: occupancy per acquisition (spliced by acqNumber -- a chunk-boundary
# acquisition is written twice, as two complementary halves).
occ_by_acq, where = {}, {}
for path in paths:
    f = ROOT.TFile.Open(path)
    t = f.Get("siwecaldecoded")
    bind(t)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        n = int(nslb[0])
        g = (bad[:n] == 0) & (bcid[:n] >= 0)
        a = int(acqn[0])
        occ_by_acq[a] = occ_by_acq.get(a, 0) + int(nh[:n][g].sum())
        where.setdefault(a, []).append((path, i))
    f.Close()

acqs = np.array(sorted(occ_by_acq))
occ = np.array([occ_by_acq[a] for a in acqs], dtype=float)
ins = occ > OCC_CUT

blocks, cur, st = [], ins[0], 0
for i, v in enumerate(ins[1:], 1):
    if v != cur:
        blocks.append((cur, st, i - st))
        cur, st = v, i
blocks.append((cur, st, len(ins) - st))
isolated = [st_ for v, st_, n in blocks if v and n < MIN_TRAIN]
ntrain = sum(1 for v, _, n in blocks if v and n >= MIN_TRAIN)
print(f"{len(acqs):,} acquisitions | {ntrain} spill trains | "
      f"{len(isolated)} isolated above-cut acquisitions")
if not isolated:
    raise SystemExit("no isolated above-cut acquisitions in this run")


def detail(pos):
    """Per-slab hits and the hits-vs-corrected_bcid profile of one acquisition."""
    a = int(acqs[pos])
    per_slab = np.zeros(NSL)
    per_chip = np.zeros((NSL, NCHIP))
    prof = {}
    for path, i in where[a]:
        f = ROOT.TFile.Open(path)
        t = f.Get("siwecaldecoded")
        bind(t)
        t.GetEntry(i)
        n = int(nslb[0])
        g = (bad[:n] == 0) & (bcid[:n] >= 0)
        sl_i, ch_i, _ = np.nonzero(g)
        np.add.at(per_slab, sl_i, nh[:n][g])
        np.add.at(per_chip, (sl_i, ch_i), nh[:n][g])
        for b_, h_ in zip(cbcid[:n][g], nh[:n][g]):
            if 0 <= b_ < NBCID:
                prof[int(b_)] = prof.get(int(b_), 0) + int(h_)
        f.Close()
    return a, per_slab, per_chip, prof


shown = isolated[:3]
ncol = len(shown)
c = ROOT.TCanvas("c", "isolated bursts", 700 * ncol, 1100)
c.Divide(ncol, 2, 0.008, 0.012)
keep = []
COL = ROOT.kViolet + 1

for j, pos in enumerate(shown):
    a, per_slab, per_chip, prof = detail(pos)
    tot = per_slab.sum()
    nslab_hit = int((per_slab > 0).sum())
    top = per_chip.max() / max(tot, 1)
    bmax = max(prof.values()) if prof else 0
    frac_top_bcid = bmax / max(tot, 1)

    # --- top row: WHERE
    c.cd(j + 1)
    ROOT.gPad.SetMargin(0.13, 0.05, 0.13, 0.12)
    hs = ROOT.TH1F(f"hs{j}", f"acq {a}  --  {int(tot)} hits, {nslab_hit}/15 slabs hit;"
                             f"slab;hits", NSL, 0, NSL)
    for s in range(NSL):
        hs.SetBinContent(s + 1, per_slab[s])
    hs.SetLineColor(COL)
    hs.SetLineWidth(3)
    hs.SetFillColorAlpha(COL, 0.30)
    hs.SetMinimum(0)
    hs.Draw("hist")
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.042)
    lat.DrawLatex(0.17, 0.82, f"busiest chip: {100 * top:.0f}% of hits")
    keep += [hs, lat]

    # --- bottom row: WHEN
    c.cd(ncol + j + 1)
    ROOT.gPad.SetMargin(0.13, 0.05, 0.13, 0.12)
    bs = sorted(prof)
    span = (max(bs) - min(bs)) if len(bs) > 1 else 1
    lo = max(min(bs) - 0.05 * span, 0)
    hi = max(bs) + 0.05 * span + 1
    hb = ROOT.TH1F(f"hb{j}", f"acq {a}  --  {len(prof)} distinct BCIDs;"
                             f"corrected BCID (unwrapped);hits", 120, lo, hi)
    for b_, v in prof.items():
        hb.Fill(b_, v)
    hb.SetLineColor(COL)
    hb.SetLineWidth(3)
    hb.SetFillColorAlpha(COL, 0.30)
    hb.SetMinimum(0)
    hb.Draw("hist")
    lat2 = ROOT.TLatex()
    lat2.SetNDC()
    lat2.SetTextSize(0.042)
    lat2.DrawLatex(0.17, 0.82, f"biggest BCID: {100 * frac_top_bcid:.0f}% of hits")
    keep += [hb, lat2]

    print(f"  acq {a}: {int(tot)} hits, {nslab_hit}/15 slabs, "
          f"busiest chip {100 * top:.0f}%, {len(prof)} BCIDs, "
          f"biggest BCID {100 * frac_top_bcid:.0f}%")

path_out = os.path.join(OUT, "isolated_bursts.png")
c.SaveAs(path_out)
print(f"saved {path_out}")
