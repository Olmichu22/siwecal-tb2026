"""The isolated out-of-spill bursts as a POPULATION, over every th220 run.

plot_isolated_bursts.py looked at three acquisitions from run_000060 and showed
they are three different things -- enough to refute "they are all cosmic
showers", not enough to say what the mixture is or how often it happens. This
sweeps all nine th220 runs -- and the answer is that they are NOT a population
with a rate. The number of spill trains is stable in the occupancy cut and is a
property of the beam; the number of isolated acquisitions is not, falling by
more than an order of magnitude across cuts that are all defensible (run_000060:
169 at cut 25, 29 at 100, 3 at 200). Panel 2 is that measurement, and it is the
point of the figure. The three acquisitions characterised by hand remain real
objects; what does not survive is counting them or quoting a probability.

An acquisition is ISOLATED if it clears the run's spill cut while sitting in an
above-cut block shorter than MIN_TRAIN acquisitions, i.e. it is not part of a
spill train.

THE CUT IS DERIVED PER RUN, anchored to that run's own beam level; see
spill_cut() for the rule and for the two wrong ones that came before it. A run
whose above-cut acquisitions are mostly NOT contiguous has no spill structure at
all and is excluded rather than forced through the machinery.

Each isolated acquisition is then reduced to three numbers, the same axes the
three case studies separated on:

    nslab      how many slabs recorded a hit          (localised vs detector-wide)
    topchip    fraction of hits on the busiest chip   (one chip vs spread out)
    nbcid      distinct corrected BCIDs               (one instant vs smeared)
    topbcid    fraction of hits in the biggest BCID

and sorted by a FIXED RULE, stated here so it can be argued with rather than
reverse-engineered from the plot:

    retrigger    topchip >= 0.60 and nbcid >= 3      one chip firing repeatedly
    flash        nslab <= 3 and topbcid >= 0.60      a few slabs, one instant
    beam-like    nslab >= 8                          detector-wide activity
    other        anything else

The rule is a summary of the three cases already characterised by hand; it is
not derived from this sample, and the "other" bin is left in the plot precisely
so that a bad rule shows up instead of hiding.

BCID must be the CORRECTED one. The raw bcid is a 12-bit counter that wraps
every 4096 within an acquisition, so counting distinct raw BCIDs silently
saturates and every acquisition comes out looking equally smeared.

Usage:
  RUNS=<comma-separated run dirs> [MIN_TRAIN=10] [FORCE_CUT=<hits>] \
      python3 diagnostics/plot_isolated_bursts_population.py [outdir]

FORCE_CUT overrides the per-run cut with a fixed one, which is how to reproduce
the run_000060 numbers of section 9 (FORCE_CUT=200).
"""
import glob
import os
import re
import sys
from collections import Counter

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUTDIR, exist_ok=True)

MIN_TRAIN = int(os.environ.get("MIN_TRAIN", "10"))
FORCE_CUT = float(os.environ["FORCE_CUT"]) if os.environ.get("FORCE_CUT") else None
CACHE = os.environ.get("CACHE", "")
if CACHE and not os.path.isdir(CACHE):
    os.makedirs(CACHE, exist_ok=True)
NSL, NCHIP, NSCA, NBCID = 15, 16, 15, 65536

BEAM_PCT = 99.0    # percentile above which the beam level is measured
CUT_FRAC = 0.2     # the cut, as a fraction of that beam level


def beam_level(o):
    """Occupancy of a well-inside-the-spill acquisition: the median above p99."""
    hi = o > np.percentile(o, BEAM_PCT)
    return float(np.median(o[hi])) if hi.any() else 0.0


def spill_cut(o):
    """The run's own occupancy cut: a fixed fraction of its own beam level.

    NOT a valley finder, and that is the result of two failed attempts rather
    than a shortcut. These distributions are not cleanly bimodal -- between the
    quiet mode at a few hits and the beam mode at ~10^3 there is a filled
    continuum, so "the valley" is not a well-defined place and any argmin lands
    somewhere arbitrary in it. The first version returned cuts of 0-4 hits for
    every run (92,067 "isolated bursts", a tenth of all acquisitions, because at
    low occupancy the small integers leave exactly-empty log bins next to the
    quiet peak); the second, smoothed and floored, returned cuts scattered from
    25 to 327 hits with no physical reason for the spread.

    Anchoring to the run's OWN beam level instead is stable and reproduces what
    section 7 chose by eye on run_000060: beam level 1,134 hits, cut 227, against
    the 200 that was picked by looking at the histogram.
    """
    b = beam_level(o)
    return (b * CUT_FRAC if b > 0 else None), dict(beam=b)


runs = [r for r in os.environ["RUNS"].split(",") if r]
runs = [r for r in runs if os.path.isdir(os.path.join(r, "chunks"))]
if not runs:
    raise SystemExit("no run directories with chunks/ in RUNS")

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


CUT_SCAN = [25, 50, 100, 200, 400, 800]
MIN_TRAIN_FRACTION = 0.5   # of above-cut acquisitions, to call a run "has beam"


def train_fraction(o, cut):
    """Fraction of above-cut acquisitions that sit inside a spill train."""
    a = o > cut
    if not a.any():
        return 0.0
    blocks, cur, st = [], a[0], 0
    for i, v in enumerate(a[1:], 1):
        if v != cur:
            blocks.append((cur, st, i - st))
            cur, st = v, i
    blocks.append((cur, st, len(a) - st))
    above = sum(n for v, _, n in blocks if v)
    inside = sum(n for v, _, n in blocks if v and n >= MIN_TRAIN)
    return inside / max(above, 1)


def block_counts(o, cut):
    """(number of spill trains, positions of the isolated above-cut blocks)."""
    a = o > cut
    if not a.any():
        return 0, []
    blocks, cur, st = [], a[0], 0
    for i, v in enumerate(a[1:], 1):
        if v != cur:
            blocks.append((cur, st, i - st))
            cur, st = v, i
    blocks.append((cur, st, len(a) - st))
    return (sum(1 for v, _, n in blocks if v and n >= MIN_TRAIN),
            [s for v, s, n in blocks if v and n < MIN_TRAIN])


def classify(nslab, topchip, nbc, topbcid):
    if topchip >= 0.60 and nbc >= 3:
        return "retrigger"
    if nslab <= 3 and topbcid >= 0.60:
        return "flash"
    if nslab >= 8:
        return "beam-like"
    return "other"


CLASSES = ["beam-like", "flash", "retrigger", "other"]
COLS = {"beam-like": ROOT.kAzure + 2, "flash": ROOT.kOrange + 7,
        "retrigger": ROOT.kRed + 1, "other": ROOT.kGray + 2}

per_run = []
records = []           # one dict per isolated acquisition
occ_diag = {}          # run -> (occupancies, chosen cut, diagnostic)
no_beam = []           # runs with no spill train at their own cut
scan = {}              # run -> [(cut, ntrain, niso), ...] over CUT_SCAN
for rdir in runs:
    short = os.path.basename(rdir.rstrip("/")).replace("TB2026CERN_", "").replace("eudaq_", "")
    paths = sorted(glob.glob(os.path.join(rdir, "chunks", "chunk_*.root")))
    if not paths:
        continue

    # Pass 1: occupancy per acquisition. Acquisitions split across a chunk
    # boundary are written twice as complementary halves sharing an acqNumber,
    # so occupancy must be SUMMED per acqNumber, never taken per entry.
    occ, where = {}, {}
    for p in paths:
        f = ROOT.TFile.Open(p)
        t = f.Get("siwecaldecoded")
        bind(t)
        for i in range(t.GetEntries()):
            t.GetEntry(i)
            n = int(nslb[0])
            g = (bad[:n] == 0) & (bcid[:n] >= 0)
            a = int(acqn[0])
            occ[a] = occ.get(a, 0) + int(nh[:n][g].sum())
            where.setdefault(a, []).append((p, i))
        f.Close()

    acqs = np.array(sorted(occ))
    o = np.array([occ[a] for a in acqs], dtype=float)
    # Pass 1 is the expensive part and the cut is the part most likely to need
    # another look, so keep the occupancies. Re-deriving cuts and scans against
    # these costs seconds instead of another sweep over 1,377 chunks -- which is
    # how the cut rule below was settled, after two wrong ones.
    if CACHE:
        np.savez_compressed(os.path.join(CACHE, f"occ_{short}.npz"), acqs=acqs, occ=o)
    cut, diag = ((FORCE_CUT, dict(beam=beam_level(o))) if FORCE_CUT
                 else spill_cut(o))
    occ_diag[short] = (o, cut, diag)
    if cut is None:
        no_beam.append(short)
        print(f"{short:<12} {len(paths):>4} chunks  {len(acqs):>7,} acquisitions  "
              f"NO BEAM LEVEL -- skipped "
              f"(occupancy median {np.median(o):.0f}, max {o.max():.0f})")
        continue

    ntrain, iso_pos = block_counts(o, cut)
    # IS THERE A SPILL AT ALL? A spill is contiguous, so in a run with beam
    # almost every above-cut acquisition sits inside a train. Measured at each
    # run's own cut, the separation is not marginal: 98-100% for six of the nine
    # runs and 75.6% for run_000062, against 5.6% for run_000142 and 2.8% for
    # run_000143. Those two have no beam -- their "beam level" is just the tail
    # of a beamless run (122 and 109 hits, against 1,066-1,854 elsewhere) -- and
    # letting them through contributed 15,284 phantom bursts out of 15,624.
    #
    # This test is used instead of a floor on the occupancy because a floor in
    # hits is exactly the kind of absolute number that failed to transport
    # between runs twice already in this study.
    duty = train_fraction(o, cut)
    if duty < MIN_TRAIN_FRACTION:
        no_beam.append(short)
        print(f"{short:<12} {len(paths):>4} chunks  {len(acqs):>7,} acquisitions  "
              f"cut {cut:>6.0f}  NO SPILL STRUCTURE -- skipped "
              f"(only {100 * duty:.1f}% of above-cut acquisitions are in trains, "
              f"beam level {diag['beam']:.0f} hits)")
        continue

    # How much of this depends on where the cut went? The scan is the honest
    # answer and it belongs in the output, not in a footnote: the train count
    # turns out to be stable and the isolated count does not.
    scan[short] = [(c, ) + tuple(len(x) if isinstance(x, list) else x
                                 for x in block_counts(o, c))
                   for c in CUT_SCAN]

    # Pass 2: only the isolated ones, read in full.
    for pos in iso_pos:
        a = int(acqs[pos])
        per_slab = np.zeros(NSL)
        per_chip = np.zeros((NSL, NCHIP))
        prof = {}
        for p, i in where[a]:
            f = ROOT.TFile.Open(p)
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
        tot = max(per_slab.sum(), 1)
        rec = dict(run=short, acq=a, hits=int(per_slab.sum()),
                   nslab=int((per_slab > 0).sum()),
                   topchip=float(per_chip.max() / tot),
                   nbcid=len(prof),
                   topbcid=float(max(prof.values()) / tot) if prof else 0.0)
        rec["cls"] = classify(rec["nslab"], rec["topchip"], rec["nbcid"], rec["topbcid"])
        records.append(rec)

    per_run.append(dict(run=short, nacq=len(acqs), ntrain=ntrain, niso=len(iso_pos),
                        nchunk=len(paths), cut=cut))
    print(f"{short:<12} {len(paths):>4} chunks  {len(acqs):>7,} acquisitions  "
          f"cut {cut:>6.0f} hits  {ntrain:>4} trains  {len(iso_pos):>4} isolated")

if no_beam:
    print(f"\nrun(s) with no spill structure, excluded from every number below: "
          f"{', '.join(no_beam)}")
if not records:
    raise SystemExit("no isolated above-cut acquisitions in any run with beam")

tot_acq = sum(r["nacq"] for r in per_run)
tot_iso = len(records)
tot_train = sum(r["ntrain"] for r in per_run)
counts = Counter(r["cls"] for r in records)
print(f"\n{tot_acq:,} acquisitions over {len(per_run)} runs, {tot_train} spill trains, "
      f"{tot_iso} isolated above-cut acquisitions at each run's own cut.")
print("DO NOT quote a rate from that total: see the cut scan below -- the train "
      "count is stable in the cut and the isolated count is not.")
print("\ncut scan (trains / isolated):")
hdr = "".join(f"{c:>12}" for c in CUT_SCAN)
print(f"   {'run':<10}{hdr}")
for short, rows in sorted(scan.items()):
    print(f"   {short:<10}" + "".join(f"{nt:>6}/{ni:<6}" for _c, nt, ni in rows))
print("classification (fixed rule, see the docstring):")
for cl in CLASSES:
    k = counts.get(cl, 0)
    print(f"   {cl:<11} {k:>4}  ({100 * k / tot_iso:5.1f}%)")

# ------------------------------------------------------------------ figure
c = ROOT.TCanvas("c", "isolated bursts population", 1800, 1300)
c.Divide(2, 2, 0.010, 0.014)
keep = []

# --- 1. why the cut has to be per run: every run's occupancy distribution with
# the cut it chose. A single fixed cut cannot sit in all of these valleys.
c.cd(1)
ROOT.gPad.SetMargin(0.13, 0.05, 0.13, 0.10)
ROOT.gPad.SetLogy()
PAL = [ROOT.kAzure + 2, ROOT.kOrange + 7, ROOT.kRed + 1, ROOT.kGreen + 2,
       ROOT.kViolet + 1, ROOT.kCyan + 2, ROOT.kMagenta + 1, ROOT.kYellow + 3,
       ROOT.kGray + 2]
occ_h, occ_lines = [], []
first = True
lo1 = ROOT.TLegend(0.60, 0.46, 0.96, 0.89)
lo1.SetBorderSize(0)
lo1.SetFillStyle(0)
lo1.SetTextSize(0.028)
for i, (short, (o, cut, _d)) in enumerate(sorted(occ_diag.items())):
    h = ROOT.TH1F(f"ho{i}", "Acquisition occupancy, run by run;"
                            "log_{10}(total accepted hits + 1);"
                            "fraction of acquisitions", 60, 0, 4.2)
    for v in np.log10(np.maximum(o, 0) + 1.0):
        h.Fill(v)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    col = PAL[i % len(PAL)]
    h.SetLineColor(col)
    h.SetLineWidth(2)
    h.SetMinimum(1e-5)
    h.SetMaximum(8.0)
    h.Draw("hist" if first else "hist same")
    first = False
    # Chopping a fixed "run_0000" prefix silently does nothing for runs 142,
    # 143 and 254, which are written with one zero fewer -- they came out on the
    # legend as the full "run_000142" while the others showed a bare number.
    m_no = re.search(r"run_0*(\d+)", short)
    tag = m_no.group(1) if m_no else short
    if short in no_beam:
        lo1.AddEntry(h, f"{tag}: no spill structure", "l")
    elif cut is not None:
        lo1.AddEntry(h, f"{tag}: cut {cut:.0f}", "l")
        ln = ROOT.TLine(np.log10(cut + 1), 1e-5, np.log10(cut + 1), 0.3)
        ln.SetLineColor(col)
        ln.SetLineStyle(2)
        ln.Draw()
        occ_lines.append(ln)
    occ_h.append(h)
lo1.Draw()
keep += occ_h + occ_lines + [lo1]

# --- 2. THE panel: how much of the answer is the cut's doing. The train count
# and the isolated count are drawn on the same axes precisely so the contrast
# is unavoidable -- one is a property of the beam, the other is not.
c.cd(2)
ROOT.gPad.SetMargin(0.13, 0.05, 0.13, 0.10)
ROOT.gPad.SetLogx()
ROOT.gPad.SetLogy()
fr2 = ROOT.TH1F("fr2", "Is it a measurement or a choice of cut?;"
                       "occupancy cut (hits);count in the run", 100, 20, 1000)
fr2.SetMinimum(0.5)
fr2.SetMaximum(3e4)
fr2.Draw()
scan_g = []
ls2 = ROOT.TLegend(0.15, 0.66, 0.55, 0.88)
ls2.SetBorderSize(0)
ls2.SetFillStyle(0)
ls2.SetTextSize(0.030)
for i, (short, rows) in enumerate(sorted(scan.items())):
    col = PAL[i % len(PAL)]
    for what, style, width in (("trains", 1, 3), ("isolated", 2, 2)):
        g = ROOT.TGraph()
        j = 0
        for cutv, nt, ni in rows:
            y = nt if what == "trains" else ni
            if y > 0:
                g.SetPoint(j, cutv, y)
                j += 1
        if j == 0:
            continue
        g.SetLineColor(col)
        g.SetLineStyle(style)
        g.SetLineWidth(width)
        g.Draw("lsame")
        scan_g.append(g)
lat2 = ROOT.TLatex()
lat2.SetNDC()
lat2.SetTextSize(0.031)
lat2.DrawLatex(0.16, 0.83, "solid: spill trains -- flat, so it is a property of the beam")
lat2.DrawLatex(0.16, 0.78, "dashed: isolated acquisitions -- falls by 20x or more")
lat2.DrawLatex(0.16, 0.73, "no rate can be quoted for the isolated ones")
keep += scan_g + [fr2, lat2]

# --- 3. the two discriminating axes, one point per burst
c.cd(3)
ROOT.gPad.SetMargin(0.13, 0.05, 0.13, 0.10)
frame = ROOT.TH2F("frame", "Every isolated out-of-spill acquisition;"
                           "slabs hit;fraction of hits on the busiest chip",
                  16, 0, 16, 20, 0, 1.05)
frame.Draw()
graphs = {}
for cl in CLASSES:
    pts = [r for r in records if r["cls"] == cl]
    if not pts:
        continue
    g = ROOT.TGraph(len(pts))
    for i, r in enumerate(pts):
        # jitter x a little so overlapping integer slab counts stay countable
        g.SetPoint(i, r["nslab"] + np.random.uniform(-0.22, 0.22), r["topchip"])
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.1)
    g.SetMarkerColor(COLS[cl])
    g.Draw("psame")
    graphs[cl] = g
l1 = ROOT.TLegend(0.45, 0.60, 0.94, 0.86)
l1.SetBorderSize(0); l1.SetFillStyle(0)
for cl, g in graphs.items():
    l1.AddEntry(g, f"{cl} ({counts[cl]})", "p")
l1.Draw()
keep += [frame, l1] + list(graphs.values())

# --- 4. the mixture and the rate
c.cd(4)
ROOT.gPad.SetMargin(0.13, 0.05, 0.13, 0.10)
hc = ROOT.TH1F("hc", "What the isolated acquisitions are;;isolated acquisitions",
               len(CLASSES), 0, len(CLASSES))
for i, cl in enumerate(CLASSES):
    hc.SetBinContent(i + 1, counts.get(cl, 0))
    hc.GetXaxis().SetBinLabel(i + 1, cl)
hc.GetXaxis().SetLabelSize(0.055)
hc.SetLineColor(ROOT.kBlack); hc.SetLineWidth(2)
hc.SetFillColorAlpha(ROOT.kAzure + 2, 0.40)
hc.SetMinimum(0)
hc.SetMaximum(max(counts.values()) * 1.55)
hc.Draw("hist")
lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.038)
lat.DrawLatex(0.17, 0.86, f"{tot_iso} bursts in {tot_acq:,} acquisitions "
                          f"({len(per_run)} runs)")
lat.DrawLatex(0.17, 0.81, f"at each run's own cut; against {tot_train} spill trains")
lat.DrawLatex(0.17, 0.76, "the SPLIT is meaningful, the TOTAL is not -- see panel 2")
keep += [hc, lat]

out = os.path.join(OUTDIR, "isolated_bursts_population.png")
c.SaveAs(out)
print(f"\nsaved {out}")
