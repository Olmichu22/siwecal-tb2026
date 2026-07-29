"""Multi-peak pedestal fraction with the STATISTICS EQUALISED across thresholds.

The raw numbers in pedestal_multipeak_summary.txt are not comparable between
thresholds, and the summary says so: th220 quotes 20.6% multi-peak against
th230's 13.2%, but th230 has far fewer runs, loses 46% of its cells to the
<200-entry cut, and what survives is biased toward the SCAs that fill first.

Two distinct biases are folded into that gap:

  1. CELL SET.  Each threshold scans whatever cells cleared its own cut, so
     the three fractions are averages over different regions of the detector
     and different SCA mixes. Since the multi-peak fraction rises steeply with
     SCA number, a sample biased toward low SCA reports a lower fraction for
     free.

  2. DETECTION POWER.  TSpectrum needs counts to resolve a second population.
     A cell with 5000 entries shows a shoulder that the same cell with 250
     entries hides in Poisson noise. More statistics therefore MANUFACTURES
     multi-peak cells, and th220 has more statistics.

Both are removed here by comparing like with like:

  - the cell set is the INTERSECTION: only (slab, chip, channel, SCA) cells
    that clear N entries in EVERY threshold are used, so the geometry and SCA
    mix are identical by construction;
  - every histogram is RESAMPLED to exactly N entries by a multinomial draw
    over its own bin contents, so each threshold is peak-searched with the
    same detection power regardless of how much data it really had.

Scanning N then separates the two readings: if the gap survives at every N,
it is physical; if it closes as N is equalised, it was a statistics artefact.

Usage:
    TH230=<merged.root> TH220=<merged.root> [TH210=<merged.root>] \
    [NS=200,500,1000,2000] [GAIN=low] [SEED=0] \
    python3 diagnostics/compare_multipeak_stat_matched.py [outdir]
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUT, exist_ok=True)

GAIN = os.environ.get("GAIN", "low")
NS = [int(x) for x in os.environ.get("NS", "200,500,1000,2000").split(",")]
SEED = int(os.environ.get("SEED", "0"))

# Same peak-search parameters as plot_pedestal_multipeak.py, so the numbers
# here sit on the same scale as the ones in the summary.
SIGMA, THRESH = 2.0, 0.30

SOURCES = [(th, os.environ[f"TH{th}"]) for th in ("210", "220", "230")
           if os.environ.get(f"TH{th}")]
if len(SOURCES) < 2:
    raise SystemExit("need at least two of TH210/TH220/TH230 to compare")

rng = np.random.default_rng(SEED)


def cell_entries(path):
    """Entries per ped_<gain>_* cell, keyed by the cell name suffix."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {path}")
    prefix = f"ped_{GAIN}_"
    out = {}
    for k in f.GetListOfKeys():
        name = k.GetName()
        if name.startswith(prefix):
            out[name[len(prefix):]] = k.ReadObj().GetEntries()
    f.Close()
    return out


def resample(h, n):
    """A fresh histogram holding exactly n entries drawn from h's shape.

    Multinomial over the bin contents: preserves the pedestal's shape and its
    genuine second population, while fixing the count that TSpectrum sees.
    Underflow/overflow are ignored -- the pedestal window contains the peaks.
    """
    nb = h.GetNbinsX()
    c = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    tot = c.sum()
    if tot <= 0:
        return None
    drawn = rng.multinomial(n, c / tot)
    hs = ROOT.TH1F(f"{h.GetName()}_rs", "", nb,
                   h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    hs.SetDirectory(0)
    for i, v in enumerate(drawn, start=1):
        hs.SetBinContent(i, float(v))
    return hs


print(f"gain={GAIN}  seed={SEED}  peak search: sigma={SIGMA} thresh={THRESH}")
print("counting entries per cell ...")
entries = {th: cell_entries(path) for th, path in SOURCES}
for th, e in entries.items():
    print(f"  th{th}: {len(e):,} cells present")

files = {th: ROOT.TFile.Open(path) for th, path in SOURCES}
sp = ROOT.TSpectrum(10)
rows = []

for n in NS:
    # Cells clearing n entries in EVERY threshold: identical geometry and SCA
    # mix for all of them.
    common = None
    for th in entries:
        ok = {c for c, v in entries[th].items() if v >= n}
        common = ok if common is None else (common & ok)
    common = sorted(common)
    if not common:
        print(f"N={n}: no common cells, skipping")
        continue

    per_th = {}
    for th in entries:
        f = files[th]
        multi = 0
        per_sca = np.zeros((15, 2))
        for cell in common:
            h = f.Get(f"ped_{GAIN}_{cell}")
            hs = resample(h, n)
            if hs is None:
                continue
            npk = sp.Search(hs, SIGMA, "nobackground", THRESH)
            sca = int(cell.split("_sca")[1])
            per_sca[sca, 0] += 1
            if npk >= 2:
                multi += 1
                per_sca[sca, 1] += 1
        per_th[th] = (multi, len(common), per_sca)

    rows.append((n, len(common), per_th))
    line = f"N={n:<5} common cells={len(common):>7,}   "
    line += "   ".join(f"th{th}: {100 * m / max(t, 1):5.2f}%"
                       for th, (m, t, _) in per_th.items())
    print(line)

if not rows:
    raise SystemExit("nothing to plot")

# One canvas: fraction vs equalised N, one line per threshold.
c = ROOT.TCanvas("c", "stat-matched multipeak", 1400, 900)
c.SetMargin(0.12, 0.05, 0.12, 0.10)
keep = []
colors = {"210": ROOT.kAzure + 1, "220": ROOT.kOrange + 7, "230": ROOT.kRed + 1}
frame = ROOT.TH1F("frame", f"multi-peak fraction at matched statistics "
                           f"({GAIN} gain, identical cell set);"
                           f"entries per cell (equalised);multi-peak fraction",
                  len(rows), 0, len(rows))
for i, (n, _, _) in enumerate(rows, start=1):
    frame.GetXaxis().SetBinLabel(i, str(n))
allv = [100 * m / max(t, 1) for _, _, pt in rows for m, t, _ in pt.values()]
frame.SetMinimum(0)
frame.SetMaximum(max(allv) * 1.35 / 100)
frame.Draw("axis")
keep.append(frame)

leg = ROOT.TLegend(0.62, 0.72, 0.93, 0.88)
leg.SetBorderSize(0)
for th in entries:
    g = ROOT.TGraph(len(rows))
    for i, (n, _, pt) in enumerate(rows):
        m, t, _ = pt[th]
        g.SetPoint(i, i + 0.5, m / max(t, 1))
    g.SetLineColor(colors[th])
    g.SetMarkerColor(colors[th])
    g.SetLineWidth(3)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.5)
    g.Draw("lp same")
    leg.AddEntry(g, f"th{th}", "lp")
    keep.append(g)
leg.Draw()
keep.append(leg)

out = os.path.join(OUT, f"multipeak_stat_matched_{GAIN}.png")
c.SaveAs(out)
print(f"saved {out}")

# The per-SCA breakdown at the largest N, where the equalisation is strongest.
n, ncom, per_th = rows[-1]
print(f"\nper-SCA multi-peak fraction at N={n} ({ncom:,} common cells)")
print("  SCA  " + "  ".join(f"th{th}" for th in entries))
for sca in range(15):
    vals = []
    for th in entries:
        _, _, ps = per_th[th]
        vals.append(f"{100 * ps[sca, 1] / ps[sca, 0]:5.1f}" if ps[sca, 0] else "    -")
    if any(v.strip() != "-" for v in vals):
        print(f"  {sca:>3}  " + "  ".join(vals))
