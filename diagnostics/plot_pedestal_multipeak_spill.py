"""Are the multi-peaked pedestals BEAM-INDUCED? In-spill vs out-of-spill.

pedestal_multipeak_summary.txt established WHAT the second population looks
like (one characteristic ~7-8 ADC step, equal in ADC in both gains, so it
enters downstream of the gain split) but not WHEN it happens: the histograms
are time-integrated, so they cannot say whether a cell's baseline jumps while
beam is in the detector or all the time.

This splits the fill by beam. The acquisitions of a run separate cleanly by
their total occupancy -- median ~4 hits when quiet against a 500-2000 hit
bump, an almost empty valley between ~50 and ~500, and the high-occupancy
acquisitions arrive in contiguous blocks of ~40 separated by long quiet
stretches, i.e. the SPS spill structure. PedestalMipCalibrator's SpillSelect
property cuts on exactly that (see its doc), so `in` and `out` grids are two
runs of the SAME fill over the SAME chunks.

TWO THINGS MAKE THE COMPARISON HONEST, and both are needed:

  1. NSlabsHit=1.  The production fill's coincidence tagger (NSlabsHit=8) is
     itself a beam trigger: it wants >=8 slabs lit at the same BCID, so 99.1%
     of the SCA slices it accepts already come from in-spill acquisitions and
     an "out" grid built with it is empty. Both grids here are filled with the
     coincidence relaxed so that ONLY the beam differs between them.

  2. The same cells at the same detection power.  Raw fractions are not
     comparable -- TSpectrum resolves a shoulder that Poisson noise hides at
     lower counts, so more statistics MANUFACTURES multi-peak cells (the
     argument is spelled out in compare_multipeak_stat_matched.py). Every
     number below is computed on the INTERSECTION of cells clearing N entries
     on both sides, with each histogram resampled to exactly N entries.

The pedestal RMS comparison is the one that needs no peak-finding at all, and
is therefore the result to trust if TSpectrum's settings are ever questioned.

SCOPE: out of spill a chip is barely triggered, so essentially only SCA0
fills. The common cell set is ~96% SCA0 and these numbers speak about SCA0.
Since the multi-peak fraction rises with SCA number, they are a LOWER bound on
the whole-detector effect, not an average over it.

Usage:  IN=<in_grid.root> OUT=<out_grid.root> [NS=200,500,1000] [GAIN=low] \
        [SEED=0] [SCA=0] python3 diagnostics/plot_pedestal_multipeak_spill.py [outdir]

SCA=0 restricts everything to SCA0 and writes to *_sca0.png. Worth running:
on a pure SCA0 set the ratio is 26.8x rather than 14.5x, because the SCA1-4
sliver that survives out of spill carries multi-peak cells of its own.

Build the two grids with (per chunk group, then hadd as usual):
  CALIB_MODE=Fill CALIB_NSLABS_HIT=1 CALIB_SPILL_SELECT=in|out \
  CALIB_INPUT_FILES=... CALIB_OUTPUT_HISTOGRAM_FILE=... \
      k4run gaudi_source/options/run_pedestal_mip.py
"""
import os
import re
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUT, exist_ok=True)

IN_PATH = os.environ["IN"]
OUT_PATH = os.environ["OUT"]
GAIN = os.environ.get("GAIN", "low")
NS = [int(x) for x in os.environ.get("NS", "200,500,1000").split(",")]
SEED = int(os.environ.get("SEED", "0"))
MIN_ENTRIES = int(os.environ.get("MIN_ENTRIES", "200"))

# Same peak-search parameters as plot_pedestal_multipeak.py, so these numbers
# sit on the same scale as the ones in the tables.
SIGMA, THRESH = 2.0, 0.30

rng = np.random.default_rng(SEED)
sp = ROOT.TSpectrum(10)
COL = {"in": ROOT.kOrange + 7, "out": ROOT.kAzure + 2}
LBL = {"in": "in spill", "out": "out of spill"}


def cell_entries(f):
    prefix = f"ped_{GAIN}_"
    return {k.GetName()[len(prefix):]: k.ReadObj().GetEntries()
            for k in f.GetListOfKeys() if k.GetName().startswith(prefix)}


def resample(h, n):
    """A fresh histogram with exactly n entries drawn from h's shape.

    Multinomial over the bin contents, as in compare_multipeak_stat_matched.py:
    keeps the pedestal's shape and any genuine second population, fixes only
    the count TSpectrum gets to work with.
    """
    nb = h.GetNbinsX()
    c = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    if c.sum() <= 0:
        return None
    drawn = rng.multinomial(n, c / c.sum())
    hs = ROOT.TH1F(f"{h.GetName()}_rs{n}", "", nb,
                   h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    hs.SetDirectory(0)
    for i, v in enumerate(drawn, start=1):
        hs.SetBinContent(i, float(v))
    return hs


files = {"in": ROOT.TFile.Open(IN_PATH), "out": ROOT.TFile.Open(OUT_PATH)}
for tag, f in files.items():
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {tag} grid")
entries = {tag: cell_entries(f) for tag, f in files.items()}

# SCA=0 (or 0,1,...) restricts every number and plot below to those SCAs and
# tags the output files _sca0. Worth doing explicitly even though the common
# set is already ~96% SCA0: it removes the last doubt that the effect is an
# SCA-to-SCA pedestal difference rather than one SCA jumping between levels,
# and it keeps the higher-N rows (which drift to 75-82% SCA0) on one SCA too.
SCAS = os.environ.get("SCA", "")
if SCAS:
    want = {int(x) for x in SCAS.split(",")}
    entries = {tag: {c: v for c, v in e.items() if int(c.split("_sca")[1]) in want}
               for tag, e in entries.items()}
    print(f"restricted to SCA {sorted(want)}")
TAG = f"_sca{SCAS.replace(',', '')}" if SCAS else ""

print(f"gain={GAIN}  seed={SEED}  peak search: sigma={SIGMA} thresh={THRESH}")
for tag in ("in", "out"):
    e = np.array(list(entries[tag].values()))
    print(f"  {LBL[tag]:<13}: {len(e):>8,} cells present, "
          f"{int(e.sum()):>13,} entries, {int((e >= MIN_ENTRIES).sum()):>7,} cells >={MIN_ENTRIES}")

# --- raw fractions, each grid on its own cells (NOT comparable; shown because
#     that is the number the naive split reports, and it is 14x off) ---
raw = {}
for tag in ("in", "out"):
    f = files[tag]
    n_multi = n_tot = 0
    for cell, e in entries[tag].items():
        if e < MIN_ENTRIES:
            continue
        n_tot += 1
        if sp.Search(f.Get(f"ped_{GAIN}_{cell}"), SIGMA, "nobackground", THRESH) >= 2:
            n_multi += 1
    raw[tag] = (n_multi, n_tot)
    print(f"  RAW {LBL[tag]:<13}: {n_multi:>7,}/{n_tot:<7,} = {100 * n_multi / max(n_tot, 1):5.2f}%  "
          f"(own cell set, own statistics)")

# --- stat-matched: identical cells, identical detection power ---
rows = []
for n in NS:
    common = sorted({c for c, v in entries["in"].items() if v >= n} &
                    {c for c, v in entries["out"].items() if v >= n})
    if not common:
        print(f"N={n}: no common cells, skipping")
        continue
    per = {}
    for tag in ("in", "out"):
        f = files[tag]
        multi = 0
        for cell in common:
            hs = resample(f.Get(f"ped_{GAIN}_{cell}"), n)
            if hs is not None and sp.Search(hs, SIGMA, "nobackground", THRESH) >= 2:
                multi += 1
        per[tag] = multi
    rows.append((n, len(common), per))
    scas = [int(c.split("_sca")[1]) for c in common]
    frac0 = 100 * sum(s == 0 for s in scas) / len(scas)
    print(f"N={n:<5} common={len(common):>6,} (SCA0 {frac0:4.1f}%)   "
          f"in {100 * per['in'] / len(common):5.2f}%   out {100 * per['out'] / len(common):5.2f}%   "
          f"ratio {per['in'] / max(per['out'], 1):5.1f}x")

if not rows:
    raise SystemExit("no common cells at any N -- nothing to compare")

# --- RMS on the same cells: no peak-finding involved ---
common0 = sorted({c for c, v in entries["in"].items() if v >= MIN_ENTRIES} &
                 {c for c, v in entries["out"].items() if v >= MIN_ENTRIES})
rms = {tag: np.array([files[tag].Get(f"ped_{GAIN}_{c}").GetRMS() for c in common0])
       for tag in ("in", "out")}
d = rms["in"] - rms["out"]
print(f"\npedestal RMS on the same {len(common0):,} cells (no peak-finding):")
for tag in ("in", "out"):
    print(f"  {LBL[tag]:<13}: median {np.median(rms[tag]):5.2f} ADC  "
          f"IQR [{np.percentile(rms[tag], 25):.2f}, {np.percentile(rms[tag], 75):.2f}]")
print(f"  median widening in spill: {np.median(d):+.2f} ADC "
      f"(x{np.median(rms['in']) / np.median(rms['out']):.2f});  "
      f"cells wider in spill: {100 * (d > 0).mean():.1f}%")

# A representative cell for the overlay: multi-peak in spill, single out of it,
# and with out-of-spill statistics near the median so it is not one of the
# pathological always-noisy channels.
#
# This search MUST happen before the canvas exists: TSpectrum::Search draws the
# histogram it searched (plus its peak markers) into the current pad, so running
# it between two Draw() calls silently overwrites a finished panel.
med_out = np.median([entries["out"][c] for c in common0])
best = None
for cell in common0:
    hi, ho = files["in"].Get(f"ped_{GAIN}_{cell}"), files["out"].Get(f"ped_{GAIN}_{cell}")
    if sp.Search(hi, SIGMA, "nobackground", THRESH) < 2:
        continue
    if sp.Search(ho, SIGMA, "nobackground", THRESH) >= 2:
        continue
    score = abs(entries["out"][cell] - med_out)
    if best is None or score < best[0]:
        best = (score, cell)

# --------------------------------------------------------------- the figure ---
c = ROOT.TCanvas("c", "multipeak in vs out of spill", 2000, 1600)
c.Divide(2, 2, 0.008, 0.012)
keep = []

# 1. stat-matched fraction vs equalised N
c.cd(1)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
frame = ROOT.TH1F("frame", "multi-peak fraction, identical cells and statistics;"
                           "entries per cell (equalised);multi-peak fraction",
                  len(rows), 0, len(rows))
for i, (n, ncom, _) in enumerate(rows, start=1):
    frame.GetXaxis().SetBinLabel(i, f"{n}  ({ncom:,} cells)")
frame.SetMinimum(0)
frame.SetMaximum(1.45 * max(p[tag] / ncom for _, ncom, p in rows for tag in p))
frame.Draw("axis")
keep.append(frame)
leg = ROOT.TLegend(0.55, 0.74, 0.93, 0.88)
leg.SetBorderSize(0)
for tag in ("in", "out"):
    g = ROOT.TGraph(len(rows))
    for i, (_, ncom, p) in enumerate(rows):
        g.SetPoint(i, i + 0.5, p[tag] / ncom)
    g.SetLineColor(COL[tag])
    g.SetMarkerColor(COL[tag])
    g.SetLineWidth(3)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.6)
    g.Draw("PL same")
    leg.AddEntry(g, LBL[tag], "lp")
    keep.append(g)
leg.Draw()
keep.append(leg)

# 2. RMS distributions, same cells
c.cd(2)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
ROOT.gPad.SetLogy()
hr = {}
for tag in ("in", "out"):
    h = ROOT.TH1F(f"rms_{tag}", "pedestal RMS on the same cells "
                                "(no peak-finding);pedestal RMS (ADC);cells", 100, 0, 10)
    for v in rms[tag]:
        h.Fill(v)
    h.SetLineColor(COL[tag])
    h.SetLineWidth(3)
    hr[tag] = h
    keep.append(h)
hr["in"].SetMaximum(3 * max(hr["in"].GetMaximum(), hr["out"].GetMaximum()))
hr["in"].Draw("hist")
hr["out"].Draw("hist same")
leg2 = ROOT.TLegend(0.52, 0.72, 0.93, 0.88)
leg2.SetBorderSize(0)
for tag in ("in", "out"):
    leg2.AddEntry(hr[tag], f"{LBL[tag]} (med {np.median(rms[tag]):.2f} ADC)", "l")
leg2.Draw()
keep.append(leg2)

# 3. per-cell RMS difference -- paired, so it removes the cell-to-cell spread
c.cd(3)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
hd = ROOT.TH1F("dr", "per-cell RMS difference (paired, same cell both sides);"
                     "RMS(in spill) - RMS(out of spill)  [ADC];cells", 120, -3, 6)
for v in d:
    hd.Fill(v)
hd.SetLineColor(ROOT.kBlack)
hd.SetLineWidth(3)
hd.SetFillColorAlpha(ROOT.kOrange + 7, 0.35)
hd.Draw("hist")
zero = ROOT.TLine(0, 0, 0, hd.GetMaximum() * 1.05)
zero.SetLineStyle(2)
zero.SetLineWidth(2)
zero.Draw()
txt = ROOT.TLatex()
txt.SetNDC()
txt.SetTextSize(0.038)
txt.DrawLatex(0.52, 0.84, f"median {np.median(d):+.2f} ADC")
txt.DrawLatex(0.52, 0.79, f"{100 * (d > 0).mean():.1f}% of cells wider in spill")
keep += [hd, zero, txt]

# 4. the representative cell picked above
c.cd(4)
ROOT.gPad.SetMargin(0.12, 0.05, 0.12, 0.10)
if best is None:
    print("no representative single-out/multi-in cell found")
else:
    cell = best[1]
    hi = files["in"].Get(f"ped_{GAIN}_{cell}").Clone("ex_in")
    ho = files["out"].Get(f"ped_{GAIN}_{cell}").Clone("ex_out")
    for h in (hi, ho):
        h.SetDirectory(0)
        h.Scale(1.0 / h.Integral())
    lo = min(hi.GetMean(), ho.GetMean()) - 22
    hi.GetXaxis().SetRangeUser(lo, lo + 55)
    hi.SetTitle(f"representative cell {cell};pedestal ADC ({GAIN} gain);fraction of entries")
    for tag, h in (("in", hi), ("out", ho)):
        h.SetLineColor(COL[tag])
        h.SetLineWidth(3)
    hi.SetMaximum(1.4 * max(hi.GetMaximum(), ho.GetMaximum()))
    hi.Draw("hist")
    ho.Draw("hist same")
    leg4 = ROOT.TLegend(0.5, 0.72, 0.93, 0.88)
    leg4.SetBorderSize(0)
    leg4.AddEntry(hi, f"in spill ({int(entries['in'][cell]):,} entries)", "l")
    leg4.AddEntry(ho, f"out of spill ({int(entries['out'][cell]):,})", "l")
    leg4.Draw()
    keep += [hi, ho, leg4]
    print(f"\nrepresentative cell drawn: {cell}")

path = os.path.join(OUT, f"pedestal_multipeak_spill_{GAIN}{TAG}.png")
c.SaveAs(path)
print(f"saved {path}")

# ------------------------------------------------- second figure: six cells ---
# One cell per slab, drawn from the BULK (out-of-spill entries in [200,600], so
# not the pathological always-noisy chips), each multi-peak in spill and single
# out of it. This is the picture the fractions above are a summary of.
#
# Counter-examples exist and are worth knowing: cells on the chips that
# self-trigger hardest with no beam (s6_c11 above all, 21k out-of-spill entries)
# overlay almost exactly in/out, same satellite comb on both. They are the
# intrinsic non-beam population, ~0.85% of the common cells.
cands = []
for cell in common0:
    eo = entries["out"][cell]
    if not (200 <= eo <= 600):
        continue
    if sp.Search(files["in"].Get(f"ped_{GAIN}_{cell}"), SIGMA, "nobackground", THRESH) < 2:
        continue
    if sp.Search(files["out"].Get(f"ped_{GAIN}_{cell}"), SIGMA, "nobackground", THRESH) >= 2:
        continue
    slab = int(re.match(r"s(\d+)_", cell).group(1))
    cands.append((cell, entries["in"][cell], eo, slab))

seen, sel = set(), []
for cell, ei, eo, slab in sorted(cands, key=lambda x: -x[1]):
    if slab in seen:
        continue
    seen.add(slab)
    sel.append((cell, ei, eo))
    if len(sel) == 6:
        break

if sel:
    c2 = ROOT.TCanvas("c2", "representative cells", 1800, 1100)
    c2.Divide(3, 2, 0.006, 0.010)
    keep2 = []
    for i, (cell, ei, eo) in enumerate(sel, start=1):
        c2.cd(i)
        hi = files["in"].Get(f"ped_{GAIN}_{cell}").Clone(f"ov_i{i}")
        ho = files["out"].Get(f"ped_{GAIN}_{cell}").Clone(f"ov_o{i}")
        for h in (hi, ho):
            h.SetDirectory(0)
            h.Scale(1.0 / h.Integral())
        lo = min(hi.GetMean(), ho.GetMean()) - 20
        hi.GetXaxis().SetRangeUser(lo, lo + 50)
        hi.SetLineColor(COL["in"])
        hi.SetLineWidth(3)
        ho.SetLineColor(COL["out"])
        ho.SetLineWidth(3)
        hi.SetTitle(f"{cell};pedestal ADC ({GAIN} gain);fraction of entries")
        hi.SetMaximum(1.45 * max(hi.GetMaximum(), ho.GetMaximum()))
        hi.Draw("hist")
        ho.Draw("hist same")
        lg = ROOT.TLegend(0.50, 0.72, 0.95, 0.88)
        lg.SetBorderSize(0)
        lg.AddEntry(hi, f"in spill ({int(ei):,})", "l")
        lg.AddEntry(ho, f"out of spill ({int(eo):,})", "l")
        lg.Draw()
        keep2 += [hi, ho, lg]
    path2 = os.path.join(OUT, f"pedestal_multipeak_spill_cells_{GAIN}{TAG}.png")
    c2.SaveAs(path2)
    print(f"saved {path2}  ({len(cands):,} candidate cells, one per slab)")
