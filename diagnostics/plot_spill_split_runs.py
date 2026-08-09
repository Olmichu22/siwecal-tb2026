"""Does the in-spill / out-of-spill result hold across runs, or was it run_000060?

Everything in section 7 of diagnostics/pedestal_multipeak_tables.txt comes from a
single th220 run. Two claims in particular are exposed to that:

  1. the size of the beam-induced excess (26.8x on SCA0 at matched N=200), and
  2. "the intrinsic, non-beam multi-peak tail is TWO WHOLE CHIPS, s13_c9 and
     s6_c11" -- which is either a property of those two chips, in which case
     the same two come back run after run, or a coincidence of one run, in
     which case the two-populations reading needs rewriting.

WHAT IT ANSWERS (9 th220 runs, SCA0, N=200, low gain):

  1. HOLDS, with the emphasis moved. The IN-SPILL fraction is the stable
     quantity: 18-24% in the six runs with enough cells. The ratio is not --
     2.2x to 108x -- because what varies between runs is the out-of-spill
     fraction (0.22% to 8.30%), i.e. the denominator. Quote the in-spill
     fraction and the RMS widening (median +1.00 ADC, 87-98% of cells wider in
     spill); quote the ratio only with its run.
  2. DOES NOT HOLD as stated. s13_c9 does not appear at all once "intrinsic" is
     measured rather than proxied by entry count -- it is a chip that
     self-triggers often, not one that is multi-peak without beam. s6_c11 shows
     up in 3 runs of 8, and every other chip in at most 2. run_000254 is a case
     of its own: 8.30% of its cells are multi-peak out of spill, an order of
     magnitude above every other run, spread over 81 chips.

Two runs are too short to say anything: run_000061 and run_000087 have 189 and
127 common cells against 8,000-15,000 elsewhere, and run_000087 even comes out
with a NARROWER in-spill RMS, which is a statistics artefact rather than a
reversal. They are marked with * on the axes. run_000143 has no common cells at
all at N=200.

This runs the section-7 comparison over every th220 run that has an in/out grid
pair and puts the answers side by side. Everything is on SCA0 only and at matched
N, for the reasons argued in plot_pedestal_multipeak_spill.py -- raw fractions
across runs would mostly measure how long each run is.

"INTRINSIC" MEANS MEASURED, NOT PROXIED. Section 7 identified the non-beam tail
as the cells holding at least 4,000 out-of-spill entries. That proxy does not
survive contact with other runs: it is an absolute threshold, so run_000142 --
503,707 acquisitions, and barely any beam -- returns every one of its 5,229
common cells on 85 chips as "intrinsic", which is an artefact of the run being
long rather than of any chip being special. Here a cell is intrinsic if it is
actually flagged MULTI-PEAK OUT OF SPILL at the matched N, which is the property
the two-populations argument is about and which carries its own detection-power
control. Results are reported per chip because an SCA slice is written for all
64 channels at once, making the underlying quantity a chip property.

Usage:
  GRIDS=<dir with *_ns1_in.root / *_ns1_out.root> [N=200] [GAIN=low] [SCA=0] \
      python3 diagnostics/plot_spill_split_runs.py [outdir]
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

GRIDS = os.environ["GRIDS"]
GAIN = os.environ.get("GAIN", "low")
N = int(os.environ.get("N", "200"))
SEED = int(os.environ.get("SEED", "0"))
SCAS = {int(x) for x in os.environ.get("SCA", "0").split(",") if x != ""}
# A chip counts as intrinsic in a run only if at least this many of its cells are
# multi-peak out of spill. TSpectrum has a small false-positive rate, so without
# a floor the tally fills up with single scattered cells from dozens of chips and
# says nothing about which chips are actually special.
MIN_CELLS_PER_CHIP = int(os.environ.get("MIN_CELLS_PER_CHIP", "8"))
SIGMA, THRESH = 2.0, 0.30            # same peak search as everywhere else

rng = np.random.default_rng(SEED)
sp = ROOT.TSpectrum(10)
CELL_RE = re.compile(r"^s(\d+)_c(\d+)_ch(\d+)_sca(\d+)$")


def resample(h, n):
    nb = h.GetNbinsX()
    c = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    if c.sum() <= 0:
        return None
    drawn = rng.multinomial(n, c / c.sum())
    hs = ROOT.TH1F(f"{h.GetName()}_rs", "", nb,
                   h.GetXaxis().GetXmin(), h.GetXaxis().GetXmax())
    hs.SetDirectory(0)
    for i, v in enumerate(drawn, start=1):
        hs.SetBinContent(i, float(v))
    return hs


def cells(f):
    prefix = f"ped_{GAIN}_"
    out = {}
    for k in f.GetListOfKeys():
        nm = k.GetName()
        if not nm.startswith(prefix):
            continue
        cell = nm[len(prefix):]
        m = CELL_RE.match(cell)
        if not m or int(m.group(4)) not in SCAS:
            continue
        out[cell] = k.ReadObj().GetEntries()
    return out


pairs = []
for p in sorted(glob.glob(os.path.join(GRIDS, f"*_ns1_in.root"))):
    q = p.replace("_ns1_in.root", "_ns1_out.root")
    if os.path.exists(q):
        pairs.append((os.path.basename(p).replace("_ns1_in.root", ""), p, q))
if not pairs:
    raise SystemExit(f"no in/out grid pairs under {GRIDS}")
print(f"{len(pairs)} runs with a grid pair, SCA {sorted(SCAS)}, N={N}, gain={GAIN}\n")

results = []
chip_hits = Counter()          # (slab, chip) -> in how many runs it is intrinsic
chip_per_run = {}
for run, pin, pout in pairs:
    fi, fo = ROOT.TFile.Open(pin), ROOT.TFile.Open(pout)
    if not fi or fi.IsZombie() or not fo or fo.IsZombie():
        print(f"{run}: cannot open a grid, skipping")
        continue
    ei, eo = cells(fi), cells(fo)
    common = sorted({c for c, v in ei.items() if v >= N} & {c for c, v in eo.items() if v >= N})
    if not common:
        print(f"{run}: no common cells at N={N}, skipping")
        fi.Close(); fo.Close()
        continue

    multi = {}
    multi_out_cells = []
    for tag, f in (("in", fi), ("out", fo)):
        k = 0
        for cell in common:
            hs = resample(f.Get(f"ped_{GAIN}_{cell}"), N)
            if hs is not None and sp.Search(hs, SIGMA, "nobackground", THRESH) >= 2:
                k += 1
                if tag == "out":
                    multi_out_cells.append(cell)
        multi[tag] = k
    rin = 100 * multi["in"] / len(common)
    rout = 100 * multi["out"] / len(common)

    rms_i = np.array([fi.Get(f"ped_{GAIN}_{c}").GetRMS() for c in common])
    rms_o = np.array([fo.Get(f"ped_{GAIN}_{c}").GetRMS() for c in common])
    wider = 100 * ((rms_i - rms_o) > 0).mean()

    # the intrinsic tail: which CHIPS carry cells that are multi-peak WITHOUT beam
    tail = multi_out_cells
    tchips = Counter()
    for c in tail:
        m = CELL_RE.match(c)
        tchips[(int(m.group(1)), int(m.group(2)))] += 1
    tchips = Counter({ch: n for ch, n in tchips.items() if n >= MIN_CELLS_PER_CHIP})
    for ch in tchips:
        chip_hits[ch] += 1
    chip_per_run[run] = tchips

    results.append(dict(run=run, ncommon=len(common), rin=rin, rout=rout,
                        ratio=rin / rout if rout > 0 else float("inf"),
                        med_i=float(np.median(rms_i)), med_o=float(np.median(rms_o)),
                        wider=wider, ntail=len(tail), nchips=len(tchips)))
    short = run.replace("TB2026CERN_", "").replace("eudaq_", "")
    top = ", ".join(f"s{s}_c{c}({n})" for (s, c), n in tchips.most_common(4))
    print(f"{short:<12} cells {len(common):>6,}   in {rin:5.2f}%  out {rout:5.2f}%  "
          f"ratio {results[-1]['ratio']:6.1f}x   RMS {np.median(rms_i):.2f}/{np.median(rms_o):.2f} "
          f"({wider:4.1f}% wider)   multi-peak w/o beam: {len(tail):>4} cells, "
          f"{len(tchips)} chip(s) above the floor: {top}")
    fi.Close()
    fo.Close()

if not results:
    raise SystemExit("nothing to plot")

print(f"\nchips carrying cells that are multi-peak WITHOUT beam "
      f"(>= {MIN_CELLS_PER_CHIP} such cells in the run), and in how many runs:")
for (s, ch), k in chip_hits.most_common(15):
    print(f"   s{s}_c{ch}: {k}/{len(results)} runs")

ratios = [r["ratio"] for r in results if np.isfinite(r["ratio"])]
print(f"\nratio across runs: median {np.median(ratios):.1f}x, "
      f"range {min(ratios):.1f}-{max(ratios):.1f}x")
wid = [r["med_i"] - r["med_o"] for r in results]
print(f"RMS widening across runs: median {np.median(wid):+.2f} ADC, "
      f"range {min(wid):+.2f} to {max(wid):+.2f} ADC")

# ------------------------------------------------------------------ figure
RUNNO_RE = re.compile(r"run_0*(\d+)$")


def runno(name):
    m = RUNNO_RE.search(name)
    return int(m.group(1)) if m else 0


# Sort by run number, and strip the label down to the number itself. Chopping a
# fixed "run_0000" prefix does not work: run 254 is written run_000254, with one
# zero fewer, so it came out on the axis as the full "run_000254" and ran into
# the axis title.
results.sort(key=lambda r: runno(r["run"]))
nr = len(results)
labels = [str(runno(r["run"])) for r in results]
# Two of the nine runs are far too short to mean anything (127 and 189 common
# cells against 8,000-15,000). Marking them on the axis is the difference
# between a reader seeing eight comparable runs and seeing six.
LOWSTAT = 1000
labels = [f"{l}*" if r["ncommon"] < LOWSTAT else l for l, r in zip(labels, results)]
C_IN, C_OUT = ROOT.kOrange + 7, ROOT.kAzure + 2

c = ROOT.TCanvas("c", "spill split across runs", 1800, 1300)
c.Divide(2, 2, 0.010, 0.014)
keep = []


def axis_labels(h):
    for i, lab in enumerate(labels):
        h.GetXaxis().SetBinLabel(i + 1, lab)
    h.GetXaxis().SetLabelSize(0.055)


# --- 1. multi-peak fraction in vs out, per run
c.cd(1)
ROOT.gPad.SetMargin(0.12, 0.04, 0.12, 0.10)
ROOT.gPad.SetGridy()
hi = ROOT.TH1F("hi", f"Multi-peak cells at matched N={N}, SCA0;"
                     f"th220 run  (* = fewer than {LOWSTAT:,} common cells);"
                     f"cells (%)", nr, 0, nr)
ho = ROOT.TH1F("ho", "", nr, 0, nr)
for i, r in enumerate(results):
    hi.SetBinContent(i + 1, r["rin"])
    ho.SetBinContent(i + 1, r["rout"])
axis_labels(hi)
hi.SetLineColor(C_IN); hi.SetLineWidth(3); hi.SetFillColorAlpha(C_IN, 0.35)
ho.SetLineColor(C_OUT); ho.SetLineWidth(3); ho.SetFillColorAlpha(C_OUT, 0.35)
hi.SetMinimum(0)
hi.SetMaximum(max(r["rin"] for r in results) * 1.45)
hi.Draw("hist")
ho.Draw("hist same")
l1 = ROOT.TLegend(0.60, 0.76, 0.95, 0.89)
l1.SetBorderSize(0); l1.SetFillStyle(0)
l1.AddEntry(hi, "in spill", "f")
l1.AddEntry(ho, "out of spill", "f")
l1.Draw()
keep += [hi, ho, l1]

# --- 2. the ratio, run by run
c.cd(2)
ROOT.gPad.SetMargin(0.12, 0.04, 0.12, 0.10)
ROOT.gPad.SetGridy()
hr = ROOT.TH1F("hr", "Beam-induced excess;th220 run;in-spill / out-of-spill", nr, 0, nr)
for i, r in enumerate(results):
    hr.SetBinContent(i + 1, min(r["ratio"], 999))
axis_labels(hr)
hr.SetLineColor(ROOT.kViolet + 1); hr.SetLineWidth(3)
hr.SetFillColorAlpha(ROOT.kViolet + 1, 0.30)
hr.SetMinimum(0)
hr.SetMaximum(max(ratios) * 1.35)
hr.Draw("hist")
med = ROOT.TLine(0, np.median(ratios), nr, np.median(ratios))
med.SetLineStyle(2); med.SetLineWidth(2); med.Draw()
lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.038)
lat.DrawLatex(0.16, 0.84, f"median {np.median(ratios):.1f}x over {nr} runs")
keep += [hr, med, lat]

# --- 3. RMS widening, the result that needs no peak-finding
c.cd(3)
ROOT.gPad.SetMargin(0.12, 0.04, 0.12, 0.10)
ROOT.gPad.SetGridy()
hw = ROOT.TH1F("hw", "Pedestal RMS on the same cells;th220 run;median RMS (ADC)", nr, 0, nr)
hw2 = ROOT.TH1F("hw2", "", nr, 0, nr)
for i, r in enumerate(results):
    hw.SetBinContent(i + 1, r["med_i"])
    hw2.SetBinContent(i + 1, r["med_o"])
axis_labels(hw)
hw.SetLineColor(C_IN); hw.SetLineWidth(3); hw.SetFillColorAlpha(C_IN, 0.35)
hw2.SetLineColor(C_OUT); hw2.SetLineWidth(3); hw2.SetFillColorAlpha(C_OUT, 0.35)
hw.SetMinimum(0)
hw.SetMaximum(max(r["med_i"] for r in results) * 1.45)
hw.Draw("hist")
hw2.Draw("hist same")
lat3 = ROOT.TLatex(); lat3.SetNDC(); lat3.SetTextSize(0.036)
lat3.DrawLatex(0.15, 0.84, f"median widening {np.median(wid):+.2f} ADC")
lat3.DrawLatex(0.15, 0.79, f"{np.mean([r['wider'] for r in results]):.0f}% of cells wider in spill")
keep += [hw, hw2, lat3]

# --- 4. the intrinsic chips: is it the same two every run?
c.cd(4)
ROOT.gPad.SetMargin(0.16, 0.13, 0.12, 0.10)
top_chips = [ch for ch, _ in chip_hits.most_common(12)]
hc = ROOT.TH2F("hc", "Chips multi-peak WITHOUT beam, at matched N;"
                     "th220 run;chip", nr, 0, nr, max(len(top_chips), 1), 0,
               max(len(top_chips), 1))
for i, r in enumerate(results):
    per = chip_per_run.get(r["run"], {})
    for j, ch in enumerate(top_chips):
        hc.SetBinContent(i + 1, j + 1, per.get(ch, 0))
axis_labels(hc)
for j, (s, ch) in enumerate(top_chips):
    hc.GetYaxis().SetBinLabel(j + 1, f"s{s}_c{ch}")
hc.GetYaxis().SetLabelSize(0.050)
hc.GetZaxis().SetTitle("cells multi-peak out of spill")
hc.SetMinimum(-1e-9)          # so genuine zeros are painted, not left white
hc.Draw("colz")
keep += [hc]

out = os.path.join(OUTDIR, f"spill_split_runs_{GAIN}.png")
c.SaveAs(out)
print(f"\nsaved {out}")
