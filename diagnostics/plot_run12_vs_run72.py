"""THE VALIDATION. Same beam, two thresholds -- do they now agree?

Everything else in this analysis is internal: it shows that the old calibration was
inconsistent with itself, and that the new one is consistent. Consistency is not
correctness. Nothing so far says the energy scale is RIGHT, only that it is less
wrong, and no amount of internal cross-checking can close that gap.

This can. run_000012 and run_000072 are the SAME 74 GeV electron beam:

    run_000012   74 GeV   ThresholdDAC 230
    run_000072   74 GeV   ThresholdDAC 220

The trigger threshold is a property of the DISCRIMINATOR. It decides which hits get
recorded. It has no business changing how much energy a 74 GeV electron deposits. So
the two runs MUST reconstruct to the same mean energy -- and under the old
calibration they did not, because th230's MIP table was inflated 1.47x by its own
trigger and run12 was calibrated against it.

Now both take their MIP from th220 and their low gain is anchored to their high gain,
so neither calibration depends on the threshold any more. The prediction is sharp:

    THE TWO MEANS SHOULD NOW AGREE.

If they do, the fix is validated against something outside itself. If they do not,
something is still wrong and we need to know before any of this goes on a slide.

The plot shows both, old calibration and new, so the convergence (or its absence) is
visible rather than asserted. hit_energy_nocalib carries the old value on the same
events, so no re-reconstruction is needed to draw the comparison.

Usage:  python3 diagnostics/plot_run12_vs_run72.py
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

RECO = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Reconstruction"
RUNS = [
    ("TB2026CERN_run_000012", "th230", ROOT.kRed + 1),
    ("TB2026CERN_run_000072", "th220", ROOT.kBlue + 1),
]
BEAM_GEV = 74

data = {}
for run, th, colour in RUNS:
    path = f"{RECO}/{run}/ecal_{run}.root"
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"ERROR: cannot open {path}")
    t = f.Get("ecal")
    if not t or not t.GetBranch("sum_energy_nocalib"):
        raise SystemExit(f"ERROR: {run} has no sum_energy_nocalib -- reconstruct it again")
    n = t.GetEntries()
    new = np.zeros(n)
    old = np.zeros(n)
    for i in range(n):
        t.GetEntry(i)
        new[i] = t.sum_energy
        old[i] = t.sum_energy_nocalib
    f.Close()
    data[run] = dict(th=th, colour=colour, new=new, old=old, n=n)
    print(f"[{run}] {th}: {n:,} events   "
          f"old mu={old.mean():.1f}   new mu={new.mean():.1f}")

a, b = data[RUNS[0][0]], data[RUNS[1][0]]      # a = run12 (th230), b = run72 (th220)


def agree(x, y):
    """Ratio of the two means, and how far from 1 it is."""
    r = x.mean() / y.mean()
    # error on each mean, then on the ratio
    ea = x.std(ddof=1) / np.sqrt(len(x))
    eb = y.std(ddof=1) / np.sqrt(len(y))
    er = r * np.hypot(ea / x.mean(), eb / y.mean())
    return r, er


r_old, e_old = agree(a["old"], b["old"])
r_new, e_new = agree(a["new"], b["new"])

c = ROOT.TCanvas("c", "run12 vs run72", 1700, 780)
c.Divide(2, 1, 0.005, 0.005)
keep = []


def panel(pad_id, key, title, ratio, err):
    global keep
    pad = c.cd(pad_id)
    pad.SetMargin(0.13, 0.05, 0.13, 0.10)
    allv = np.concatenate([data[r][key] for r, _, _ in RUNS])
    hi = float(np.percentile(allv, 99.5)) * 1.15
    hs = []
    for run, th, colour in RUNS:
        d = data[run]
        h = ROOT.TH1F(f"h{pad_id}_{th}", f"{title};sum_energy  [MIP];events (norm.)",
                      90, 0, hi)
        h.SetDirectory(0)
        for v in d[key]:
            h.Fill(v)
        if h.Integral():
            h.Scale(1.0 / h.Integral())
        h.SetLineColor(colour)
        h.SetLineWidth(3)
        h.SetFillColorAlpha(colour, 0.18)
        hs.append((h, run, th, d[key]))
        keep.append(h)

    ymax = max(h.GetMaximum() for h, _, _, _ in hs)
    leg = ROOT.TLegend(0.44, 0.55, 0.96, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.80)
    leg.SetTextSize(0.026)
    for i, (h, run, th, arr) in enumerate(hs):
        h.GetYaxis().SetRangeUser(0, ymax * 1.45)
        h.Draw("hist" if i == 0 else "hist same")
        mu, sd = arr.mean(), arr.std(ddof=1)
        leg.AddEntry(h, f"{run.split('_')[-1]}  ({th}, {BEAM_GEV} GeV)", "f")
        leg.AddEntry(0, f"    #mu = {mu:.0f}   #sigma/#mu = {sd / mu:.4f}", "")
    leg.AddEntry(0, "", "")
    leg.AddEntry(0, f"#mu(th230) / #mu(th220) = {ratio:.4f} #pm {err:.4f}", "")
    off = 100 * abs(ratio - 1)
    leg.AddEntry(0, f"   -> {off:.1f}% apart "
                    f"({'AGREE' if off < 2 else 'DO NOT AGREE'})", "")
    leg.Draw()
    keep.append(leg)

    # Both means, marked. Same beam -> they must land on top of each other.
    for h, run, th, arr in hs:
        ln = ROOT.TLine(arr.mean(), 0, arr.mean(), ymax * 1.10)
        ln.SetLineColor(h.GetLineColor())
        ln.SetLineWidth(3)
        ln.SetLineStyle(2)
        ln.Draw()
        keep.append(ln)


panel(1, "old", "OLD calibration: th230 MIP + MIP_lg switch", r_old, e_old)
panel(2, "new", "NEW: th220 MIP everywhere + low gain anchored", r_new, e_new)

out = os.path.join(OUT, "run12_vs_run72_validation.png")
c.SaveAs(out)

print("\n" + "=" * 76)
print(f"  THE VALIDATION: same {BEAM_GEV} GeV beam, two trigger thresholds")
print("=" * 76)
print(f"\n  {'':22s} {'run12 (th230)':>16s} {'run72 (th220)':>16s} {'ratio':>16s}")
print("  " + "-" * 72)
print(f"  {'OLD calibration':22s} {a['old'].mean():16.0f} {b['old'].mean():16.0f} "
      f"{r_old:10.4f} +/- {e_old:.4f}")
print(f"  {'NEW calibration':22s} {a['new'].mean():16.0f} {b['new'].mean():16.0f} "
      f"{r_new:10.4f} +/- {e_new:.4f}")
print(f"\n  The threshold is a property of the DISCRIMINATOR. It cannot change how much")
print(f"  energy a {BEAM_GEV} GeV electron deposits. These two ratios MUST be 1.")
print(f"\n     old:  {100 * abs(r_old - 1):5.1f}% apart")
print(f"     new:  {100 * abs(r_new - 1):5.1f}% apart")
if abs(r_new - 1) < abs(r_old - 1):
    print(f"\n  -> the gap CLOSED by a factor {abs(r_old - 1) / max(abs(r_new - 1), 1e-9):.1f}")
else:
    print(f"\n  -> the gap did NOT close. Something is still wrong; do not put this on a slide.")
print(f"\nsaved {out}")
