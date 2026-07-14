"""Does anchoring the low gain to the high gain give a better energy distribution?

hit_energy now takes saturated hits (raw adc_high >= 1900) from the low-gain branch
ANCHORED to the high gain -- adc_high_equiv = (adc_low - ped_lg - c)/k, energy =
adc_high_equiv / MIP_hg -- so MIP_lg never enters. hit_energy_nocalib keeps the old
value, the low gain on its own MIP_lg scale, on the SAME events.

MIP_lg is inflated 1.5x-2x (a low-gain MIP is 2.08 ADC against ~1 ADC of noise, S/N
~ 2, so the Landau fit lands at 3.2-4.5 instead). Dividing by an inflated MIP
UNDER-reads the energy. So the anchored version should push the saturated hits up,
and the prediction is specific:

  - the mean goes UP,
  - the resolution sigma/mu gets BETTER, because the old way scattered every
    saturated hit by a channel-dependent wrong factor, and that scatter was pure
    noise added to the energy sum,
  - events with NO saturated hit do not move at all. That is the control: if they
    move, something other than the intended change is going on.

Panels:
  1. sum_energy, both ways, whole sample.
  2. Only events WITH at least one saturated hit -- where the change actually acts.
  3. THE CONTROL: only events with NO saturated hit. Must be identical.
  4. Per-hit: the two energies against each other for saturated hits.

Usage:  RUN=TB2026CERN_run_000012 python3 diagnostics/plot_energy_calib.py [outdir]
"""
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.environ.get("RUN", "TB2026CERN_run_000012")
TH = os.environ.get("TH", "230")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)

ECAL = os.environ.get(
    "ECAL",
    f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Reconstruction/{RUN}/ecal_{RUN}.root")

f = ROOT.TFile.Open(ECAL)
if not f or f.IsZombie():
    raise SystemExit(f"ERROR: cannot open {ECAL}")
t = f.Get("ecal")
if not t:
    raise SystemExit(f"ERROR: no 'ecal' tree in {ECAL}")
n = t.GetEntries()
print(f"[{RUN}] {n:,} reconstructed events")

if not t.GetBranch("hit_energy_nocalib"):
    raise SystemExit("ERROR: no hit_energy_nocalib branch -- rebuild and re-run the event builder")

# Per event: the two sums, and whether any hit was saturated. A hit is saturated
# exactly when the two energies differ (below the switch both read the high gain and
# are bit-identical), which is more robust than re-deriving the threshold here.
e_new, e_old, has_sat, n_sat = [], [], [], []
for i in range(n):
    t.GetEntry(i)
    nc = t.nhit_chan
    a = np.frombuffer(t.hit_energy, dtype=np.float32, count=nc).astype(float)
    b = np.frombuffer(t.hit_energy_nocalib, dtype=np.float32, count=nc).astype(float)
    d = np.abs(a - b) > 1e-6
    e_new.append(float(t.sum_energy))
    e_old.append(float(t.sum_energy_nocalib))
    has_sat.append(bool(d.any()))
    n_sat.append(int(d.sum()))
    if i and i % 20000 == 0:
        print(f"  ... {i:,}/{n:,}")

e_new = np.array(e_new)
e_old = np.array(e_old)
has_sat = np.array(has_sat)
n_sat = np.array(n_sat)
print(f"  events with >=1 saturated hit: {int(has_sat.sum()):,} "
      f"({100 * has_sat.mean():.1f}%)   saturated hits: {int(n_sat.sum()):,}")

c = ROOT.TCanvas("c", "energy", 1700, 1200)
c.Divide(2, 2, 0.005, 0.005)
keep = []
BLUE, RED = ROOT.kBlue + 1, ROOT.kRed + 1


def stat(v):
    if len(v) == 0:
        return float("nan"), float("nan"), float("nan")
    mu, sd = float(v.mean()), float(v.std(ddof=1)) if len(v) > 1 else 0.0
    return mu, sd, (sd / mu if mu else float("nan"))


def compare(pad_id, mask, title, note):
    pad = c.cd(pad_id)
    pad.SetMargin(0.13, 0.05, 0.13, 0.10)
    a, b = e_new[mask], e_old[mask]
    if len(a) == 0:
        return
    lo = 0.0
    hi = float(np.percentile(np.concatenate([a, b]), 99.5)) * 1.15
    hs = []
    for arr, lab, col in ((b, "hit_energy_nocalib (low gain on MIP_lg)", RED),
                          (a, "hit_energy (low gain anchored to high gain)", BLUE)):
        h = ROOT.TH1F(f"h{pad_id}_{lab[:6]}", f"{title};sum_energy  [MIP];events", 90, lo, hi)
        h.SetDirectory(0)
        for v in arr:
            h.Fill(v)
        h.SetLineColor(col)
        h.SetLineWidth(3)
        h.SetFillColorAlpha(col, 0.18)
        hs.append((h, lab, arr))
        keep.append(h)
    ymax = max(h.GetMaximum() for h, _, _ in hs)
    leg = ROOT.TLegend(0.42, 0.58, 0.96, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.80)
    leg.SetTextSize(0.024)
    for i, (h, lab, arr) in enumerate(hs):
        h.GetYaxis().SetRangeUser(0, ymax * 1.45)
        h.Draw("hist" if i == 0 else "hist same")
        mu, sd, res = stat(arr)
        leg.AddEntry(h, f"{lab}", "f")
        leg.AddEntry(0, f"    #mu = {mu:.1f}   #sigma = {sd:.1f}   #sigma/#mu = {res:.4f}", "")
    leg.AddEntry(0, "", "")
    for ln in note:
        leg.AddEntry(0, ln, "")
    leg.AddEntry(0, f"N = {len(a):,} events", "")
    leg.Draw()
    keep.append(leg)


mu_n, sd_n, r_n = stat(e_new)
mu_o, sd_o, r_o = stat(e_old)
compare(1, np.ones(len(e_new), dtype=bool), f"{RUN}: all events",
        [f"mean shifts by {100 * (mu_n / mu_o - 1):+.2f}%",
         f"resolution {r_o:.4f} -> {r_n:.4f} "
         f"({100 * (r_n / r_o - 1):+.1f}%)"])

if has_sat.any():
    a, b = e_new[has_sat], e_old[has_sat]
    ms, _, rs = stat(a)
    mo, _, ro = stat(b)
    compare(2, has_sat, f"{RUN}: events WITH a saturated hit -- where the change acts",
            [f"mean shifts by {100 * (ms / mo - 1):+.2f}%",
             f"resolution {ro:.4f} -> {rs:.4f}"])

ctrl = ~has_sat
if ctrl.any():
    dmax = float(np.abs(e_new[ctrl] - e_old[ctrl]).max())
    compare(3, ctrl, f"{RUN}: CONTROL -- events with NO saturated hit",
            ["these MUST be identical",
             f"largest difference: {dmax:.2e} MIP"])

# ---- 4. per-hit, the two energies against each other ---------------------
pad = c.cd(4)
pad.SetMargin(0.13, 0.14, 0.13, 0.10)
pad.SetLogz()
xs, ys = [], []
for i in range(min(n, 40000)):
    t.GetEntry(i)
    nc = t.nhit_chan
    a = np.frombuffer(t.hit_energy, dtype=np.float32, count=nc).astype(float)
    b = np.frombuffer(t.hit_energy_nocalib, dtype=np.float32, count=nc).astype(float)
    d = np.abs(a - b) > 1e-6
    if d.any():
        xs.append(b[d])
        ys.append(a[d])
if xs:
    X = np.concatenate(xs)
    Y = np.concatenate(ys)
    hi = float(np.percentile(np.concatenate([X, Y]), 99.5)) * 1.1
    h4 = ROOT.TH2F("h4", f"{RUN}: saturated hits, the two calibrations against each other;"
                         f"hit_energy_nocalib  [MIP];hit_energy (anchored)  [MIP]",
                   90, 0, hi, 90, 0, hi)
    for x, y in zip(X, Y):
        h4.Fill(x, y)
    h4.Draw("colz")
    diag = ROOT.TLine(0, 0, hi, hi)
    diag.SetLineColor(ROOT.kBlack)
    diag.SetLineWidth(3)
    diag.SetLineStyle(2)
    diag.Draw()
    ratio = float(np.median(Y / np.where(X > 0, X, np.nan)))
    leg4 = ROOT.TLegend(0.16, 0.72, 0.62, 0.89)
    leg4.SetBorderSize(0)
    leg4.SetFillColorAlpha(ROOT.kWhite, 0.80)
    leg4.SetTextSize(0.026)
    leg4.AddEntry(diag, "equal (no change)", "l")
    leg4.AddEntry(0, f"median ratio new/old = {ratio:.2f}", "")
    leg4.AddEntry(0, f"{len(X):,} saturated hits", "")
    leg4.Draw()
    keep += [h4, diag, leg4]

out = os.path.join(OUT, f"energy_calib_{RUN.split('_')[-1]}.png")
c.SaveAs(out)

print(f"\n{RUN}: anchored low gain vs MIP_lg low gain\n")
print(f"  {'sample':34s} {'N':>9s} {'mu':>9s} {'sigma':>8s} {'sigma/mu':>9s}")
print("  " + "-" * 74)
for lab, mask in (("all events", np.ones(len(e_new), dtype=bool)),
                  ("with a saturated hit", has_sat),
                  ("CONTROL: no saturated hit", ctrl)):
    if not mask.any():
        continue
    for tag, arr in (("old (MIP_lg)", e_old[mask]), ("new (anchored)", e_new[mask])):
        mu, sd, r = stat(arr)
        print(f"  {lab + '  ' + tag:34s} {int(mask.sum()):9,} {mu:9.1f} {sd:8.1f} {r:9.4f}")
    print()
print(f"  overall: mu {mu_o:.1f} -> {mu_n:.1f} ({100 * (mu_n / mu_o - 1):+.2f}%), "
      f"sigma/mu {r_o:.4f} -> {r_n:.4f} ({100 * (r_n / r_o - 1):+.1f}%)")
print(f"\nsaved {out}")
