"""Where does the high gain stop being trustworthy? Measure it, don't eyeball it.

The adc_high:adc_low band (plot_saturation.py) shows the high-gain preamp bending
over somewhere north of 1900 ADC, but "somewhere" is not a number you can defend
in a threshold.

This is the direct test. For every hit, compute the energy BOTH gains claim:

    E_hg = (adc_high - pedestal_hg) / mip_hg
    E_lg = (adc_low  - pedestal_lg) / mip_lg

While the high gain is linear the two must agree, so E_lg/E_hg = 1. Once the high
gain saturates it under-reads, so E_hg shrinks and the ratio climbs away from 1.
**The adc_high at which the ratio departs from 1 IS the switch point** -- switching
earlier than that replaces a good high-gain measurement with a noisier low-gain
one, which is why moving the threshold DOWN to 1200 made run_000013's fitted mu
go DOWN: recovering genuinely saturated hits can only ever add energy.

Output: the ratio plane with its median profile, the 1% and 5% departure crossings
marked, and the old (1200) / new (1900) thresholds. Plus, printed: the fraction of
hits that sit between the two thresholds -- the hits whose energy the old value was
degrading.

Usage:  python3 diagnostics/plot_gain_agreement.py [outdir]
"""
import array
import glob
import os
import sys

import numpy as np
import ROOT

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from siwecal_eventbuilder.cli import resolve_gaudi_calib_files  # noqa: E402

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.environ.get("RUN", "TB2026CERN_run_000013")
CHUNKS = f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi/{RUN}/chunks"
TH = int(os.environ.get("TH", "230"))
OLD_SWITCH, NEW_SWITCH = 1200, 1900
N_ACQ = int(os.environ.get("N_ACQ", "4000"))
N_SCA, N_CHAN, N_CHIP = 15, 64, 16


def read_pedestals(path):
    """{(slab, chip, chan): [mean per SCA]} from a Pedestal_*.txt table."""
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 4:
            continue
        key = (int(f[0]), int(f[1]), int(f[2]))
        # 3 columns per SCA: mean, width, ... -- take every 3rd from index 3
        means = [float(x) for x in f[3::3]]
        out[key] = means
    return out


def read_mips(path):
    """{(slab, chip, chan): mpv} from a MIP_*.txt table (0 = masked)."""
    out = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        if len(f) < 4:
            continue
        out[(int(f[0]), int(f[1]), int(f[2]))] = float(f[3])
    return out


ped_hg_p, mip_hg_p, ped_lg_p, mip_lg_p = resolve_gaudi_calib_files(TH)
if not ped_lg_p:
    raise SystemExit(f"ERROR: no low-gain tables for th{TH} -- this plot needs both gains")
ped_hg, mip_hg = read_pedestals(ped_hg_p), read_mips(mip_hg_p)
ped_lg, mip_lg = read_pedestals(ped_lg_p), read_mips(mip_lg_p)
print(f"[calib] th{TH}: {len(mip_hg)} HG channels, {len(mip_lg)} LG channels")

chain = ROOT.TChain("siwecaldecoded")
for chunk in sorted(glob.glob(os.path.join(CHUNKS, "chunk_*.root")))[:20]:
    chain.AddFile(chunk)
if chain.GetEntries() == 0:
    raise SystemExit(f"ERROR: no decoded chunks under {CHUNKS}")

h = ROOT.TH2F("h", f"{RUN}: do the two gains agree?;adc_high (raw ADC);E(low gain) / E(high gain)",
              120, 200, 2600, 120, 0.0, 3.0)

n_hits = n_between = 0
for entry in range(min(N_ACQ, chain.GetEntries())):
    chain.GetEntry(entry)
    nsl = int(chain.n_slboards)
    sid = chain.GetLeaf("slboard_id")
    cid = chain.GetLeaf("chipid")
    hb = chain.GetLeaf("hitbit_high")
    ah = chain.GetLeaf("adc_high")
    al = chain.GetLeaf("adc_low")
    for slb in range(nsl):
        slab = int(sid.GetValue(slb))
        if slab < 0:
            continue
        for chip in range(N_CHIP):
            cc = int(cid.GetValue(slb * N_CHIP + chip))
            if cc < 0:
                continue
            for sca in range(N_SCA):
                base = ((slb * N_CHIP + chip) * N_SCA + sca) * N_CHAN
                for chn in range(N_CHAN):
                    if int(hb.GetValue(base + chn)) != 1:
                        continue
                    key = (slab, chip, chn)
                    mh, ml = mip_hg.get(key, 0.), mip_lg.get(key, 0.)
                    if mh <= 0 or ml <= 0:
                        continue  # masked in one gain: nothing to compare
                    ph = ped_hg.get(key, [])
                    pl = ped_lg.get(key, [])
                    if len(ph) <= sca or len(pl) <= sca:
                        continue
                    a_hi = float(ah.GetValue(base + chn))
                    a_lo = float(al.GetValue(base + chn))
                    e_hi = (a_hi - ph[sca]) / mh
                    e_lo = (a_lo - pl[sca]) / ml
                    if e_hi <= 0.5:      # below ~half a MIP: ratio is noise/noise
                        continue
                    n_hits += 1
                    if OLD_SWITCH <= a_hi < NEW_SWITCH:
                        n_between += 1
                    h.Fill(a_hi, e_lo / e_hi)

print(f"[hits] {n_hits:,} usable hits; {n_between:,} ({100 * n_between / max(n_hits, 1):.2f}%) "
      f"have adc_high in [{OLD_SWITCH}, {NEW_SWITCH}) -- these are the hits the old threshold "
      f"was handing to the low gain while the high gain was still linear")

c = ROOT.TCanvas("c", "gain agreement", 1250, 900)
ROOT.gPad.SetMargin(0.12, 0.14, 0.12, 0.08)
ROOT.gPad.SetLogz()
h.Draw("colz")

# Median profile of the ratio, per adc_high column.
prof = ROOT.TGraph()
q = array.array("d", [0.0])
prob = array.array("d", [0.5])
xs, ys = [], []
for ix in range(1, h.GetNbinsX() + 1):
    col = h.ProjectionY(f"c{ix}", ix, ix)
    if col.GetEntries() < 40:
        continue
    col.GetQuantiles(1, q, prob)
    x = h.GetXaxis().GetBinCenter(ix)
    prof.SetPoint(prof.GetN(), x, q[0])
    xs.append(x)
    ys.append(q[0])
prof.SetMarkerStyle(20)
prof.SetMarkerSize(0.7)
prof.SetMarkerColor(ROOT.kBlack)
prof.Draw("P same")

# Departure is measured from the PLATEAU, not from 1.
#
# The ratio does NOT sit at 1 where the high gain is linear -- it sits flat at
# ~0.76, i.e. the two gains' calibrations disagree by ~24% in scale. That is a
# real and separate problem (every saturation-recovered hit enters ~24% light),
# but it is a CALIBRATION offset, not a saturation effect: it shifts the whole
# curve up or down without moving the adc_high at which it starts to bend.
#
# So the switch point is where the curve LEAVES ITS OWN PLATEAU. An earlier
# version of this script measured |ratio - 1| instead and reported nonsense --
# the ratio is never 1, so the test fired in the very first bin it looked at.
xs, ys = np.array(xs), np.array(ys)
ref = (xs > 500) & (xs < 1200)          # high gain is unambiguously linear here
plateau = float(np.median(ys[ref]))
print(f"[plateau] E_low/E_high = {plateau:.3f} where the high gain is linear "
      f"(adc_high 500-1200). It should be 1.000 -- the two gains' calibrations "
      f"disagree by {100 * abs(1 - plateau):.0f}%, so every hit handed to the low "
      f"gain enters {100 * abs(1 - plateau):.0f}% light. Separate bug; see the README.")

keep = []
for frac, colour in ((0.05, ROOT.kOrange + 7), (0.20, ROOT.kMagenta + 2)):
    over = np.nonzero((ys > plateau * (1 + frac)) & (xs > 1200))[0]
    if over.size:
        x0 = xs[over[0]]
        ln = ROOT.TLine(x0, 0.0, x0, 3.0)
        ln.SetLineColor(colour)
        ln.SetLineWidth(2)
        ln.SetLineStyle(3)
        ln.Draw()
        keep.append(ln)
        print(f"[departure] the high gain has lost {100 * frac:.0f}% by adc_high = {x0:.0f}")

one = ROOT.TLine(200, plateau, 2600, plateau)
one.SetLineColor(ROOT.kGray + 2)
one.SetLineWidth(2)
one.SetLineStyle(2)
one.Draw()

leg = ROOT.TLegend(0.15, 0.66, 0.50, 0.90)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.026)
leg.AddEntry(prof, "median E_low / E_high", "p")
leg.AddEntry(one, "plateau (linear high gain)", "l")
for y, colour, style, label in ((OLD_SWITCH, ROOT.kRed + 1, 2, "previous switch (1200)"),
                                (NEW_SWITCH, ROOT.kGreen + 2, 1, "new switch (1900)")):
    ln = ROOT.TLine(y, 0.0, y, 3.0)
    ln.SetLineColor(colour)
    ln.SetLineWidth(3)
    ln.SetLineStyle(style)
    ln.Draw()
    leg.AddEntry(ln, label, "l")
    keep.append(ln)
leg.Draw()

# Per-threshold subfolder, always. A flat directory means the th220 run silently
# overwrites the th230 one -- which is exactly what happened once.
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)
out = os.path.join(OUT, f"gain_agreement_th{TH}.png")
c.SaveAs(out)
print(f"saved {out}")
