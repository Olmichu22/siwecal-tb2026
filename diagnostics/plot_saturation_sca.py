"""High gain vs low gain, hit by hit, ONE PANEL PER SCA.

The all-hits version (plot_saturation.py) shows a single band and tells you where
the high-gain preamp starts to bend. What it cannot tell you is whether that bend
is the same in every SCA column. The SKIROC's 15 SCA cells are separate storage
capacitors read through the same preamp, so a per-SCA difference in the HG/LG
relation would mean the saturation threshold (EVBLD_ADC_SATURATION, 1900) is not
one number but fifteen -- and SCA0 in particular is known to behave differently
(it carries the pedestal-tail structure that made a channel look double-peaked).

Each panel is the same adc_high:adc_low plane with the old (1200) and new (1900)
switch points drawn, plus the per-column ridge (the MEDIAN of adc_high in each
adc_low column, over all hits) so the bend is visible rather than inferred.

Usage:  python3 diagnostics/plot_saturation_sca.py [outdir]
"""
import array
import glob
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
RUN = os.environ.get("RUN", "TB2026CERN_run_000013")
TH = os.environ.get("TH", "230")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, f"th{TH}")
os.makedirs(OUT, exist_ok=True)
CHUNKS = f"/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi/{RUN}/chunks"
OLD_SWITCH, NEW_SWITCH = 1200, 1900
N_ACQ = 6000  # each panel gets ~1/15 of the hits; the band is dense enough at this depth
N_SCA = 15

t = ROOT.TChain("siwecaldecoded")
for chunk in sorted(glob.glob(os.path.join(CHUNKS, "chunk_*.root")))[:12]:
    t.AddFile(chunk)
if t.GetEntries() == 0:
    raise SystemExit(f"ERROR: no decoded chunks under {CHUNKS}")

c = ROOT.TCanvas("c", "saturation per SCA", 2000, 1500)
c.Divide(4, 4, 0.001, 0.001)

keep = []
for sca in range(N_SCA):
    pad = c.cd(sca + 1)
    pad.SetMargin(0.13, 0.13, 0.12, 0.08)
    pad.SetLogz()

    # The ADC arrays are [slab][chip][sca][channel]; Iteration$ walks the flattened
    # array, so the SCA index of an element is (Iteration$ / 64) % 15.
    sca_cut = f"(Iteration$/64)%{N_SCA}=={sca}"

    h = ROOT.TH2F(f"h_{sca}", f"SCA {sca};adc_low (raw ADC);adc_high (raw ADC)",
                  110, 150, 680, 130, 0, 2600)
    t.Draw(f"adc_high:adc_low>>h_{sca}", f"hitbit_high==1 && {sca_cut}", "goff", N_ACQ)
    h.Draw("colz")
    keep.append(h)

    # Ridge = MEDIAN of adc_high in each adc_low column, over ALL hits. No cut:
    # see plot_saturation.py -- cutting on adc_high to "isolate the band" makes the
    # ridge pin itself to the cut edge wherever the band has not risen yet, drawing
    # points where the data says otherwise.
    ridge = ROOT.TGraph()
    q = array.array("d", [0.0])
    prob = array.array("d", [0.5])
    for ix in range(1, h.GetNbinsX() + 1):
        col = h.ProjectionY(f"col_{sca}_{ix}", ix, ix)
        if col.GetEntries() < 20:
            continue
        col.GetQuantiles(1, q, prob)
        ridge.SetPoint(ridge.GetN(), h.GetXaxis().GetBinCenter(ix), q[0])
    ridge.SetMarkerStyle(20)
    ridge.SetMarkerSize(0.5)
    ridge.SetMarkerColor(ROOT.kBlack)
    ridge.Draw("P same")
    keep.append(ridge)

    for y, colour, style in ((OLD_SWITCH, ROOT.kRed + 1, 2), (NEW_SWITCH, ROOT.kGreen + 2, 1)):
        ln = ROOT.TLine(150, y, 680, y)
        ln.SetLineColor(colour)
        ln.SetLineWidth(2)
        ln.SetLineStyle(style)
        ln.Draw()
        keep.append(ln)

# Legend in the 16th, unused pad.
c.cd(16)
leg = ROOT.TLegend(0.05, 0.35, 0.95, 0.75)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.09)
for colour, style, label in ((ROOT.kRed + 1, 2, "previous switch (1200)"),
                             (ROOT.kGreen + 2, 1, "new switch (1900)")):
    ln = ROOT.TLine()
    ln.SetLineColor(colour)
    ln.SetLineWidth(3)
    ln.SetLineStyle(style)
    leg.AddEntry(ln, label, "l")
    keep.append(ln)
leg.Draw()

out = os.path.join(OUT, "saturation_hg_vs_lg_per_sca.png")
c.SaveAs(out)
print(f"saved {out}")
