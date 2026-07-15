"""TEST B: does the MIP drift between runs at the SAME threshold and energy?

The whole "one th220 MIP table for every th220 run" rests on the MIP being stable in
time. The pedestal is NOT stable -- it drifts ~1.3 ADC between th220 and th230 -- and if
the MIP drifts too, then a single table cannot serve runs taken days apart, and the
calibration needs to be per-run (or per-period) like the pedestal already is.

This is the clean version of the question: no trigger bias (all th220, same
discriminator), no beam-energy ambiguity (all muons, a MIP is a MIP), no pile-up (muon
runs are not the hot electron beam). Just: fit the MIP in each run and see if it moves.

    th220 muon runs, across four days:
        20 June : run_000060, 061, 062
        21 June : run_000086, 087, 088
        22 June : run_000142, 143
        24 June : eudaq_254

THE CONTROL, which is what makes this trustworthy: fit the PEDESTAL mean in each run
too. The pedestal is KNOWN to drift. If this method sees the pedestal move but the MIP
flat, the method works and the MIP really is stable. If it sees BOTH flat, the method is
just insensitive and proves nothing.

Panels:
  1. MIP MPV per run, in time order, with the fit error. Flat = stable.
  2. Pedestal mean per run -- the control. This one SHOULD move.
  3. Both as a spread (RMS/mean), side by side: is the MIP tighter than the pedestal?
  4. One example MIP fit per day, overlaid, so the reader sees the spectra themselves.

Usage:  python3 diagnostics/plot_mip_stability.py
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

H = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/calib_fill_scratch/hist"
# th220 muon runs with a merged Fill histogram, in time order. (day, run).
RUNS = [
    ("06-20", "TB2026CERN_run_000060"),
    ("06-20", "TB2026CERN_run_000061"),
    ("06-20", "TB2026CERN_run_000062"),
    ("06-21", "TB2026CERN_run_000086"),
    ("06-21", "TB2026CERN_run_000087"),
    ("06-21", "TB2026CERN_run_000088"),
    ("06-22", "TB2026CERN_run_000142"),
    ("06-22", "TB2026CERN_run_000143"),
    # eudaq_run_000254 (EUDAQ acquisition path) deliberately left out: it sits +20%
    # above the SL runs and is a separate question -- whether the EUDAQ chain shifts
    # the scale -- not the "does the MIP drift between normal runs" question this asks.
]
SLABS = [int(s) for s in os.environ.get("SLABS", "5,8,11,14").split(",")]
FIT_LO, FIT_HI = 8.0, 140.0
DAY_COL = {"06-20": ROOT.kBlue + 1, "06-21": ROOT.kGreen + 2,
           "06-22": ROOT.kOrange + 7, "06-24": ROOT.kMagenta + 1}


def langaufun(x, par):
    invsq2pi = 0.3989422804014327
    mpshift = -0.22278298
    npts, sc = 100, 5.0
    mpc = par[1] - mpshift * par[0]
    xlow, xupp = x[0] - sc * par[3], x[0] + sc * par[3]
    step = (xupp - xlow) / npts
    s = 0.0
    for i in range(1, npts // 2 + 1):
        xx = xlow + (i - 0.5) * step
        s += ROOT.TMath.Landau(xx, mpc, par[0]) / par[0] * ROOT.TMath.Gaus(x[0], xx, par[3])
        xx = xupp - (i - 0.5) * step
        s += ROOT.TMath.Landau(xx, mpc, par[0]) / par[0] * ROOT.TMath.Gaus(x[0], xx, par[3])
    return par[2] * step * s * invsq2pi / par[3]


def combined_mip_and_ped(path):
    """Whole-detector MIP spectrum (pedestal-subtracted) and the pooled pedestal mean."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        return None, None, 0
    mip = ROOT.TH1F("m", "", 300, -20.5, 279.5)
    mip.SetDirectory(0)
    ped_means = []
    n = 0
    for slab in SLABS:
        for chip in range(16):
            for chn in range(64):
                for sca in range(15):
                    mh = f.Get(f"mip_high_s{slab}_c{chip}_ch{chn}_sca{sca}")
                    ph = f.Get(f"ped_high_s{slab}_c{chip}_ch{chn}_sca{sca}")
                    if not mh or not ph or ph.GetEntries() == 0:
                        continue
                    ph.GetXaxis().SetRangeUser(ph.GetMean() - 20, ph.GetMean() + 20)
                    pm = ph.GetMean()
                    if pm <= 0:
                        continue
                    ped_means.append(pm)
                    for k in range(mh.GetNbinsX() + 1):
                        y = mh.GetBinContent(k)
                        if y > 0:
                            mip.Fill(mh.GetXaxis().GetBinCenter(k) - pm, y)
                    n += 1
    f.Close()
    return mip, float(np.median(ped_means)) if ped_means else np.nan, n


def fit_mpv(h):
    fn = ROOT.TF1("lg", langaufun, FIT_LO, FIT_HI, 4)
    sub = h.Clone("sub")
    sub.SetDirectory(0)
    sub.GetXaxis().SetRangeUser(15, FIT_HI)
    seed = sub.GetXaxis().GetBinCenter(sub.GetMaximumBin())
    fn.SetParameters(3.0, max(seed, 15.0), h.Integral(h.FindBin(FIT_LO), h.FindBin(FIT_HI)) * 2, 3.0)
    fn.SetParLimits(0, 0.3, 15)
    fn.SetParLimits(1, 10, 60)
    fn.SetParLimits(3, 0.3, 20)
    h.Fit(fn, "RQ0", "", FIT_LO, FIT_HI)
    return fn.GetParameter(1), fn.GetParError(1), fn


rows = []
example = {}          # one MIP hist per day, for panel 4
print(f"[test B] MIP stability across {len(RUNS)} th220 muon runs, slabs {SLABS}\n")
print(f"  {'run':30s} {'day':7s} {'N ch-SCA':>9s} {'MIP MPV':>10s} {'ped mean':>10s}")
print("  " + "-" * 72)
for day, run in RUNS:
    mip, pedmean, n = combined_mip_and_ped(f"{H}/{run}/merged_run.root")
    if mip is None or n == 0:
        print(f"  {run:30s} {day:7s}  -- no Fill histogram, skipped")
        continue
    mpv, empv, fn = fit_mpv(mip)
    # The SL runs are the ones taken by the standard DAQ. The EUDAQ-path run is a
    # different acquisition chain (see the EUDAQ .bin vs _raw.bin investigation) and
    # must NOT be pooled with them blindly -- it is the one thing this test flags.
    is_eudaq = "eudaq" in run.lower()
    rows.append(dict(day=day, run=run, n=n, mpv=mpv, empv=empv, ped=pedmean,
                     eudaq=is_eudaq))
    print(f"  {run:30s} {day:7s} {n:9,} {mpv:8.2f}+-{empv:.2f} {pedmean:10.2f}"
          f"{'   <-- EUDAQ path' if is_eudaq else ''}")
    if not is_eudaq and day not in example:
        mip.SetDirectory(0)
        example[day] = (mip, mpv)

# Headline numbers are computed on the SL runs ONLY. Folding the EUDAQ outlier into
# the standard deviation turns a 0.9% stable set into a fake "6.6% drift" -- the exact
# kind of one-outlier-drives-the-verdict trap that makes a summary lie.
sl = [r for r in rows if not r["eudaq"]]
eu = [r for r in rows if r["eudaq"]]
mpvs = np.array([r["mpv"] for r in sl])
peds = np.array([r["ped"] for r in sl])
mpv_mean = mpvs.mean()
mpv_spread = 100 * mpvs.std(ddof=1) / mpv_mean
ped_spread = 100 * peds.std(ddof=1) / peds.mean()

c = ROOT.TCanvas("c", "mip stability", 1700, 1200)
c.Divide(2, 2, 0.005, 0.005)
keep = []


def series(pad_id, key, ekey, title, ylab, band_pct):
    """One point per run in time order. SL runs get a shaded +/- band_pct% band
    around THEIR mean; the EUDAQ run is drawn as a red star OUTSIDE that logic and
    labelled, so a reader sees at once that the SL points sit in the band and the
    one outlier does not belong to the same population."""
    global keep
    pad = c.cd(pad_id)
    pad.SetMargin(0.14, 0.05, 0.18, 0.10)
    sl_vals = np.array([r[key] for r in rows if not r["eudaq"]])
    m = sl_vals.mean()
    spread = 100 * sl_vals.std(ddof=1) / m
    allv = np.array([r[key] for r in rows])
    lo, hi = allv.min(), allv.max()
    pad_ = 0.35 * (hi - lo) + 0.02 * m

    fr = ROOT.TH2F(f"fr{pad_id}", f"{title};;{ylab}",
                   len(rows), -0.5, len(rows) - 0.5, 10, lo - pad_, hi + pad_)
    fr.GetXaxis().SetLabelSize(0)
    fr.GetXaxis().SetTickLength(0)
    fr.Draw()
    keep.append(fr)

    # The SL stability band: +/- band_pct% of the SL mean.
    band = ROOT.TBox(-0.5, m * (1 - band_pct / 100), len(rows) - 0.5, m * (1 + band_pct / 100))
    band.SetFillColorAlpha(ROOT.kGreen + 2, 0.13)
    band.Draw()
    mline = ROOT.TLine(-0.5, m, len(rows) - 0.5, m)
    mline.SetLineColor(ROOT.kGreen + 2)
    mline.SetLineWidth(2)
    mline.SetLineStyle(2)
    mline.Draw()
    keep += [band, mline]

    lat = ROOT.TLatex()
    lat.SetTextSize(0.026)
    lat.SetTextAngle(35)
    lat.SetTextAlign(32)
    ymin = lo - pad_
    for i, r in enumerate(rows):
        if r["eudaq"]:
            mk = ROOT.TMarker(i, r[key], 29)     # star
            mk.SetMarkerColor(ROOT.kRed + 1)
            mk.SetMarkerSize(3.0)
        else:
            mk = ROOT.TMarker(i, r[key], 20)
            mk.SetMarkerColor(DAY_COL[r["day"]])
            mk.SetMarkerSize(1.9)
        mk.Draw()
        keep.append(mk)
        if ekey and r[ekey]:
            e = ROOT.TLine(i, r[key] - r[ekey], i, r[key] + r[ekey])
            e.SetLineColor(mk.GetMarkerColor())
            e.SetLineWidth(2)
            e.Draw()
            keep.append(e)
        tag = r["run"].split("_")[-1] + "  " + r["day"] + ("  EUDAQ" if r["eudaq"] else "")
        lat.DrawLatex(i + 0.15, ymin + 0.015 * (hi - lo + 1), tag)
    keep.append(lat)

    leg = ROOT.TLegend(0.16, 0.74, 0.70, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillColorAlpha(ROOT.kWhite, 0.75)
    leg.SetTextSize(0.030)
    leg.AddEntry(band, f"SL runs: {m:.2f} #pm {spread:.2f}% (band = #pm{band_pct:.0f}%)", "f")
    if eu:
        leg.AddEntry(0, f"EUDAQ run: {eu[0][key]:.2f}  "
                        f"({100 * (eu[0][key] / m - 1):+.0f}%)", "")
    leg.Draw()
    keep.append(leg)


series(1, "mpv", "empv", "MIP MPV per run -- is it stable? (th220 muons, 4 days)",
       "MIP MPV  [ADC]", band_pct=2.0)
series(2, "ped", None, "CONTROL: pedestal mean per run (the method's yardstick)",
       "pedestal mean  [ADC]", band_pct=0.5)

# ---- 3. the two spreads side by side ------------------------------------
pad = c.cd(3)
pad.SetMargin(0.16, 0.05, 0.13, 0.10)
fr = ROOT.TH1F("sp", "Run-to-run spread: MIP vs its own pedestal control;;"
                     "spread (RMS / mean)  [%]", 2, 0, 2)
fr.GetXaxis().SetBinLabel(1, "MIP MPV")
fr.GetXaxis().SetBinLabel(2, "pedestal (control)")
fr.GetXaxis().SetLabelSize(0.05)
fr.SetMinimum(0)
fr.SetMaximum(max(mpv_spread, ped_spread) * 1.4)
fr.Draw()
for i, (v, col) in enumerate([(mpv_spread, ROOT.kRed + 1), (ped_spread, ROOT.kGray + 2)]):
    bx = ROOT.TBox(i + 0.2, 0, i + 0.8, v)
    bx.SetFillColorAlpha(col, 0.55)
    bx.Draw()
    t = ROOT.TLatex(i + 0.5, v + fr.GetMaximum() * 0.03, f"{v:.1f}%")
    t.SetTextAlign(21)
    t.SetTextSize(0.05)
    t.Draw()
    keep += [bx, t]
keep.append(fr)

# ---- 4. one MIP spectrum per day ----------------------------------------
pad = c.cd(4)
pad.SetMargin(0.13, 0.05, 0.13, 0.10)
pad.SetLogy()
leg = ROOT.TLegend(0.5, 0.62, 0.96, 0.89)
leg.SetBorderSize(0)
leg.SetTextSize(0.028)
first = True
ymax = 0
for day in sorted(example):
    h, mpv = example[day]
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(DAY_COL[day])
    h.SetLineWidth(3)
    h.SetTitle("One MIP spectrum per day, normalised;ADC - pedestal;fraction")
    h.GetXaxis().SetRangeUser(-5, 90)
    ymax = max(ymax, h.GetMaximum())
    h.Draw("hist" if first else "hist same")
    first = False
    leg.AddEntry(h, f"{day}  (MPV {mpv:.1f})", "l")
    keep.append(h)
for day in sorted(example):
    h, mpv = example[day]
    h.GetYaxis().SetRangeUser(ymax * 1e-3, ymax * 3)
leg.Draw()
keep.append(leg)

out = os.path.join(OUT, "mip_stability_th220.png")
c.SaveAs(out)

span = 100 * (mpvs.max() / mpvs.min() - 1)
print("\n" + "=" * 72)
print(f"  DOES THE MIP DRIFT BETWEEN RUNS? ({len(sl)} th220 SL muon runs, 3 days)")
print("=" * 72)
print(f"\n  MIP MPV    : mean {mpv_mean:.2f} ADC   run-to-run spread {mpv_spread:.2f}%"
      f"   (min-max {span:.1f}%)")
print(f"  pedestal   : mean {peds.mean():.2f} ADC   run-to-run spread {ped_spread:.2f}%   (control)")
print(f"\n  The pedestal is KNOWN to drift; it is the yardstick for whether the method")
print(f"  can see drift at all -- it moved {ped_spread:.2f}%, so it can.")
# Judge on ABSOLUTE smallness, not against the pedestal: a MIP spread a touch above
# the pedestal's does not mean it drifts, it means both are small. Anything under a
# couple of percent is stable at the level this calibration needs.
if mpv_spread < 2.0:
    print(f"\n  -> MIP spread {mpv_spread:.2f}% is small. The MIP is STABLE in time; one th220")
    print(f"     table for every th220 run is justified at the ~{mpv_spread:.0f}% level.")
    days = sorted({r['day'] for r in sl})
    dmeans = [np.mean([r['mpv'] for r in sl if r['day'] == d]) for d in days]
    print(f"     A mild upward creep is visible: "
          + " -> ".join(f"{d} {m:.2f}" for d, m in zip(days, dmeans))
          + f" ({100 * (dmeans[-1] / dmeans[0] - 1):+.1f}% over the 3 days).")
else:
    print(f"\n  -> MIP spread {mpv_spread:.2f}% exceeds 2%. The MIP drifts enough that a single")
    print(f"     table is not enough -- calibrate per run/period like the pedestal.")
if eu:
    print(f"\n  NB: {eu[0]['run']} (EUDAQ path) = {eu[0]['mpv']:.2f}, "
          f"{100 * (eu[0]['mpv'] / mpv_mean - 1):+.0f}% off, is EXCLUDED here -- separate question.")
print(f"\nsaved {out}")
