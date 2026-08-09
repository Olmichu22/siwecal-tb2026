"""WHAT moves the pedestal in spill: the chip's own activity, or the detector's?

plot_pedestal_multipeak_spill.py establishes THAT the second pedestal population
is beam-induced (26.8x more multi-peak SCA0 cells in spill at matched
statistics). It cannot say what the beam does to make it happen. This does, and
it is the figure that carries the mechanism -- the result existed only as a
table before 2026-08-09.

METHOD. Only SCA slices with nhits<=1 enter, i.e. the same cut the pedestal fill
uses, so the chip was quiet in the very slice being measured. For each such
slice, every NON-HIT channel contributes

    dev = adc_low - (that cell's OUT-OF-SPILL mean)

-- how far the baseline has moved from where it sits with no beam at all. Each
measurement is then binned two ways:

    chipbusy = hits on the SAME chip, all its SCAs, in this acquisition
    detbusy  = total hits in the acquisition (the whole detector)

The two are correlated, so the discriminating question is which one the
deviation keeps tracking. If the satellites live at high chipbusy and are absent
at chipbusy<=1, the disturbance is chip-local. If they are already there at
chipbusy<=1 and grow with detbusy, it is global.

MEASURED ON run_000060, all 59 chunks, 18.1M measurements: |dev|>4 ADC goes
4.65% -> 53.22% across the chip-activity bins, while across detector-activity
bins it only creeps 24.7% -> 32.5% and saturates. Chip-local. The median also
moves DOWN, +0.14 to -2.44 ADC, which is what a common-mode sag looks like.

PANEL 4 ASKS A DIFFERENT QUESTION -- can any correction give the beam-off width
back? -- and the answer is mostly no. Two candidates, both tested out of sample:

  * an occupancy-dependent pedestal, i.e. one offset per chip-activity bin,
    fitted on EVEN acquisitions and applied to ODD ones: RMS -0.2%, tail -0.5%.
    It buys NOTHING, and in hindsight it could not: what grows with occupancy is
    the WIDTH, and an offset only moves the mean.
  * a per-slice common mode, the median deviation of the chip's OTHER quiet
    channels in the same slice (channels split even/odd so each half corrects
    the other, never itself): -4.0% overall, but that average hides the result.
    With the chip almost quiet (1 hit) it removes 46% and takes the RMS to
    1.66 ADC -- the out-of-spill pedestal width of these same cells is 1.70. At
    2-10 hits it removes nothing at all, so there the disturbance is not common
    to the 64 channels and a chip-wide median cannot see it.

THE ADC WINDOW IS LOAD-BEARING; see the comment at the cut. Without it the RMS
of every bin is garbage while the medians and tail fractions look fine, which is
the worst kind of bug -- the first version of this script had it.

CAVEAT: one run, low gain, th220. And note the reference means come from the
out-of-spill grid, so a cell only enters if it has >=200 out-of-spill entries --
which, out of spill, essentially means SCA0 (see the SCOPE note in section 7 of
diagnostics/pedestal_multipeak_tables.txt).

Usage:
  REF=<run_ns1_out.root> RUN=<dir with chunks/> [OCC_CUT=200] [MAXCHUNKS=0] \
      python3 diagnostics/plot_common_mode.py [outdir]
"""
import glob
import os
import re
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

_HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_HERE, "compare")
os.makedirs(OUTDIR, exist_ok=True)

REF = os.environ["REF"]
RUN = os.environ["RUN"]
OCC_CUT = float(os.environ.get("OCC_CUT", "200"))
MAXCHUNKS = int(os.environ.get("MAXCHUNKS", "0"))
MIN_REF_ENTRIES = int(os.environ.get("MIN_REF_ENTRIES", "200"))

NSL, NCHIP, NSCA, NCHN = 15, 16, 15, 64
# The pedestal grids are booked over 100.5-500.5 ADC; anything outside is not a
# baseline measurement and must not enter (see the note at the cut below).
ADC_LO, ADC_HI = 100.5, 500.5

# ped_low_s<slab>_c<chip>_ch<channel>_sca<sca>. Parse it with a regex: chained
# str.replace() does NOT work here, because replacing "_c" also eats the "_c"
# inside "_ch" and the channel token comes out as "h5".
CELL_RE = re.compile(r"^ped_low_s(\d+)_c(\d+)_ch(\d+)_sca(\d+)$")

ref = np.full((NSL, NCHIP, NSCA, NCHN), np.nan)
ref_rms = []
fref = ROOT.TFile.Open(REF)
if not fref or fref.IsZombie():
    raise SystemExit(f"cannot open reference grid: {REF}")
for k in fref.GetListOfKeys():
    m = CELL_RE.match(k.GetName())
    if not m:
        continue
    h = k.ReadObj()
    if h.GetEntries() < MIN_REF_ENTRIES:
        continue
    sl, ch, cn, sc = (int(x) for x in m.groups())
    ref[sl, ch, sc, cn] = h.GetMean()
    ref_rms.append(h.GetRMS())
fref.Close()
n_ref = int(np.isfinite(ref).sum())
print(f"reference cells with an out-of-spill mean (>= {MIN_REF_ENTRIES} entries): {n_ref:,}")
if n_ref == 0:
    raise SystemExit("no reference cells -- wrong grid?")
# The out-of-spill pedestal width of the very same cells. This is the target any
# in-spill correction is trying to get back to, so it belongs on the plot.
OUT_RMS = float(np.median(ref_rms)) if ref_rms else 0.0
print(f"their out-of-spill pedestal RMS: median {OUT_RMS:.2f} ADC")

paths = sorted(glob.glob(os.path.join(RUN, "chunks", "chunk_*.root")))
if MAXCHUNKS:
    paths = paths[:MAXCHUNKS]
if not paths:
    raise SystemExit(f"no chunks under {RUN}")
print(f"scanning {len(paths)} chunks of {os.path.basename(RUN.rstrip('/'))}")

bcid = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
bad = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
nh = np.zeros((NSL, NCHIP, NSCA), dtype=np.int32)
adc = np.zeros((NSL, NCHIP, NSCA, NCHN), dtype=np.int32)
hitb = np.zeros((NSL, NCHIP, NSCA, NCHN), dtype=np.int32)
nslb = np.zeros(1, dtype=np.int32)
acqn = np.zeros(1, dtype=np.int32)

# Half of the chip's quiet channels estimate a per-slice common mode for the
# other half. Both halves are kept, each carrying the estimate made from the
# channels it did NOT contribute to, so the test is out of sample by
# construction -- subtracting an estimate that includes the channel being
# corrected would reduce its own scatter no matter what the electronics does.
CM_MIN_CH = 8

dev_l, cb_l, db_l, par_l, cm_l = [], [], [], [], []
for ip, path in enumerate(paths):
    f = ROOT.TFile.Open(path)
    t = f.Get("siwecaldecoded")
    t.SetBranchStatus("*", 0)
    for b in ("acqNumber", "n_slboards", "bcid", "badbcid", "nhits", "adc_low", "hitbit_low"):
        t.SetBranchStatus(b, 1)
    t.SetBranchAddress("acqNumber", acqn)
    t.SetBranchAddress("n_slboards", nslb)
    t.SetBranchAddress("bcid", bcid)
    t.SetBranchAddress("badbcid", bad)
    t.SetBranchAddress("nhits", nh)
    t.SetBranchAddress("adc_low", adc)
    t.SetBranchAddress("hitbit_low", hitb)
    for i in range(t.GetEntries()):
        t.GetEntry(i)
        nsl = int(nslb[0])
        good = (bad[:nsl] == 0) & (bcid[:nsl] >= 0)
        detbusy = int(nh[:nsl][good].sum())
        if detbusy <= OCC_CUT:            # in-spill acquisitions only
            continue
        parity = int(acqn[0]) & 1
        chipbusy = np.where(good, nh[:nsl], 0).sum(axis=2)          # (nsl, 16)
        sl_i, ch_i, sc_i = np.nonzero(good & (nh[:nsl] <= 1))
        for sl, ch, sc in zip(sl_i, ch_i, sc_i):
            r = ref[sl, ch, sc]
            a = adc[sl, ch, sc]
            # The ADC window is NOT cosmetic. Without it the RMS is meaningless:
            # measured on all 59 chunks, the chipbusy==0 bin came out at RMS
            # 215 ADC and every detector-activity bin at ~38, while the medians
            # and tail fractions were unaffected -- i.e. a handful of invalid
            # readings per cell, all of them in slices where the chip recorded
            # nothing at all. This is the same window the pedestal histograms
            # are booked with (100.5-500.5), so restricting to it compares like
            # with like against the out-of-spill mean, which was computed inside it.
            m = np.isfinite(r) & (hitb[sl, ch, sc] == 0) & (a > ADC_LO) & (a < ADC_HI)
            k = int(m.sum())
            if not k:
                continue
            d = (a[m] - r[m]).astype(np.float32)
            dev_l.append(d)
            cb_l.append(np.full(k, chipbusy[sl, ch], dtype=np.int32))
            db_l.append(np.full(k, detbusy, dtype=np.int32))
            par_l.append(np.full(k, parity, dtype=np.int8))
            # per-slice common mode, estimated across the chip's OTHER channels
            chn = np.nonzero(m)[0]
            even = (chn % 2) == 0
            cm = np.full(k, np.nan, dtype=np.float32)
            if int(even.sum()) >= CM_MIN_CH and int((~even).sum()) >= CM_MIN_CH:
                cm[even] = np.median(d[~even])       # even channels corrected by odd
                cm[~even] = np.median(d[even])       # and the other way round
            cm_l.append(cm)
    f.Close()
    if (ip + 1) % 10 == 0:
        print(f"  {ip + 1}/{len(paths)} chunks, {sum(len(x) for x in dev_l):,} measurements")

dev = np.concatenate(dev_l)
cb = np.concatenate(cb_l)
db = np.concatenate(db_l)
par = np.concatenate(par_l)
cm = np.concatenate(cm_l)
del dev_l, cb_l, db_l, par_l, cm_l
print(f"non-hit channel measurements in quiet slices, in spill: {len(dev):,}\n")

CB_EDGES = [0, 1, 2, 5, 10, 20, 50, 1e9]
DB_EDGES = [OCC_CUT, 500, 800, 1100, 1400, 1e9]


def profile(var, edges):
    """(label, n, median, rms, frac|dev|>4) per bin, skipping starved bins."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (var >= lo) & (var < hi)
        n = int(m.sum())
        if n < 50:
            continue
        d = dev[m]
        hl = "inf" if hi > 1e8 else f"{int(hi)}"
        rows.append((f"[{int(lo)},{hl})", n, float(np.median(d)), float(d.std()),
                     float((np.abs(d) > 4).mean())))
    return rows


cb_rows = profile(cb, CB_EDGES)
db_rows = profile(db, DB_EDGES)
# The chipbusy==0 bin is kept in the printed table but left OFF the curves. A
# slice recorded on a chip that saw no hits at all in the acquisition is not a
# measurement of what the beam does to that chip, and it behaves unlike every
# other bin (median +0.59, tail 23.5%), so plotting it turns a monotonic trend
# into a spurious V.
cb_curve = [r for r in cb_rows if not r[0].startswith("[0,")]
for name, rows in (("CHIP activity (same chip, all SCAs, this acquisition)", cb_rows),
                   ("DETECTOR activity (whole acquisition)", db_rows)):
    print(f"baseline deviation from the out-of-spill mean, binned by {name}:")
    print("   bin              cases       median    RMS    |dev|>4 ADC")
    for lab, n, med, rms, fr in rows:
        print(f"   {lab:<14} {n:>11,}   {med:+6.2f}  {rms:5.2f}    {100 * fr:5.2f}%")
    print()

# --- panel 4 input: an occupancy-dependent pedestal, trained out of sample.
# Offset per chipbusy value (capped), fitted on EVEN acquisitions only and
# applied to ODD ones. Anything fitted and evaluated on the same events would
# improve by construction and would prove nothing.
CAP = 40
cbc = np.minimum(cb, CAP)
train, test = par == 0, par == 1
off = np.zeros(CAP + 1)
for v in range(CAP + 1):
    m = train & (cbc == v)
    if m.sum() >= 200:
        off[v] = np.median(dev[m])
raw = dev[test]
cor = raw - off[cbc[test]]
gain_rms = 100 * (1 - cor.std() / raw.std())
gain_tail = 100 * (1 - (np.abs(cor) > 4).mean() / max((np.abs(raw) > 4).mean(), 1e-12))
print(f"occupancy-dependent pedestal, trained on even acquisitions and tested on "
      f"odd ({test.sum():,} measurements):")
print(f"   raw        RMS {raw.std():.3f} ADC   |dev|>4 {100 * (np.abs(raw) > 4).mean():.2f}%")
print(f"   corrected  RMS {cor.std():.3f} ADC   |dev|>4 {100 * (np.abs(cor) > 4).mean():.2f}%")
print(f"   => RMS -{gain_rms:.1f}%, tail -{gain_tail:.1f}%\n")

# --- the other candidate correction, and the one the chip-local result actually
# points at: a per-slice common mode measured on the chip's own quiet channels.
# An offset per occupancy bin can only move the mean; if what grows with
# occupancy is the WIDTH, only something measured slice by slice can remove it.
have_cm = np.isfinite(cm)
raw_cm = dev[have_cm]
cor_cm = raw_cm - cm[have_cm]
cm_rms = 100 * (1 - cor_cm.std() / raw_cm.std())
cm_tail = 100 * (1 - (np.abs(cor_cm) > 4).mean() / max((np.abs(raw_cm) > 4).mean(), 1e-12))
print(f"per-slice common mode from the chip's OTHER quiet channels "
      f"({have_cm.sum():,} measurements with an estimate, "
      f"{100 * have_cm.mean():.1f}% of the sample):")
print(f"   raw        RMS {raw_cm.std():.3f} ADC   "
      f"|dev|>4 {100 * (np.abs(raw_cm) > 4).mean():.2f}%")
print(f"   corrected  RMS {cor_cm.std():.3f} ADC   "
      f"|dev|>4 {100 * (np.abs(cor_cm) > 4).mean():.2f}%")
print(f"   => RMS -{cm_rms:.1f}%, tail -{cm_tail:.1f}%")
# How much of the disturbance is coherent should GROW with chip activity if the
# chip-local reading is right.
print("   by chip activity:")
for lo, hi in zip(CB_EDGES[:-1], CB_EDGES[1:]):
    msk = have_cm & (cb >= lo) & (cb < hi)
    if msk.sum() < 200:
        continue
    a0 = dev[msk]
    a1 = a0 - cm[msk]
    hl = "inf" if hi > 1e8 else f"{int(hi)}"
    print(f"     [{int(lo)},{hl})".ljust(16) +
          f"{int(msk.sum()):>11,}   RMS {a0.std():5.2f} -> {a1.std():5.2f}   "
          f"({100 * (1 - a1.std() / max(a0.std(), 1e-9)):5.1f}% removed)")
print()

# ------------------------------------------------------------------ figure
C_QUIET, C_BUSY, C_CHIP, C_DET = ROOT.kAzure + 2, ROOT.kRed + 1, ROOT.kRed + 1, ROOT.kAzure + 2
c = ROOT.TCanvas("c", "common mode", 1800, 1300)
c.Divide(2, 2, 0.010, 0.014)
keep = []


def style_pad():
    ROOT.gPad.SetMargin(0.13, 0.05, 0.13, 0.10)
    ROOT.gPad.SetGridy()


def band(lo, hi):
    return (cb >= lo) & (cb < hi)


# --- 1. the deviation itself, quiet chip vs busy chip
c.cd(1)
style_pad()
ROOT.gPad.SetLogy()
h_q = ROOT.TH1F("h_q", "Baseline deviation of quiet, non-hit channels;"
                       "adc_low - out-of-spill mean of that cell (ADC);"
                       "fraction of measurements", 160, -20, 20)
h_b = ROOT.TH1F("h_b", "", 160, -20, 20)
for h, m, col in ((h_q, band(0, 2), C_QUIET), (h_b, band(20, 1e9), C_BUSY)):
    d = dev[m]
    for v in d:
        h.Fill(v)
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral())
    h.SetLineColor(col)
    h.SetLineWidth(3)
h_q.SetMaximum(max(h_q.GetMaximum(), h_b.GetMaximum()) * 3)
h_q.SetMinimum(1e-6)
h_q.Draw("hist")
h_b.Draw("hist same")
leg = ROOT.TLegend(0.58, 0.72, 0.94, 0.88)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(h_q, f"chip quiet (<2 hits), {int(band(0, 2).sum()):,}", "l")
leg.AddEntry(h_b, f"chip busy (#geq20 hits), {int(band(20, 1e9).sum()):,}", "l")
leg.Draw()
keep += [h_q, h_b, leg]

# --- 2. THE panel: tail fraction against chip activity and detector activity,
# both on one x axis of "hits in the acquisition" so they can be compared.
c.cd(2)
style_pad()
ROOT.gPad.SetLogx()


def graph(rows, edges, col, marker):
    g = ROOT.TGraph()
    lo_edges = [e for e, hi in zip(edges[:-1], edges[1:])]
    j = 0
    for (lab, n, med, rms, fr) in rows:
        lo = float(lab.strip("[)").split(",")[0])
        hi = lab.strip("[)").split(",")[1]
        hi = lo * 3 if hi == "inf" else float(hi)
        g.SetPoint(j, max(np.sqrt(max(lo, 0.5) * hi), 0.5), 100 * fr)
        j += 1
    g.SetLineColor(col)
    g.SetMarkerColor(col)
    g.SetLineWidth(3)
    g.SetMarkerStyle(marker)
    g.SetMarkerSize(1.8)
    return g


g_cb = graph(cb_curve, CB_EDGES, C_CHIP, 20)
g_db = graph(db_rows, DB_EDGES, C_DET, 21)
fr2 = ROOT.TH1F("fr2", "What the deviation tracks;"
                       "hits in the acquisition (bin centre, log);"
                       "measurements with |deviation| > 4 ADC (%)", 100, 0.4, 4000)
fr2.SetMinimum(0)
fr2.SetMaximum(60)
fr2.Draw()
g_cb.Draw("plsame")
g_db.Draw("plsame")
leg2 = ROOT.TLegend(0.16, 0.70, 0.62, 0.87)
leg2.SetBorderSize(0)
leg2.SetFillStyle(0)
leg2.AddEntry(g_cb, "binned by the CHIP's own hits", "pl")
leg2.AddEntry(g_db, "binned by the whole DETECTOR's hits", "pl")
leg2.Draw()
lat2 = ROOT.TLatex()
lat2.SetNDC()
lat2.SetTextSize(0.033)
lat2.DrawLatex(0.40, 0.24, "chip-local: the detector curve saturates near 32%,")
lat2.DrawLatex(0.40, 0.19, "the chip curve keeps climbing to 53%")
keep += [fr2, g_cb, g_db, leg2, lat2]

# --- 3. the median: a sag, not a symmetric widening
c.cd(3)
style_pad()
hm = ROOT.TH1F("hm", "The baseline sags with the chip's own activity;"
                     "hits on the chip in the acquisition;"
                     "median deviation (ADC)", len(cb_curve), 0, len(cb_curve))
for i, (lab, n, med, rms, fr) in enumerate(cb_curve):
    hm.SetBinContent(i + 1, med)
    hm.GetXaxis().SetBinLabel(i + 1, lab)
hm.GetXaxis().SetLabelSize(0.048)
hm.SetLineColor(C_CHIP)
hm.SetLineWidth(3)
hm.SetFillColorAlpha(C_CHIP, 0.30)
lo_m = min(r[2] for r in cb_curve)
hm.SetMinimum(min(lo_m * 1.35, -0.5))
hm.SetMaximum(max(0.6, max(r[2] for r in cb_curve) * 1.4))
hm.Draw("hist")
z = ROOT.TLine(0, 0, len(cb_curve), 0)
z.SetLineStyle(2)
z.Draw()
keep += [hm, z]

# --- 4. can either correction give the beam-off width back? Per chip-activity
# bin, because the overall figure hides the only regime where one of them works.
c.cd(4)
style_pad()
nb4 = len(cb_curve)
h_r = ROOT.TH1F("h_r", "Can a correction recover the beam-off width?;"
                       "hits on the chip in the acquisition;"
                       "RMS of the deviation (ADC)", nb4, 0, nb4)
h_m = ROOT.TH1F("h_m", "", nb4, 0, nb4)
for i, (lab, n, med, rms, fr) in enumerate(cb_curve):
    lo = float(lab.strip("[)").split(",")[0])
    hi_s = lab.strip("[)").split(",")[1]
    hi = 1e9 if hi_s == "inf" else float(hi_s)
    msk = have_cm & (cb >= lo) & (cb < hi)
    h_r.SetBinContent(i + 1, rms)
    if msk.sum() >= 200:
        h_m.SetBinContent(i + 1, float((dev[msk] - cm[msk]).std()))
    h_r.GetXaxis().SetBinLabel(i + 1, lab)
h_r.GetXaxis().SetLabelSize(0.048)
h_r.SetLineColor(C_BUSY)
h_r.SetLineWidth(3)
h_m.SetLineColor(ROOT.kGreen + 2)
h_m.SetLineWidth(3)
h_r.SetMinimum(0)
h_r.SetMaximum(max(r[3] for r in cb_curve) * 1.9)
h_r.Draw("hist")
h_m.Draw("hist same")
tgt = ROOT.TLine(0, OUT_RMS, nb4, OUT_RMS)
tgt.SetLineStyle(2)
tgt.SetLineWidth(2)
tgt.Draw()
leg4 = ROOT.TLegend(0.16, 0.74, 0.62, 0.89)
leg4.SetBorderSize(0)
leg4.SetFillStyle(0)
leg4.AddEntry(h_r, "as calibrated", "l")
leg4.AddEntry(h_m, "minus per-slice common mode", "l")
leg4.AddEntry(tgt, f"out-of-spill pedestal RMS ({OUT_RMS:.2f} ADC)", "l")
leg4.Draw()
lat4 = ROOT.TLatex()
lat4.SetNDC()
lat4.SetTextSize(0.030)
# Keep every annotation ABOVE the curves: the one bin that matters here is the
# leftmost, where the green drops to the beam-off line, and text at the bottom
# of the pad lands exactly on it.
lat4.DrawLatex(0.63, 0.86, "common mode: median of the")
lat4.DrawLatex(0.63, 0.82, "chip's OTHER quiet channels,")
lat4.DrawLatex(0.63, 0.78, f"same slice -- overall RMS -{cm_rms:.1f}%")
lat4.DrawLatex(0.63, 0.72, "occupancy-dependent OFFSET")
lat4.DrawLatex(0.63, 0.68, f"instead: RMS -{gain_rms:.1f}% -- nothing")
keep += [h_r, h_m, tgt, leg4, lat4]

out = os.path.join(OUTDIR, "common_mode.png")
c.SaveAs(out)
print(f"saved {out}")
