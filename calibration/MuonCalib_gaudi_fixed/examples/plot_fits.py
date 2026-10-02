"""One example pedestal fit and one MIP fit, with the calibrator's own fit code (PedestalMipCalib.h),
on the Fill grid of muon run eudaq 152 (th210, fixed chunks)."""
import sys
import numpy as np
import ROOT
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = "/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026"
ROOT.gInterpreter.AddIncludePath(f"{REPO}/gaudi_source/include")
ROOT.gInterpreter.Declare('#include "k4SiWEcalReco/PedestalMipCalib.h"')
k = ROOT.k4siwecal
# fitLangau keeps its TF1 local; refit-free way to get the full curve: run fitLangau, then rebuild the same TF1
# from the same start values and limits and fit again (identical procedure, so identical parameters)
ROOT.gInterpreter.Declare('''
std::vector<double> langauParams(TH1F& h) {
  TH1F* c = (TH1F*)h.Clone("mip_clone"); c->SetDirectory(nullptr);
  k4siwecal::fitLangau(*c, true);
  double pllo[4] = {0.1, 14., 1.0, 0.5}, plhi[4] = {20.0, 220.0, 1.0e8, 10.0}, fr[2] = {0., 100.};
  TH1F& hist = h;
  hist.GetXaxis()->SetRangeUser(fr[0], fr[1]);
  TF1 fitlandau("fl2", "landau", fr[0], fr[1]);
  hist.Fit(&fitlandau, "QRE0");
  double sv[4];
  sv[0] = std::min(std::max(fitlandau.GetParameter(2), pllo[0]), plhi[0]);
  sv[1] = std::min(std::max(fitlandau.GetParameter(1), pllo[1]), plhi[1]);
  sv[2] = hist.Integral("width") * 1.2;
  sv[3] = std::min(std::max(hist.GetRMS() * 0.1, pllo[3]), plhi[3]);
  TF1 langau("lg2", k4siwecal::langaufun, fr[0], fr[1], 4);
  langau.SetParameters(sv);
  for (int i = 0; i < 4; ++i) langau.SetParLimits(i, pllo[i], plhi[i]);
  hist.Fit(&langau, "RBQM0");
  return {langau.GetParameter(0), langau.GetParameter(1), langau.GetParameter(2), langau.GetParameter(3)};
}
double langauEval(double x, std::vector<double> p) { return k4siwecal::langaufun(&x, p.data()); }
''')

fin, outdir, slab, chip, chan, sca = sys.argv[1], sys.argv[2], *map(int, sys.argv[3:7])
f = ROOT.TFile.Open(fin)
names = set(x.GetName() for x in f.GetListOfKeys())

def get(prefix, s=sca):
    for n in (f"{prefix}_s{slab}_c{chip}_ch{chan}_sca{s}", f"{prefix}_s{slab}_c{chip}_ch{chan}"):
        if n in names:
            h = f.Get(n); h.SetDirectory(0); return h
    raise SystemExit(f"no {prefix} histogram for {slab} {chip} {chan} {s}; ")

def arrays(h):
    n = h.GetNbinsX()
    x = np.array([h.GetBinCenter(i) for i in range(1, n + 1)]); y = np.array([h.GetBinContent(i) for i in range(1, n + 1)])
    return x, y, h.GetBinWidth(1)

plt.rcParams.update({"font.size": 13, "axes.spines.top": False, "axes.spines.right": False})

# pedestal: dual-window Gaussian, narrower of the two kept
hp = get("ped_high")
r = k.fitPedestalSca(hp)
x, y, w = arrays(hp)
mu, sg = r.mean, r.width
sel = (x > mu - 13 * sg) & (x < mu + 7 * sg)
xx = np.linspace(mu - 4 * sg, mu + 4 * sg, 400)
amp = y[sel].sum() * w / (np.sqrt(2 * np.pi) * sg)
fig, ax = plt.subplots(figsize=(6.4, 4.6))
ax.step(x[sel], y[sel], where="mid", color="#47505E", lw=1.4, label="data (4 muon runs)")
ax.plot(xx, amp * np.exp(-0.5 * ((xx - mu) / sg) ** 2), color="#2B6CB0", lw=2.2, label="Gaussian fit")
ax.set_xlabel("high-gain ADC"); ax.set_ylabel("entries / ADC")
ax.set_title(f"Pedestal · slab {slab}, chip {chip}, channel {chan}, SCA {sca}", fontsize=13)
ax.text(0.03, 0.74, f"{int(hp.GetEntries())} entries\nmean = {mu:.2f} ± {r.error:.2f} ADC\nwidth = {sg:.3f} ± {r.widthError:.3f} ADC",
        transform=ax.transAxes, ha="left", va="top", fontsize=12)
ax.legend(loc="upper left", frameon=False, fontsize=11)
fig.tight_layout(); fig.savefig(f"{outdir}/pedestal_fit.png", dpi=160)

# MIP: Landau (x) Gaussian on pedestal-subtracted ADC
# the MIP fit is per channel: the 15 SCA histograms combined, each minus its own truncated-mean pedestal
# (PedestalMipCalibrator::writeMipTable, quirks included)
hm = ROOT.TH1F("hmips", "", 900, -100.5, 799.5); hm.SetDirectory(0)
for s_ in range(15):
    n = f"mip_high_s{slab}_c{chip}_ch{chan}_sca{s_}"; pn = f"ped_high_s{slab}_c{chip}_ch{chan}_sca{s_}"
    if n not in names or pn not in names: continue
    mh, ph = f.Get(n), f.Get(pn)
    ph.GetXaxis().SetRangeUser(ph.GetMean() - 20, ph.GetMean() + 20); ped = ph.GetMean()
    if ped <= 0: continue
    for kb in range(900):
        y = mh.GetBinContent(kb)
        if y > 0: hm.Fill(int(mh.GetXaxis().GetBinCenter(kb) - ped), y)
m = k.fitLangau(hm.Clone("m0"), True)
par = ROOT.langauParams(hm)
x, y, w = arrays(hm)
sel = (x > -5) & (x < 80)
xx = np.linspace(0.5, 80, 500)
fig, ax = plt.subplots(figsize=(6.4, 4.6))
ax.step(x[sel], y[sel], where="mid", color="#47505E", lw=1.4, label=f"data, {int(hm.Integral())} hits (4 muon runs)")
ax.plot(xx, [ROOT.langauEval(float(v), par) for v in xx], color="#B5541A", lw=2.2, label="Landau ⊗ Gaussian fit")
ax.axvline(m.mpv, color="#B5541A", lw=1, ls="--")
ax.set_xlabel("high-gain ADC − pedestal"); ax.set_ylabel(f"entries / {w:g} ADC")
ax.set_title(f"MIP · slab {slab}, chip {chip}, channel {chan} (all SCAs)", fontsize=13)
ax.text(0.97, 0.95, f"MPV = {m.mpv:.2f} ± {m.empv:.2f} ADC\nLandau width = {m.width:.2f} ADC\n"
        f"Gauss σ = {par[3]:.2f} ADC",
        transform=ax.transAxes, ha="right", va="top", fontsize=12)
ax.legend(loc="center right", frameon=False, fontsize=11)
fig.tight_layout(); fig.savefig(f"{outdir}/mip_fit_s{slab}_c{chip}_ch{chan}.png", dpi=160)
print("refit MPV", par[1]); print(f"ped {mu:.2f} {sg:.2f} status {int(r.status)}  mip {m.mpv:.2f} {m.width:.2f} chi2 {m.chi2ndf:.2f}")
