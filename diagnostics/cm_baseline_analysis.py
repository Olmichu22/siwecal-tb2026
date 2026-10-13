"""Figures and numbers of the common-mode study, from diagnostics/cm_baseline_collect.py outputs (cm_<run>.npz).

    python diagnostics/cm_baseline_analysis.py <indir> <figdir>

Baseline = ADC - pedestal of the channels WITHOUT a hit bit. Two estimators:
  mode    peak of the per-channel residual histogram (+/- 1.5 ADC centroid), per nbits bin: the baseline itself,
          not pulled up by sub-threshold signal
  median  median over chip-SCAs of the per-chip-SCA median of the 'far' channels (>= 2 pads from any hit), used
          where the data are split in two variables (nbits and charge) and a histogram per cell would be too thin
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

RUNS = [("e151", "muons (th210)"), ("e293", "7.5 GeV"), ("e292", "10 GeV"), ("e291", "20 GeV"), ("e288", "34 GeV"),
        ("e287", "52 GeV"), ("e286", "74 GeV"), ("e166", "99 GeV")]
# sequential blue ramp, light -> dark with energy; muons in neutral grey
COL = {"e151": "#8a8f98", "e293": "#9ecae9", "e292": "#78b3df", "e291": "#4f97cf", "e288": "#2f7bbb",
       "e287": "#1d60a0", "e286": "#134682", "e166": "#0b2c5c"}
INK, MUTED, GRID = "#1b2227", "#56626b", "#d7dde2"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7,
                     "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})


def mode_of(h, bins):
    ctr = 0.5 * (bins[1:] + bins[:-1])
    if h.sum() < 200:
        return np.nan
    i = np.argmax(h)
    w = slice(max(i - 6, 0), i + 7)
    return float(np.sum(ctr[w] * h[w]) / max(h[w].sum(), 1))


def load(indir):
    out = {}
    for r, lab in RUNS:
        p = f"{indir}/cm_{r}.npz"
        if not os.path.exists(p):
            continue
        z = np.load(p)
        cols = list(z["cols"])
        rec = {c: z["rec"][:, i] for i, c in enumerate(cols)}
        out[r] = dict(label=lab, rec=rec, nb=z["nb_edges"], hb=z["hbins"], n_acq=int(z["n_acq"]),
                      h={k: z[f"h_{k}"] for k in ("all", "near", "far", "lg_far")})
    return out


def nb_centres(nb):
    return np.array([0.5 * (lo + hi - 1) for lo, hi in zip(nb[:-1], nb[1:])])


def binned_median(x, y, edges, nmin=30):
    m, c = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = (x >= lo) & (x < hi) & np.isfinite(y)
        m.append(np.median(y[k]) if k.sum() >= nmin else np.nan)
        c.append(np.median(x[k]) if k.sum() >= nmin else np.nan)
    return np.array(c), np.array(m)


def main():
    indir, figdir = sys.argv[1], sys.argv[2]
    os.makedirs(figdir, exist_ok=True)
    D = load(indir)
    lines = []
    def out(s=""):
        print(s); lines.append(s)

    # ---------------------------------------------------------------- F1: mode vs nbits, every energy
    fig, ax = plt.subplots(figsize=(8.6, 5.4))
    out("Baseline (mode, far channels, high gain) [ADC] by hit bits in the chip-SCA")
    nbc = None
    for r, d in D.items():
        nbc = nb_centres(d["nb"])
        mo = np.array([mode_of(h, d["hb"]) for h in d["h"]["far"]])
        ax.plot(nbc, mo, "o-", color=COL[r], lw=1.8, ms=5, label=d["label"])
        out(f"  {d['label']:<14}" + " ".join(f"{m:+5.1f}" if np.isfinite(m) else "   . " for m in mo))
    out("  bins: " + " ".join(f"{int(lo)}-{int(hi) - 1}" for lo, hi in zip(nbc is not None and D[next(iter(D))]["nb"][:-1], D[next(iter(D))]["nb"][1:])))
    ax.axhline(0, color=MUTED, lw=0.8)
    ax.set_xlabel("Channels with a hit bit in the chip-SCA (of 64)")
    ax.set_ylabel("Baseline of the non-fired channels [ADC]\n(mode of ADC − pedestal, ≥ 2 pads from a hit)")
    ax.legend(ncol=2, fontsize=9.5, loc="lower left")
    fig.tight_layout(); fig.savefig(f"{figdir}/f1_mode_vs_nbits.png", dpi=130); plt.close(fig)

    # ---------------------------------------------------------------- F2: charge or channel count?
    allrec = {k: np.concatenate([d["rec"][k] for r, d in D.items() if r != "e151"]) for k in D[next(iter(D))]["rec"]}
    nbits, q, bf = allrec["nbits"], allrec["q"], allrec["base_far"]
    QE = np.array([0, 500, 1000, 2000, 4000, 8000, 16000, 1e9])
    NE = np.array([1, 5, 10, 15, 20, 30, 45])
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2))
    qcol = plt.cm.Blues(np.linspace(0.35, 1.0, len(QE) - 1))
    ncol = plt.cm.Oranges(np.linspace(0.35, 1.0, len(NE) - 1))
    nb_fine = np.arange(1, 50, 3)
    out("\nSplit by charge (median of per-chip-SCA far medians) [ADC]")
    for i, (lo, hi) in enumerate(zip(QE[:-1], QE[1:])):
        k = (q >= lo) & (q < hi)
        c, m = binned_median(nbits[k], bf[k], nb_fine)
        lab = f"Q {lo / 1000:g}–{hi / 1000:g} k ADC" if hi < 1e8 else f"Q > {lo / 1000:g} k ADC"
        axs[0].plot(c, m, "o-", color=qcol[i], ms=4, lw=1.6, label=lab)
    q_fine = np.geomspace(100, 60000, 18)
    for i, (lo, hi) in enumerate(zip(NE[:-1], NE[1:])):
        k = (nbits >= lo) & (nbits < hi)
        c, m = binned_median(q[k], bf[k], q_fine)
        axs[1].plot(c, m, "s-", color=ncol[i], ms=4, lw=1.6, label=f"{lo}–{hi - 1} hit bits")
        out(f"  {lo:2d}-{hi - 1:2d} bits: " + " ".join(f"{x:.0f}:{y:+.1f}" for x, y in zip(c, m) if np.isfinite(y)))
    axs[0].set_xlabel("Channels with a hit bit"); axs[0].set_ylabel("Baseline, far channels [ADC] (median)")
    axs[0].set_title("At fixed charge, more fired channels", fontsize=11.5)
    axs[1].set_xscale("log")
    axs[1].set_xlabel("Charge of the fired channels [high-gain ADC, from the low gain]")
    axs[1].set_title("At fixed fired channels, more charge", fontsize=11.5)
    for a in axs:
        a.axhline(0, color=MUTED, lw=0.8); a.legend(fontsize=8.8)
    fig.tight_layout(); fig.savefig(f"{figdir}/f2_charge_vs_count.png", dpi=130); plt.close(fig)

    # linear model base = a + b*q + c*nbits on busy-ish chip-SCAs (fit quality reported)
    k = np.isfinite(bf) & (nbits >= 1) & (q > 0) & (q < 60000)
    X = np.stack([np.ones(k.sum()), q[k] / 1000, nbits[k]], 1)
    coef, *_ = np.linalg.lstsq(X, bf[k], rcond=None)
    out(f"\nLinear fit base_far = a + b*Q[k ADC] + c*nbits on {k.sum()} chip-SCAs: a = {coef[0]:+.2f}, "
        f"b = {coef[1]:+.3f} ADC per 1000 ADC of charge, c = {coef[2]:+.3f} ADC per fired channel")
    X1 = X[:, :2]; c1, *_ = np.linalg.lstsq(X1, bf[k], rcond=None)
    X2 = X[:, [0, 2]]; c2, *_ = np.linalg.lstsq(X2, bf[k], rcond=None)
    rms = lambda c, X: float(np.sqrt(np.mean((bf[k] - X @ c) ** 2)))  # noqa: E731
    out(f"  residual RMS: charge+count {rms(coef, X):.3f}, charge only {rms(c1, X1):.3f} (b = {c1[1]:+.3f}), "
        f"count only {rms(c2, X2):.3f} (c = {c2[1]:+.3f}) [ADC]")

    # ---------------------------------------------------------------- F3: near vs far, and the low gain
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.0))
    out("\nNear (adjacent pad) minus far baseline, and low gain (far, in high-gain ADC via /k), mode [ADC]")
    for r in ("e291", "e286"):
        if r not in D:
            continue
        d = D[r]; nbc = nb_centres(d["nb"])
        mf = np.array([mode_of(h, d["hb"]) for h in d["h"]["far"]])
        mn = np.array([mode_of(h, d["hb"]) for h in d["h"]["near"]])
        ml = np.array([mode_of(h, d["hb"]) for h in d["h"]["lg_far"]])
        axs[0].plot(nbc, mf, "o-", color=COL[r], lw=1.8, label=f"{d['label']}: ≥ 2 pads from a hit")
        axs[0].plot(nbc, mn, "o--", color=COL[r], lw=1.4, mfc="white", label=f"{d['label']}: adjacent to a hit")
        pass
        out(f"  {d['label']:<8} near-far " + " ".join(f"{a - b:+.1f}" if np.isfinite(a - b) else "  . " for a, b in zip(mn, mf)))
    # high gain against low gain, both against the chip charge; LG values are integers, so per-chip-SCA medians
    # averaged over chip-SCAs (a mode would sit on the integer steps)
    pool = [r for r in ("e291", "e287", "e286", "e166") if r in D]
    qq = np.concatenate([D[r]["rec"]["q"] for r in pool]); nbq = np.concatenate([D[r]["rec"]["nbits"] for r in pool])
    hgq = np.concatenate([D[r]["rec"]["base_far"] for r in pool]); lgq = np.concatenate([D[r]["rec"]["lg_far"] for r in pool])
    QB = np.array([0, 250, 500, 1000, 2000, 4000, 8000, 16000])
    kk = (nbq >= 1) & np.isfinite(hgq) & np.isfinite(lgq)
    xc, yh, yl = [], [], []
    for lo, hi in zip(QB[:-1], QB[1:]):
        m = kk & (qq >= lo) & (qq < hi)
        xc.append(np.median(qq[m])); yh.append(np.mean(hgq[m])); yl.append(np.mean(lgq[m]) / 0.0924)
    yh, yl = np.array(yh), np.array(yl)
    axs[1].plot(xc, yh - yh[0], "o-", color=COL["e286"], lw=1.8, label="high gain")
    axs[1].plot(xc, yl - yl[0], "D--", color="#c0392b", lw=1.6, mfc="white", label="low gain ÷ k (in high-gain ADC)")
    axs[1].set_xscale("log"); axs[1].set_xlabel("Charge of the fired channels [high-gain ADC]")
    for a, t in zip(axs, ("Common mode or crosstalk: distance to the nearest hit",
                          "High gain against low gain (20–99 GeV pooled)")):
        a.axhline(0, color=MUTED, lw=0.8); a.set_title(t, fontsize=11.5); a.legend(fontsize=8.8, loc="lower left")
    axs[0].set_xlabel("Channels with a hit bit in the chip-SCA")
    axs[0].set_ylabel("Baseline [ADC] (mode)")
    axs[1].set_ylabel("Baseline change from the quietest bin\n[high-gain ADC] (mean of chip-SCA medians)")
    fig.tight_layout(); fig.savefig(f"{figdir}/f3_near_far_lg.png", dpi=130); plt.close(fig)

    # ---------------------------------------------------------------- F4: SCA and slab, busy chips at 74 GeV
    if "e286" in D:
        rec = D["e286"]["rec"]
        busy = (rec["nbits"] >= 15) & (rec["nbits"] < 35)
        quiet = (rec["nbits"] >= 1) & (rec["nbits"] < 4)
        fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
        out("\n74 GeV, median far baseline [ADC] by SCA (busy 15-34 bits | quiet 1-3 bits) and by slab (busy)")
        for sel, lab, col in ((busy, "15–34 hit bits", COL["e286"]), (quiet, "1–3 hit bits", COL["e293"])):
            c, m = binned_median(rec["sca"][sel], rec["base_far"][sel], np.arange(-0.5, 15.5, 1), nmin=20)
            axs[0].plot(np.arange(15), m, "o-", color=col, label=lab)
            c2, m2 = binned_median(rec["slab"][sel], rec["base_far"][sel], np.arange(-0.5, 15.5, 1), nmin=20)
            axs[1].plot(np.arange(15), m2, "o-", color=col, label=lab)
            out(f"  {lab:<14} SCA : " + " ".join(f"{x:+.1f}" if np.isfinite(x) else " . " for x in m))
            out(f"  {lab:<14} slab: " + " ".join(f"{x:+.1f}" if np.isfinite(x) else " . " for x in m2))
        axs[0].set_xlabel("SCA (memory cell, in order of filling)"); axs[1].set_xlabel("Slab")
        for a in axs:
            a.axhline(0, color=MUTED, lw=0.8); a.set_ylabel("Baseline, far channels [ADC] (median)"); a.legend(fontsize=9)
        fig.tight_layout(); fig.savefig(f"{figdir}/f4_sca_slab.png", dpi=130); plt.close(fig)

    # ---------------------------------------------------------------- F5: the residual distributions themselves
    if "e286" in D:
        d = D["e286"]; hb = d["hb"]; ctr = 0.5 * (hb[1:] + hb[:-1])
        fig, ax = plt.subplots(figsize=(8.6, 5.0))
        pick = [0, 2, 5, 8, 10, 11]
        shades = plt.cm.Blues(np.linspace(0.35, 1.0, len(pick)))
        for j, i in enumerate(pick):
            h = d["h"]["far"][i]
            if h.sum() < 500:
                continue
            lo, hi = d["nb"][i], d["nb"][i + 1] - 1
            h2 = h[: len(h) // 4 * 4].reshape(-1, 4).sum(1); c2 = ctr[: len(h) // 4 * 4].reshape(-1, 4).mean(1)   # 1 ADC bins
            ax.plot(c2, h2 / h2.sum(), color=shades[j], lw=1.8, drawstyle="steps-mid",
                    label=f"{lo}–{hi} hit bits (mode {mode_of(h, hb):+.1f})" if hi > lo else f"{lo} hit bits (mode {mode_of(h, hb):+.1f})")
        ax.set_xlim(-15, 15); ax.set_xlabel("ADC − pedestal, non-fired channels ≥ 2 pads from a hit [ADC]")
        ax.set_ylabel("Channels (unit area)"); ax.axvline(0, color=MUTED, lw=0.8); ax.legend(fontsize=9)
        fig.tight_layout(); fig.savefig(f"{figdir}/f5_distributions_74.png", dpi=130); plt.close(fig)

    # ---------------------------------------------------------------- F6: the LG->HG anchor against occupancy
    K0, C0 = 0.0924, 2.35
    AN = [(r, lab) for r, lab in RUNS if os.path.exists(f"{indir}/anchor_{r}.npy")]
    if AN:
        OB = [(1, 5), (5, 10), (10, 20), (20, 30), (30, 65)]
        sl_edges = np.arange(150, 1201, 50)   # the official fiducial band (plot_gain_anchor_fit.py)
        fig, axs = plt.subplots(1, 2, figsize=(13, 5.0))
        out("\nLG->HG anchor by hit bits in the chip-SCA: line through the median adc_low in 50-ADC slices of "
            "adc_high (150-1200, the official fiducial band), fired channels")
        out(f"  reconstruction uses k = {K0}, c = {C0}; bias = energy of a saturated hit read through (k0, c0) "
            "when the chip's real line is (k, c), at 2000 / 4000 HG ADC")
        mk = {"e291": "o", "e287": "s", "e286": "D"}
        for r, lab in AN:
            arr = np.load(f"{indir}/anchor_{r}.npy")
            hg, lg, nb = arr[:, 0], arr[:, 1], arr[:, 2]
            cs, ks, xs = [], [], []
            for lo, hi in OB:
                m = (nb >= lo) & (nb < hi)
                c, med = binned_median(hg[m], lg[m], sl_edges, nmin=40)
                ok = np.isfinite(med)
                if ok.sum() < 5:
                    cs.append(np.nan); ks.append(np.nan); xs.append(0.5 * (lo + hi - 1)); continue
                k_, c_ = np.polyfit(c[ok], med[ok], 1)
                cs.append(c_); ks.append(k_); xs.append(0.5 * (lo + hi - 1))
                b2 = ((k_ * 2000 + c_ - C0) / K0) / 2000 - 1
                b4 = ((k_ * 4000 + c_ - C0) / K0) / 4000 - 1
                out(f"  {lab:<8} {lo:2d}-{hi - 1:2d} bits: k = {k_:.4f}  c = {c_:+.2f} LG ADC  (N = {m.sum():7d})  "
                    f"bias {100 * b2:+.1f} % / {100 * b4:+.1f} %")
            axs[0].plot(xs, cs, mk.get(r, "o") + "-", color=COL[r], ms=6, lw=1.6, label=lab)
            axs[1].plot(xs, np.array(ks), mk.get(r, "o") + "-", color=COL[r], ms=6, lw=1.6, label=lab)
        axs[0].axhline(C0, color=MUTED, lw=1, ls="--", label=f"reconstruction, c = {C0}")
        axs[1].axhline(K0, color=MUTED, lw=1, ls="--", label=f"reconstruction, k = {K0}")
        axs[0].set_ylabel("Anchor intercept c [LG ADC]"); axs[1].set_ylabel("Anchor slope k [LG ADC / HG ADC]")
        for a_ in axs:
            a_.set_xlabel("Channels with a hit bit in the chip-SCA"); a_.legend(fontsize=9)
        fig.tight_layout(); fig.savefig(f"{figdir}/f6_anchor_vs_occupancy.png", dpi=130); plt.close(fig)

    out("\nAcquisitions per run: " + ", ".join(f"{d['label']} {d['n_acq']}" for d in D.values()))
    open(f"{figdir}/cm_summary.txt", "w").write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
