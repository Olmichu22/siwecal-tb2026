#!/usr/bin/env python
"""
Compare two pedestal tables channel by channel and SCA by SCA (e.g. th210 before and after the SCA column-order
fix): mean, width (RMS of the Gaussian fit), per slab and per SCA, both gains.

    python diagnostics/compare_pedestal_tables.py --old calibration/MuonCalib_gaudi \
        --new calibration/MuonCalib_gaudi_fixed --th 210 --out calibration/MuonCalib_gaudi_fixed/compare

Table format: '#layer chip channel' then (ped_mean, ped_error, ped_width) per SCA; ped_error <= 0 marks an
unfitted SCA. Only (channel, SCA) cells valid in BOTH tables are compared.

Writes pedestal_{compare,width,sca}_th<N>_<gain>gain.png and pedestal_compare_th<N>.txt (tables: per slab, per SCA, global).
"""
import argparse
import glob
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

C_OLD, C_NEW = "#8c8c8c", "#29a3dc"


def read_table(path):
    a = np.loadtxt(path, comments="#")
    key = a[:, :3].astype(int)
    v = a[:, 3:]
    nsca = v.shape[1] // 3
    mean, err, width = v[:, 0::3][:, :nsca], v[:, 1::3][:, :nsca], v[:, 2::3][:, :nsca]
    ok = (err > 0) & (width > 0) & (mean > 0)
    return key, mean, width, ok


def find(base, th, gain):
    g = sorted(glob.glob(os.path.join(base, "pedestals", f"th{th}", f"Pedestal_*_{gain}gain.txt")))
    if not g:
        raise SystemExit(f"no {gain}-gain pedestal table in {base}/pedestals/th{th}")
    return g[-1]


def stats(d):
    d = d[np.isfinite(d)]
    if not len(d):
        return (0, np.nan, np.nan, np.nan, np.nan)
    q16, q50, q84 = np.percentile(d, [16, 50, 84])
    return (len(d), q50, 0.5 * (q84 - q16), float(np.sqrt(np.mean(d ** 2))), float(np.mean(np.abs(d) > 0.5)))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--old", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--th", default="210")
    p.add_argument("--out", required=True)
    p.add_argument("--label-old", default="legacy decoding")
    p.add_argument("--label-new", default="fixed SCA pairing")
    args = p.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    txt = open(os.path.join(args.out, f"pedestal_compare_th{args.th}.txt"), "w")

    def out(s=""):
        print(s); txt.write(s + "\n")

    for gain in ("high", "low"):
        fo, fn = find(args.old, args.th, gain), find(args.new, args.th, gain)
        ko, mo, wo, oko = read_table(fo)
        kn, mn, wn, okn = read_table(fn)
        assert np.array_equal(ko, kn), "tables list different channels"
        ok = oko & okn
        nsca = mo.shape[1]
        slab = ko[:, 0]
        dm = np.where(ok, mn - mo, np.nan)
        dw = np.where(ok, wn - wo, np.nan)
        # SCA structure: each cell's pedestal minus its channel's mean over its valid SCAs
        def sca_offsets(m):
            mm = np.where(ok, m, np.nan)
            return mm - np.nanmean(mm, axis=1, keepdims=True)
        so, sn = sca_offsets(mo), sca_offsets(mn)

        out(f"==== th{args.th}, {gain} gain ====")
        out(f"  old: {fo}\n  new: {fn}")
        out(f"  valid (channel, SCA) cells: old {oko.sum()}, new {okn.sum()}, both {ok.sum()} of {ok.size}")
        out("")
        out("  Pedestal MEAN difference new - old [ADC]")
        out(f"  {'':<10}{'N':>8}{'median':>9}{'hw68':>8}{'RMS':>8}{'|d|>0.5':>9}")
        out(f"  {'all':<10}" + "{:8d}{:9.2f}{:8.2f}{:8.2f}{:9.1%}".format(*stats(dm[ok])))
        for s in np.unique(slab):
            out(f"  {'slab ' + str(s):<10}" + "{:8d}{:9.2f}{:8.2f}{:8.2f}{:9.1%}".format(*stats(dm[slab == s][ok[slab == s]])))
        out("")
        out("  Pedestal WIDTH (RMS) [ADC]: median old, median new, median ratio new/old")
        out(f"  {'':<10}{'N':>8}{'old':>9}{'new':>9}{'new/old':>9}")

        def wrow(lbl, sel):
            a, b = wo[sel], wn[sel]
            out(f"  {lbl:<10}{sel.sum():8d}{np.median(a):9.3f}{np.median(b):9.3f}{np.median(b / a):9.3f}")
        wrow("all", ok)
        for s in np.unique(slab):
            wrow(f"slab {s}", ok & (slab[:, None] == s))
        out("")
        out("  Per SCA: spread of the SCA offsets (pedestal - channel's mean over SCAs), and widths")
        out(f"  {'SCA':<5}{'N':>8}{'offset old':>12}{'offset new':>12}{'|dmean|':>9}{'w old':>8}{'w new':>8}")
        for k in range(nsca):
            sel = ok[:, k]
            if sel.sum() < 10:
                continue
            out(f"  {k:<5}{sel.sum():8d}{np.median(so[sel, k]):12.2f}{np.median(sn[sel, k]):12.2f}"
                f"{np.median(np.abs(dm[sel, k])):9.2f}{np.median(wo[sel, k]):8.3f}{np.median(wn[sel, k]):8.3f}")
        out(f"  RMS of the SCA offsets over all cells: old {np.nanstd(so[ok]):.2f} ADC, new {np.nanstd(sn[ok]):.2f} ADC")
        out("")

        # ---- figure 1: mean -------------------------------------------------------------------------
        fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
        d = dm[ok]
        lo, hi = np.percentile(d, [0.5, 99.5])
        ax[0].hist(d, bins=np.linspace(lo, hi, 120), color=C_NEW, histtype="stepfilled", alpha=0.5)
        ax[0].hist(d, bins=np.linspace(lo, hi, 120), color=C_NEW, histtype="step", lw=1.5)
        ax[0].set_xlabel("Pedestal mean, new − old [ADC]"); ax[0].set_ylabel("(channel, SCA) cells")
        ax[0].set_yscale("log"); ax[0].grid(alpha=0.25)
        n, med, hw, rms, frac = stats(d)
        ax[0].set_title(f"median {med:+.2f}, RMS {rms:.2f}, |Δ| > 0.5 ADC: {frac:.0%}", fontsize=10)
        slabs = np.unique(slab)
        for lbl, m, c in ((args.label_old, mo, C_OLD), (args.label_new, mn, C_NEW)):
            med_s = [np.median(m[(slab == s)[:, None] & ok]) for s in slabs]
            ax[1].plot(slabs, med_s, "o-", color=c, lw=2, ms=7, label=lbl)
        ax[1].set_xlabel("Slab"); ax[1].set_ylabel("Median pedestal mean [ADC]"); ax[1].grid(alpha=0.25)
        ax[1].legend(fontsize=9)
        b = np.linspace(-8, 8, 161)
        for lbl, so_, c, st in ((args.label_old, so, C_OLD, "stepfilled"), (args.label_new, sn, C_NEW, "step")):
            ax[2].hist(so_[ok], bins=b, color=c, histtype=st, alpha=0.6 if st == "stepfilled" else 1, lw=2,
                       label=f"{lbl}: RMS {np.nanstd(so_[ok]):.2f} ADC")
        ax[2].set_xlabel("SCA offset: pedestal − channel mean over SCAs [ADC]"); ax[2].grid(alpha=0.25)
        ax[2].legend(fontsize=9)
        fig.suptitle(f"th{args.th} pedestal tables, {gain} gain: {args.label_new} vs {args.label_old}", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(os.path.join(args.out, f"pedestal_compare_th{args.th}_{gain}gain.png"), dpi=110); plt.close(fig)

        # ---- figure 2: width (RMS) ------------------------------------------------------------------
        fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))
        hiw = np.percentile(np.concatenate([wo[ok], wn[ok]]), 99.5)
        bw = np.linspace(0, hiw, 120)
        ax[0].hist(wo[ok], bins=bw, color=C_OLD, histtype="stepfilled", alpha=0.6,
                   label=f"{args.label_old}: median {np.median(wo[ok]):.2f}")
        ax[0].hist(wn[ok], bins=bw, color=C_NEW, histtype="step", lw=2,
                   label=f"{args.label_new}: median {np.median(wn[ok]):.2f}")
        ax[0].set_xlabel("Pedestal width (RMS) [ADC]"); ax[0].set_ylabel("(channel, SCA) cells")
        ax[0].legend(fontsize=9); ax[0].grid(alpha=0.25)
        r = (wn / wo)[ok]
        ax[1].hist(r, bins=np.linspace(0, 2, 121), color=C_NEW, histtype="stepfilled", alpha=0.5)
        ax[1].hist(r, bins=np.linspace(0, 2, 121), color=C_NEW, histtype="step", lw=1.5)
        ax[1].axvline(1, color="0.4", ls="--", lw=1)
        ax[1].set_xlabel("Width ratio new / old"); ax[1].grid(alpha=0.25)
        ax[1].set_title(f"median {np.median(r):.3f}", fontsize=10)
        for lbl, w, c in ((args.label_old, wo, C_OLD), (args.label_new, wn, C_NEW)):
            ax[2].plot(slabs, [np.median(w[(slab == s)[:, None] & ok]) for s in slabs], "o-", color=c, lw=2, ms=7,
                       label=lbl)
        ax[2].set_xlabel("Slab"); ax[2].set_ylabel("Median pedestal width [ADC]"); ax[2].grid(alpha=0.25)
        ax[2].legend(fontsize=9)
        fig.suptitle(f"th{args.th} pedestal widths (RMS), {gain} gain", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(os.path.join(args.out, f"pedestal_width_th{args.th}_{gain}gain.png"), dpi=110); plt.close(fig)

        # ---- figure 3: per SCA ----------------------------------------------------------------------
        fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
        scas = [k for k in range(nsca) if ok[:, k].sum() >= 10]
        for lbl, so_, w, c in ((args.label_old, so, wo, C_OLD), (args.label_new, sn, wn, C_NEW)):
            q = np.array([np.percentile(so_[ok[:, k], k], [16, 50, 84]) for k in scas])
            ax[0].fill_between(scas, q[:, 0], q[:, 2], color=c, alpha=0.2, lw=0)
            ax[0].plot(scas, q[:, 1], "o-", color=c, lw=2, ms=6, label=f"{lbl} (median, 16-84 % band)")
            ax[1].plot(scas, [np.median(w[ok[:, k], k]) for k in scas], "o-", color=c, lw=2, ms=6, label=lbl)
        ax[0].set_xlabel("SCA"); ax[0].set_ylabel("SCA offset [ADC]"); ax[0].legend(fontsize=9); ax[0].grid(alpha=0.25)
        ax[1].set_xlabel("SCA"); ax[1].set_ylabel("Median pedestal width [ADC]"); ax[1].legend(fontsize=9)
        ax[1].grid(alpha=0.25)
        fig.suptitle(f"th{args.th} pedestals per SCA, {gain} gain", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        fig.savefig(os.path.join(args.out, f"pedestal_sca_th{args.th}_{gain}gain.png"), dpi=110); plt.close(fig)
    txt.close()


if __name__ == "__main__":
    main()
