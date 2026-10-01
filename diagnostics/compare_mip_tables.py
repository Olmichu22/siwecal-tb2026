#!/usr/bin/env python
"""
Compare two sets of calibration tables (e.g. before and after the SCA column-order fix) for every threshold:
MIP MPV and Landau width per channel, pedestal mean and width (RMS) per (channel, SCA), both gains.

    python diagnostics/compare_mip_tables.py --old calibration/MuonCalib_gaudi \
        --new calibration/MuonCalib_gaudi_fixed --out calibration/MuonCalib_gaudi_fixed/compare [--th 210,220,230]

MIP table: '#layer chip channel mpv empv widthmpv chi2ndf nentries' (empv <= 0 or mpv <= 0: no fit).
Writes, per threshold and gain, mip_compare_th<N>_<gain>gain.png, and for all thresholds together
calib_overview_<gain>gain.png (per-slab medians of MPV, MIP width, pedestal, pedestal width, old vs new),
plus mip_compare.txt and calib_overview.txt.
"""
import argparse
import glob
import os
import sys

import numpy as np
np.seterr(divide="ignore", invalid="ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_pedestal_tables import read_table as read_ped  # noqa: E402

C_OLD, C_NEW = "#8c8c8c", "#29a3dc"
TH_MARK = {"210": "o", "220": "s", "230": "^"}


def read_mip(path):
    a = np.loadtxt(path, comments="#")
    key = a[:, :3].astype(int)
    mpv, empv, width, nent = a[:, 3], a[:, 4], a[:, 5], a[:, 7]
    ok = (mpv > 0) & (empv > 0) & (width > 0)
    return key, mpv, width, nent, ok


def table(base, kind, th, gain):
    pre = "MIP_pedestalsubmode1_" if kind == "mips" else "Pedestal_"
    g = sorted(glob.glob(os.path.join(base, kind, f"th{th}", f"{pre}*_{gain}gain.txt")))
    return g[-1] if g else None


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--old", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--th", default="210,220,230")
    p.add_argument("--label-old", default="legacy decoding")
    p.add_argument("--label-new", default="fixed SCA pairing")
    args = p.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    ths = args.th.split(",")
    txt = open(os.path.join(args.out, "mip_compare.txt"), "w")
    ov = open(os.path.join(args.out, "calib_overview.txt"), "w")

    def out(s="", f=txt):
        print(s); f.write(s + "\n")

    overview = {}   # (gain, th, which) -> dict of per-slab medians
    for gain in ("high", "low"):
        for th in ths:
            fo, fn = table(args.old, "mips", th, gain), table(args.new, "mips", th, gain)
            if not (fo and fn):
                out(f"==== th{th}, {gain} gain: MIP table missing ({'old' if not fo else 'new'}), skipped\n")
                continue
            ko, mo, wo, no, oko = read_mip(fo)
            kn, mn, wn, nn, okn = read_mip(fn)
            assert np.array_equal(ko, kn)
            ok = oko & okn
            slab = ko[:, 0]
            r, rw = mn / mo, wn / wo
            out(f"==== th{th}, {gain} gain: MIP ====\n  old: {fo}\n  new: {fn}")
            out(f"  fitted channels: old {oko.sum()}, new {okn.sum()}, both {ok.sum()} of {len(ok)}")
            out(f"  entries per channel (median): old {np.median(no[oko]):.0f}, new {np.median(nn[okn]):.0f}")
            out(f"  {'':<9}{'N':>6}{'MPV old':>9}{'MPV new':>9}{'new/old':>9}{'hw68':>7}{'|r-1|>3%':>10}"
                f"{'w old':>8}{'w new':>8}{'w n/o':>7}")

            def row(lbl, sel):
                q16, q50, q84 = np.percentile(r[sel], [16, 50, 84])
                out(f"  {lbl:<9}{sel.sum():6d}{np.median(mo[sel]):9.2f}{np.median(mn[sel]):9.2f}{q50:9.3f}"
                    f"{0.5 * (q84 - q16):7.3f}{np.mean(np.abs(r[sel] - 1) > 0.03):10.1%}"
                    f"{np.median(wo[sel]):8.2f}{np.median(wn[sel]):8.2f}{np.median(rw[sel]):7.3f}")
            row("all", ok)
            for s in np.unique(slab):
                if (ok & (slab == s)).sum() >= 10:
                    row(f"slab {s}", ok & (slab == s))
            out("")
            slabs = np.unique(slab)
            for which, m, w, okx in (("old", mo, wo, oko), ("new", mn, wn, okn)):
                overview[(gain, th, which)] = {
                    "mpv": [np.median(m[okx & (slab == s)]) if (okx & (slab == s)).sum() else np.nan for s in slabs],
                    "mipw": [np.median(w[okx & (slab == s)]) if (okx & (slab == s)).sum() else np.nan for s in slabs]}

            fig, ax = plt.subplots(2, 3, figsize=(18, 10))
            hi = np.percentile(np.concatenate([mo[ok], mn[ok]]), 99.5)
            b = np.linspace(0, hi, 100)
            for a_, v_o, v_n, xl in ((ax[0, 0], mo, mn, "MIP MPV [ADC]"), (ax[1, 0], wo, wn, "MIP Landau width [ADC]")):
                hi = np.percentile(np.concatenate([v_o[ok], v_n[ok]]), 99.5)
                bb = np.linspace(0, hi, 100)
                a_.hist(v_o[ok], bins=bb, color=C_OLD, histtype="stepfilled", alpha=0.6,
                        label=f"{args.label_old}: median {np.median(v_o[ok]):.2f}")
                a_.hist(v_n[ok], bins=bb, color=C_NEW, histtype="step", lw=2,
                        label=f"{args.label_new}: median {np.median(v_n[ok]):.2f}")
                a_.set_xlabel(xl); a_.set_ylabel("channels"); a_.legend(fontsize=9); a_.grid(alpha=0.25)
            for a_, rr, xl in ((ax[0, 1], r, "MPV ratio new / old"), (ax[1, 1], rw, "Width ratio new / old")):
                bb = np.linspace(0.7, 1.3, 121)
                a_.hist(np.clip(rr[ok], 0.7, 1.3), bins=bb, color=C_NEW, histtype="stepfilled", alpha=0.5)
                a_.hist(np.clip(rr[ok], 0.7, 1.3), bins=bb, color=C_NEW, histtype="step", lw=1.5)
                a_.axvline(1, color="0.4", ls="--", lw=1); a_.set_xlabel(xl); a_.grid(alpha=0.25)
                q16, q50, q84 = np.percentile(rr[ok], [16, 50, 84])
                a_.set_title(f"median {q50:.3f}, 68 % half-width {0.5 * (q84 - q16):.3f} (under/overflow clipped)",
                             fontsize=10)
            for a_, key, yl in ((ax[0, 2], "mpv", "Median MIP MPV [ADC]"), (ax[1, 2], "mipw", "Median Landau width [ADC]")):
                for which, c, lbl in (("old", C_OLD, args.label_old), ("new", C_NEW, args.label_new)):
                    a_.plot(slabs, overview[(gain, th, which)][key], "o-", color=c, lw=2, ms=7, label=lbl)
                a_.set_xlabel("Slab"); a_.set_ylabel(yl); a_.legend(fontsize=9); a_.grid(alpha=0.25)
            fig.suptitle(f"th{th} MIP tables, {gain} gain: {args.label_new} vs {args.label_old}", fontsize=13)
            fig.tight_layout(rect=[0, 0, 1, 0.95])
            fig.savefig(os.path.join(args.out, f"mip_compare_th{th}_{gain}gain.png"), dpi=105); plt.close(fig)

        # pedestals for the overview
        for th in ths:
            fo, fn = table(args.old, "pedestals", th, gain), table(args.new, "pedestals", th, gain)
            if not (fo and fn):
                continue
            ko, mo, wo, oko = read_ped(fo)
            kn, mn, wn, okn = read_ped(fn)
            slab = ko[:, 0]
            slabs = np.unique(slab)
            for which, m, w, okx in (("old", mo, wo, oko), ("new", mn, wn, okn)):
                d = overview.setdefault((gain, th, which), {})
                d["ped"] = [np.median(m[okx & (slab == s)[:, None]]) for s in slabs]
                d["pedw"] = [np.median(w[okx & (slab == s)[:, None]]) for s in slabs]

    # ---- overview: per-slab medians of the four quantities, all thresholds, old vs new ----------------
    Q = [("mpv", "MIP MPV [ADC]"), ("mipw", "MIP Landau width [ADC]"), ("ped", "Pedestal mean [ADC]"),
         ("pedw", "Pedestal width (RMS) [ADC]")]
    for gain in ("high", "low"):
        fig, ax = plt.subplots(2, 2, figsize=(15, 10))
        out(f"==== per-slab medians, {gain} gain (old / new) ====", ov)
        for a_, (key, yl) in zip(ax.ravel(), Q):
            for th in ths:
                for which, c, ls in (("old", C_OLD, "--"), ("new", C_NEW, "-")):
                    d = overview.get((gain, th, which), {})
                    if key not in d:
                        continue
                    a_.plot(range(len(d[key])), d[key], ls=ls, marker=TH_MARK.get(th, "o"), color=c, lw=1.6, ms=7,
                            mfc=c if which == "new" else "white", label=f"th{th}, {args.label_old if which == 'old' else args.label_new}")
            a_.set_xlabel("Slab"); a_.set_ylabel(yl); a_.grid(alpha=0.25)
        ax[0, 0].legend(fontsize=8, ncol=2)
        fig.suptitle(f"Calibration tables, {gain} gain: per-slab medians, every threshold "
                     f"(dashed/open: {args.label_old}; solid/filled: {args.label_new})", fontsize=12.5)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(os.path.join(args.out, f"calib_overview_{gain}gain.png"), dpi=105); plt.close(fig)
        for key, yl in Q:
            out(f"  {yl}", ov)
            out("  " + f"{'slab':<6}" + "".join(f"{'th' + th + ' old':>11}{'th' + th + ' new':>11}" for th in ths), ov)
            n = max(len(overview.get((gain, th, 'old'), {}).get(key, [])) for th in ths)
            for s in range(n):
                cells = ""
                for th in ths:
                    for which in ("old", "new"):
                        v = overview.get((gain, th, which), {}).get(key, [])
                        cells += f"{v[s]:11.2f}" if s < len(v) else f"{'-':>11}"
                out(f"  {s:<6}" + cells, ov)
            out("", ov)
    txt.close(); ov.close()


if __name__ == "__main__":
    main()
