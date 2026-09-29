#!/usr/bin/env python
"""
Compare the calibration tables of the three threshold sets (th210, th220, th230) with each other, from ONE table
tree (e.g. the fixed-pairing MuonCalib_gaudi_fixed): are the pedestals and the MIPs the same across thresholds?

    python diagnostics/compare_thresholds.py --calib calibration/MuonCalib_gaudi_fixed \
        --out calibration/MuonCalib_gaudi_fixed/compare [--tag fixed]

The preamplifier and the ADC do not know the trigger DAC, so the pedestal of a (channel, SCA) should only differ
between sets by drift between the runs; the MIP MPV, fitted on hits that fired, is cut from the left by the
discriminator and so rises with the threshold (the reason one gain, 19.5 ADC/MIP, is measured at th210).

Pairs compared: th220 - th210, th230 - th210, th230 - th220, on cells valid in both tables.
Writes thresholds_<tag>.txt and thresholds_ped_<tag>_<gain>gain.png, thresholds_mip_<tag>_<gain>gain.png.
"""
import argparse
import itertools
import os
import sys

import numpy as np
np.seterr(divide="ignore", invalid="ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_pedestal_tables import read_table as read_ped  # noqa: E402
from compare_mip_tables import read_mip, table  # noqa: E402

THS = ("210", "220", "230")
COLOR = {"210": "#29a3dc", "220": "#eda100", "230": "#c0392b"}      # validated categorical order
PAIR_COLOR = {("210", "220"): "#29a3dc", ("210", "230"): "#c0392b", ("220", "230"): "#eda100"}
MARK = {"210": "o", "220": "s", "230": "^"}


def pct(d):
    q16, q50, q84 = np.percentile(d, [16, 50, 84])
    return q50, 0.5 * (q84 - q16)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calib", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--tag", default="fixed")
    p.add_argument("--label", default="fixed SCA pairing")
    args = p.parse_args(argv)
    os.makedirs(args.out, exist_ok=True)
    txt = open(os.path.join(args.out, f"thresholds_{args.tag}.txt"), "w")

    def out(s=""):
        print(s); txt.write(s + "\n")

    out(f"Calibration tables across threshold sets, {args.label} ({args.calib})")
    for gain in ("high", "low"):
        P, M = {}, {}
        for th in THS:
            fp, fm = table(args.calib, "pedestals", th, gain), table(args.calib, "mips", th, gain)
            if fp:
                P[th] = read_ped(fp)
            if fm:
                M[th] = read_mip(fm)
        have = [th for th in THS if th in P and th in M]
        out(f"\n==== {gain} gain: thresholds available {', '.join('th' + t for t in have)} ====")
        if len(have) < 2:
            continue
        slab = P[have[0]][0][:, 0]
        slabs = np.unique(slab)
        pairs = list(itertools.combinations(have, 2))

        # ---------------- pedestals ----------------
        out("\n  PEDESTAL per (channel, SCA): median of each set, and pairwise differences [ADC]")
        out(f"  {'set':<8}{'cells':>8}{'mean med':>10}{'width med':>11}")
        for th in have:
            _, m, w, ok = P[th]
            out(f"  th{th:<6}{ok.sum():8d}{np.median(m[ok]):10.2f}{np.median(w[ok]):11.3f}")
        out(f"  {'pair':<14}{'cells':>8}{'dmean med':>11}{'hw68':>7}{'RMS':>7}{'|d|>0.5':>9}{'|d|>2':>8}"
            f"{'w ratio':>9}")
        for a, b in pairs:
            _, ma, wa, oka = P[a]
            _, mb, wb, okb = P[b]
            ok = oka & okb
            d = (mb - ma)[ok]
            med, hw = pct(d)
            out(f"  th{b}-th{a:<7}{ok.sum():8d}{med:11.2f}{hw:7.2f}{np.sqrt(np.mean(d ** 2)):7.2f}"
                f"{np.mean(np.abs(d) > 0.5):9.1%}{np.mean(np.abs(d) > 2):8.1%}{np.median((wb / wa)[ok]):9.3f}")
        out("  per slab: median pedestal mean [ADC] and median width [ADC]")
        out(f"  {'slab':<6}" + "".join(f"{'th' + t + ' mean':>12}" for t in have) + "".join(f"{'th' + t + ' w':>10}" for t in have)
            + "".join(f"{'d' + b + '-' + a:>10}" for a, b in pairs))
        for s in slabs:
            row = f"  {s:<6}"
            for th in have:
                _, m, w, ok = P[th]
                row += f"{np.median(m[ok & (slab == s)[:, None]]):12.2f}"
            for th in have:
                _, m, w, ok = P[th]
                row += f"{np.median(w[ok & (slab == s)[:, None]]):10.3f}"
            for a, b in pairs:
                ok = P[a][3] & P[b][3] & (slab == s)[:, None]
                row += f"{np.median((P[b][1] - P[a][1])[ok]):10.2f}"
            out(row)

        # ---------------- MIPs ----------------
        out("\n  MIP per channel: MPV and Landau width")
        out(f"  {'set':<8}{'fitted':>8}{'MPV med':>9}{'width med':>11}{'width/MPV':>11}")
        for th in have:
            _, mpv, w, _, ok = M[th]
            out(f"  th{th:<6}{ok.sum():8d}{np.median(mpv[ok]):9.2f}{np.median(w[ok]):11.2f}{np.median((w / mpv)[ok]):11.3f}")
        out(f"  {'pair':<14}{'chans':>8}{'MPV ratio':>11}{'hw68':>7}{'width ratio':>13}")
        for a, b in pairs:
            ok = M[a][4] & M[b][4]
            r = (M[b][1] / M[a][1])[ok]
            med, hw = pct(r)
            out(f"  th{b}/th{a:<7}{ok.sum():8d}{med:11.3f}{hw:7.3f}{np.median((M[b][2] / M[a][2])[ok]):13.3f}")
        out("  per slab: median MPV [ADC], median Landau width [ADC]")
        out(f"  {'slab':<6}" + "".join(f"{'th' + t + ' MPV':>11}" for t in have) + "".join(f"{'th' + t + ' w':>9}" for t in have))
        for s in slabs:
            row = f"  {s:<6}"
            for key in (1, 2):
                for th in have:
                    sel = M[th][4] & (M[th][0][:, 0] == s)
                    row += (f"{np.median(M[th][key][sel]):{11 if key == 1 else 9}.2f}" if sel.sum()
                            else f"{'-':>{11 if key == 1 else 9}}")
            out(row)

        # ---------------- figure: pedestals ----------------
        fig, ax = plt.subplots(1, 3, figsize=(19, 5.4))
        for th in have:
            _, m, w, ok = P[th]
            ax[0].plot(slabs, [np.median(m[ok & (slab == s)[:, None]]) for s in slabs], marker=MARK[th], color=COLOR[th],
                       lw=2, ms=7, label=f"th{th}")
            ax[2].plot(slabs, [np.median(w[ok & (slab == s)[:, None]]) for s in slabs], marker=MARK[th], color=COLOR[th],
                       lw=2, ms=7, label=f"th{th}")
        bins = np.linspace(-6, 6, 121)
        for a, b in pairs:
            ok = P[a][3] & P[b][3]
            d = (P[b][1] - P[a][1])[ok]
            med, hw = pct(d)
            ax[1].hist(np.clip(d, -6, 6), bins=bins, histtype="step", lw=2, color=PAIR_COLOR[(a, b)],
                       label=f"th{b} − th{a}: median {med:+.2f}, hw68 {hw:.2f}")
        ax[0].set_ylabel("Median pedestal mean [ADC]"); ax[2].set_ylabel("Median pedestal width (RMS) [ADC]")
        ax[1].set_xlabel("Pedestal difference per (channel, SCA) [ADC] (clipped at ±6)"); ax[1].set_yscale("log")
        ax[1].set_ylabel("cells")
        for a_ in (ax[0], ax[2]):
            a_.set_xlabel("Slab")
        for a_ in ax:
            a_.grid(alpha=0.25); a_.legend(fontsize=9)
        fig.suptitle(f"Pedestals across threshold sets, {gain} gain ({args.label})", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(os.path.join(args.out, f"thresholds_ped_{args.tag}_{gain}gain.png"), dpi=105); plt.close(fig)

        # ---------------- figure: MIPs ----------------
        fig, ax = plt.subplots(1, 3, figsize=(19, 5.4))
        allm = np.concatenate([M[th][1][M[th][4]] for th in have])
        bins = np.linspace(0, np.percentile(allm, 99.5), 100)
        for th in have:
            _, mpv, w, _, ok = M[th]
            ax[0].hist(mpv[ok], bins=bins, histtype="step", lw=2, color=COLOR[th],
                       label=f"th{th}: median {np.median(mpv[ok]):.2f} ADC")
            ax[2].plot(slabs, [np.median(w[ok & (M[th][0][:, 0] == s)]) if (ok & (M[th][0][:, 0] == s)).sum() else np.nan
                               for s in slabs], marker=MARK[th], color=COLOR[th], lw=2, ms=7,
                       label=f"th{th}: median {np.median(w[ok]):.2f}")
            ax[1].plot(slabs, [np.median(mpv[ok & (M[th][0][:, 0] == s)]) if (ok & (M[th][0][:, 0] == s)).sum() else np.nan
                               for s in slabs], marker=MARK[th], color=COLOR[th], lw=2, ms=7, label=f"th{th}")
        ax[0].set_xlabel("MIP MPV [ADC]"); ax[0].set_ylabel("channels")
        ax[1].set_xlabel("Slab"); ax[1].set_ylabel("Median MIP MPV [ADC]")
        ax[2].set_xlabel("Slab"); ax[2].set_ylabel("Median MIP Landau width [ADC]")
        for a_ in ax:
            a_.grid(alpha=0.25); a_.legend(fontsize=9)
        fig.suptitle(f"MIPs across threshold sets, {gain} gain ({args.label})", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.94])
        fig.savefig(os.path.join(args.out, f"thresholds_mip_{args.tag}_{gain}gain.png"), dpi=105); plt.close(fig)
    txt.close()


if __name__ == "__main__":
    main()
