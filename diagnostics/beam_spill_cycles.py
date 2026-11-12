"""Spill cycle of the SPS beam per run: how often a spill arrives, measured from the data.

The DAQ has no wall-clock time per acquisition, only acqNumber. A spill shows up as a train of consecutive
high-occupancy acquisitions (see plot_beam_spill_structure.py); the distance between the starts of two trains, in
acquisitions, is the spill period in acquisition units. It is converted to seconds with the run's mean time per
acquisition, total acquisitions / run duration, for runs whose duration the logbook (TB2026-Run_Summary.xlsx)
records. During TB2026 the SPS alternated between one spill per supercycle and three spills per supercycle
(logbook: "1 spill per cycle", "3 spills per cycle", "spill each 13 s", "3 spills in 40 s").

    python3 diagnostics/beam_spill_cycles.py survey <outdir> RUN [RUN ...]   [MAXCHUNKS=20]
        -> <outdir>/cycles.txt (per run: acquisitions read, trains, median period in acquisitions, total acquisitions)
    python3 diagnostics/beam_spill_cycles.py analyze <outdir> RUN_A RUN_B <run_durations.txt>
        -> spill_periods.txt, spill_cycles_2x2.png, spill_period_vs_run.png (see analyze())

Occupancy = sum of `nhits` over valid (bcid >= 0, badbcid == 0) SCAs of the acquisition, spliced by acqNumber
(an acquisition split across two chunks is written twice as complementary halves).
"""
import glob
import os
import sys

import numpy as np
import uproot

DATA = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"
OCC_CUT = float(os.environ.get("OCC_CUT", "200"))
MIN_TRAIN = int(os.environ.get("MIN_TRAIN", "10"))     # acquisitions above the cut for a real train
MERGE_GAP = int(os.environ.get("MERGE_GAP", "30"))     # gaps shorter than this inside a train are bridged


def run_dir(run):
    return f"{DATA}/TB2026CERN_run_{int(run):06d}"


def occupancy(run, maxchunks):
    paths = sorted(glob.glob(f"{run_dir(run)}/chunks/chunk_*.root"))
    if maxchunks:
        paths = paths[:maxchunks]
    occ = {}
    for p in paths:
        t = uproot.open(p)["siwecaldecoded"]
        a = t.arrays(["acqNumber", "n_slboards", "nhits", "bcid", "badbcid"], library="np")
        for acq, nsl, nh, bc, bad in zip(a["acqNumber"], a["n_slboards"], a["nhits"], a["bcid"], a["badbcid"]):
            good = (bad[:nsl] == 0) & (bc[:nsl] >= 0)
            occ[int(acq)] = occ.get(int(acq), 0) + int(nh[:nsl][good].sum())
    acqs = np.array(sorted(occ))
    return acqs, np.array([occ[k] for k in acqs], dtype=float), len(paths)


def total_acquisitions(run):
    paths = sorted(glob.glob(f"{run_dir(run)}/chunks/chunk_*.root"))
    first = uproot.open(paths[0])["siwecaldecoded"]["acqNumber"].array(library="np")
    last = uproot.open(paths[-1])["siwecaldecoded"]["acqNumber"].array(library="np")
    return int(last.max() - first.min() + 1), len(paths)


def trains(acqs, occ):
    """Start and end acqNumber of each spill train (blocks above the cut, short gaps bridged)."""
    on = acqs[occ > OCC_CUT]
    if len(on) == 0:
        return []
    out, s, e, n = [], on[0], on[0], 1
    for a in on[1:]:
        if a - e <= MERGE_GAP:
            e, n = a, n + 1
        else:
            if n >= MIN_TRAIN:
                out.append((s, e, n))
            s, e, n = a, a, 1
    if n >= MIN_TRAIN:
        out.append((s, e, n))
    return out


def survey(outdir, runs, maxchunks):
    os.makedirs(outdir, exist_ok=True)
    lines = ["# run chunks_read acq_read trains median_period_acq p25 p75 median_train_len_acq total_acq total_chunks"]
    for r in runs:
        try:
            acqs, occ, nch = occupancy(r, maxchunks)
            tot, totch = total_acquisitions(r)
        except Exception as e:  # noqa: BLE001
            print(f"run {r}: {e}")
            continue
        tr = trains(acqs, occ)
        starts = np.array([t[0] for t in tr])
        per = np.diff(starts) if len(starts) > 1 else np.array([np.nan])
        lens = np.array([t[2] for t in tr]) if tr else np.array([np.nan])
        lines.append(f"{r} {nch} {len(acqs)} {len(tr)} {np.median(per):.0f} {np.percentile(per, 25):.0f} "
                     f"{np.percentile(per, 75):.0f} {np.median(lens):.0f} {tot} {totch}")
        print(lines[-1], flush=True)
        np.savez(f"{outdir}/occ_{r}.npz", acqs=acqs, occ=occ, periods=per)
    with open(f"{outdir}/cycles.txt", "a") as f:
        f.write("\n".join(lines) + "\n")


def read_survey(outdir):
    rows = {}
    for line in open(f"{outdir}/cycles.txt"):
        f = line.split()
        if not f or f[0].startswith("#"):
            continue
        rows[int(f[0])] = dict(acq_read=int(f[2]), trains=int(f[3]), period=float(f[4]), total=int(f[8]))
    return rows


def read_durations(path):
    d = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        d[int(f[0])] = float(f[1])
    return d


def analyze(outdir, run_a, run_b, durations_file, dur_err_s=30.0):
    """Spill period per mode from many trains.

    Within a run the interval between consecutive train starts is measured in acquisitions; it is very stable
    (a few acquisitions wide), so the time scale dominates: s/acq = logged duration / total acquisitions, the
    duration recorded to the minute (+-dur_err_s, uncorrelated between runs). Per run: median of the core intervals
    (within 15 % of the run's median; a missed spill gives a multiple), stat. error = 1.2533 * rms / sqrt(n).
    Per mode: weighted mean of the runs with stat (+) duration errors, inflated by sqrt(chi2/ndf) when > 1.
    Figures: spill_cycles_2x2.png ((a), (b): spill trains against time for run_a and run_b; (c), (d): accumulated
    intervals of all trains of each mode with a Gaussian fit) and spill_period_vs_run.png (per-run period)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows, dur = read_survey(outdir), read_durations(durations_file)
    # each interval is converted to seconds with its run's time scale and assigned to a mode by its length, so a run
    # that changed mode (logbook: runs 100, 105, 115) contributes to both; a missed spill (2 x 13.3 s) falls outside
    modes = {"3 spills per supercycle": (10.0, 18.0), "1 spill per supercycle": (30.0, 50.0)}
    lines = ["# spill period per mode from the intervals between consecutive train starts",
             f"# OCC_CUT={OCC_CUT} MIN_TRAIN={MIN_TRAIN} MERGE_GAP={MERGE_GAP}; duration error {dur_err_s:.0f} s per run;"
             " mode windows " + ", ".join(f"{m}: {lo}-{hi} s" for m, (lo, hi) in modes.items()),
             "# run mode total_acq duration_s ms_per_acq n_intervals median_s stat_s dur_s"]
    pooled = {m: [] for m in modes}
    aligned = {m: [] for m in modes}
    acqmed = {}
    perrun = {m: [] for m in modes}
    other = []
    for r, row in sorted(rows.items()):
        if not os.path.exists(f"{outdir}/occ_{r}.npz"):
            continue
        if r not in dur:
            other.append((r, row["period"], row["trains"], False))
            continue
        per = np.load(f"{outdir}/occ_{r}.npz")["periods"]
        spa = dur[r] / row["total"]
        sec = per[np.isfinite(per)] * spa
        used = False
        for m, (lo, hi) in modes.items():
            x = sec[(sec > lo) & (sec < hi)]
            if len(x) < 3:
                continue
            med = np.median(x)
            stat = 1.2533 * np.std(x) / np.sqrt(len(x))
            sysd = med * dur_err_s / dur[r]
            pooled[m].append(x)
            aligned[m].append(x / med)
            acqmed.setdefault(r, {})[m] = med / spa
            perrun[m].append((r, med, stat, sysd, len(x)))
            used = True
            lines.append(f"{r} {m.split()[0]}/cycle {row['total']} {dur[r]:.0f} {1e3 * spa:.2f} {len(x)} "
                         f"{med:.2f} {stat:.3f} {sysd:.3f}")
        if not used:
            other.append((r, row["period"], row["trains"], True))
    res = {}
    m3, m1 = list(modes)
    # time-scale-free check, every run (logged duration or not): a run whose intervals form two clusters (split at
    # the geometric mean of the shortest and longest) with >= 3 intervals each gives the ratio of the two periods
    ratios = []
    for r in sorted(rows):
        f = f"{outdir}/occ_{r}.npz"
        if not os.path.exists(f):
            continue
        q = np.load(f)["periods"]
        q = q[np.isfinite(q) & (q > 0)]
        if len(q) < 6 or q.max() / q.min() < 2:
            continue
        cut = np.sqrt(q.min() * q.max())
        lo_, hi_ = q[q < cut], q[q >= cut]
        if len(lo_) >= 3 and len(hi_) >= 3:
            ratios.append((r, np.median(hi_) / np.median(lo_), len(lo_), len(hi_)))
    for m in modes:
        if not perrun[m]:
            continue
        pr = np.array([p[1] for p in perrun[m]])
        P, eP, rms = pr.mean(), pr.std(ddof=1) / np.sqrt(len(pr)), pr.std(ddof=1)
        x = np.concatenate(pooled[m])
        u = np.concatenate(aligned[m])
        res[m] = (P, eP, rms, x, u, pr)
        lines.append(f"MODE {m}: {len(pr)} runs, {len(x)} intervals; period {P:.2f} +- {eP:.2f} s (mean of the run "
                     f"medians +- rms/sqrt(n)); run-to-run rms {rms:.2f} s (time scale from the logged durations); "
                     "")
    lines += [f"RATIO run {r}: long / short interval cluster in acquisitions = {q:.3f} ({a} short, {b} long)"
              for r, q, a, b in ratios]
    lines.append("# runs not used (no logged duration, too few trains, or another period): run median_acq trains has_duration")
    lines += [f"{r} {p:.0f} {n} {int(h)}" for r, p, n, h in other]
    print("\n".join(l for l in lines if l.startswith("MODE")))

    # ---- per-run period against run number
    fig, ax = plt.subplots(figsize=(11, 3.8))
    for m, mk in zip(modes, ("o", "s")):
        if m not in res:
            continue
        rr = [p[0] for p in perrun[m]]
        ax.errorbar(rr, [p[1] for p in perrun[m]], yerr=np.hypot([p[2] for p in perrun[m]], [p[3] for p in perrun[m]]),
                    fmt=mk, ms=4, capsize=2, label=m)
        ax.axhline(res[m][0], color="0.5", lw=0.8, ls="--")
    ax.set_xlabel("run")
    ax.set_ylabel("spill period [s]")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{outdir}/spill_period_vs_run.png", dpi=150)
    plt.close(fig)

    # ---- 2 x 2: (a), (b) spill trains against time, one run per mode; (c), (d) accumulated intervals of every
    # train of every run of the mode, with a Gaussian fit (Poisson errors per bin, +-12 % around the period)
    from scipy.optimize import curve_fit
    fig, axs = plt.subplots(2, 2, figsize=(13, 9))
    for ax, r, lab in zip(axs[0], (run_a, run_b), "ab"):
        z = np.load(f"{outdir}/occ_{r}.npz")
        spa = dur[r] / rows[r]["total"]
        t = (z["acqs"] - z["acqs"][0]) * spa
        k = t < 200
        ax.plot(t[k], np.maximum(z["occ"][k], 0.5), ".", ms=1.5, color="#29a3dc")
        ax.axhline(OCC_CUT, color="#c0392b", ls=":", lw=1, label=f"train threshold ({OCC_CUT:.0f} hits)")
        ax.set_yscale("log")
        ax.set_xlabel("time [s]  (acquisition index × %.1f ms)" % (1e3 * spa))
        ax.set_ylabel("accepted hits per acquisition (0 drawn at 0.5)")
        m = next((k_ for k_, (lo, hi) in modes.items() if lo < rows[r]["period"] * spa < hi), "")
        ax.set_title(f"({lab}) run {r}, {m}: one spill every {rows[r]['period'] * spa:.1f} s", fontsize=10.5)
        ax.legend(fontsize=8, loc="center right")
    for ax, (m, lab) in zip(axs[1], zip(modes, "cd")):
        if m not in res:
            continue
        P, eP, rms, x, u, pr = res[m]
        lo, hi = 0.85 * P, 1.15 * P
        bins = np.linspace(lo, hi, 61)
        ax.hist(x, bins=bins, histtype="step", color="0.55", lw=1.2,
                label=f"raw: each run on its logged time scale (run rms {rms:.2f} s)")
        y = u * P                                   # every train of every run, aligned on its run's median
        h, e = np.histogram(y, bins=bins)
        ctr = 0.5 * (e[1:] + e[:-1])
        ok = h > 0
        g = lambda tt, A, mu, sg: A * np.exp(-0.5 * ((tt - mu) / sg) ** 2)  # noqa: E731
        p, cov = curve_fit(g, ctr[ok], h[ok], p0=[h.max(), P, 0.01 * P], bounds=([0, lo, 1e-4 * P], [np.inf, hi, 0.1 * P]),
                           sigma=np.sqrt(h[ok]), absolute_sigma=True, maxfev=20000)
        ep = np.sqrt(np.diag(cov))
        chi2f = float((((h[ok] - g(ctr[ok], *p)) / np.sqrt(h[ok])) ** 2).sum())
        ndff = int(ok.sum()) - 3
        ax.errorbar(ctr[ok], h[ok], yerr=np.sqrt(h[ok]), fmt="o", ms=3.5, color="#29a3dc", capsize=0,
                    label=f"{len(y)} trains of {len(pr)} runs, aligned on their run median")
        tt = np.linspace(lo, hi, 600)
        ax.plot(tt, g(tt, *p), color="#c0392b", lw=2,
                label=f"Gaussian fit: σ = {abs(p[2]):.3f} ± {ep[2]:.3f} s ({100 * abs(p[2]) / P:.2f} %), χ²/ndf = {chi2f:.0f}/{ndff}")
        ax.axvspan(P - eP, P + eP, color="#c0392b", alpha=0.15, label=f"period {P:.2f} ± {eP:.2f} s")
        ax.set_xlabel("interval between consecutive spill trains [s]")
        ax.set_ylabel(f"intervals / {e[1] - e[0]:.2f} s")
        ax.set_title(f"({lab}) {m}: accumulated spill periods", fontsize=10.5)
        ax.legend(fontsize=7.5, loc="upper left")
        ax.set_ylim(0, 1.45 * max(h.max(), 1))
        lines.append(f"FIT {m}: aligned intervals, Gaussian sigma {abs(p[2]):.3f} +- {ep[2]:.3f} s "
                     f"({100 * abs(p[2]) / P:.2f} %), chi2/ndf {chi2f:.1f}/{ndff}, {len(y)} intervals")
    fig.tight_layout()
    fig.savefig(f"{outdir}/spill_cycles_2x2.png", dpi=150)
    plt.close(fig)
    open(f"{outdir}/spill_periods.txt", "w").write("\n".join(lines) + "\n")
    print("\n".join(l for l in lines if l.startswith("FIT")))

if __name__ == "__main__":
    mode, outdir, *runs = sys.argv[1:]
    if mode == "survey":
        survey(outdir, [int(r) for r in runs], int(os.environ.get("MAXCHUNKS", "20")))
    elif mode == "analyze":
        analyze(outdir, int(runs[0]), int(runs[1]), runs[2])
    else:
        raise SystemExit(__doc__)
