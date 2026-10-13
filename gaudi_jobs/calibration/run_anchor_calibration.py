#!/usr/bin/env python
"""
Batch driver for the LG->HG anchor calibration (GainAnchorCalibrator).

Above the high-gain switch the event builder takes a hit's energy from the low gain through the anchor line
adc_low - ped_lg = k * (adc_high - ped_hg) + c. This measures (k, c) per data-taking period from the decoded
chunks of an electron run and writes the file the reconstruction reads:

    <calib-dir>/anchor/th<N>/gain_anchor_th<N>.txt          (k, c, switch, source_run, ...)
    <calib-dir>/anchor/th<N>/gain_anchor_fit_<run>_th<N>.root/.png   (with --diagnostics, the default)

gaudi_jobs/condor/generate_reco_dag.py --calib-dir <calib-dir> [--anchor-own-th] picks it up from there.

One k4run of gaudi_source/options/run_gain_anchor.py per period. The anchor needs showers (the fiducial band is
150-1200 high-gain ADC; a MIP is ~19), so it is measured on an electron run, never on the muon runs the
pedestal and MIP tables come from. The default run per period is the one the deployed anchors were fitted on:

    th210  TB2026CERN_eudaq_run_000166   (99 GeV)
    th220  TB2026CERN_run_000072
    th230  TB2026CERN_run_000012         (74 GeV)

Pedestals: <calib-dir>/pedestals/th<ped-th>/Pedestal_*_{high,low}gain.txt, ped-th = 210 by default because every
run of the final campaign is reconstructed with the th210 pedestals, and the anchor must use the pedestals the
reconstruction subtracts.

Examples::

    # the three periods, as deployed, into a scratch folder to compare
    python gaudi_jobs/calibration/run_anchor_calibration.py --all-periods --calib-dir /tmp/anchor_check

    # one period, traced on the occupancy of the saturated hits only
    python gaudi_jobs/calibration/run_anchor_calibration.py --th 210 --min-hitbits 12 --max-hitbits 24 \\
        --calib-dir /eos/experiment/drdcalo/siw-ecal/TB2026-06/Calib/MuonCalib --force

An existing anchor file is never overwritten without --force.
"""
import argparse
import glob
import os
import subprocess
import sys

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
_STEERING = os.path.join(_REPO, "gaudi_source", "options", "run_gain_anchor.py")
_DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"
# The official tables on EOS (settings.yml calib_dir()/MuonCalib).
sys.path.insert(0, _REPO)
from siwecal_common import paths  # noqa: E402
_DEFAULT_CALIB_DIR = os.path.join(paths.calib_dir(), "MuonCalib")

# Electron run each data-taking period's anchor is measured on (the ones the deployed anchors used).
PERIOD_RUNS = {
    "210": "TB2026CERN_eudaq_run_000166",
    "220": "TB2026CERN_run_000072",
    "230": "TB2026CERN_run_000012",
}


def pedestal_files(calib_dir, ped_th):
    """The (high, low) gain pedestal tables of <calib_dir>/pedestals/th<ped_th>/ -- exactly one of each."""
    d = os.path.join(calib_dir, "pedestals", f"th{ped_th}")
    out = []
    for gain in ("high", "low"):
        found = sorted(glob.glob(os.path.join(d, f"Pedestal_*_{gain}gain.txt")))
        if len(found) != 1:
            raise SystemExit(f"ERROR: expected one Pedestal_*_{gain}gain.txt in {d}, found {len(found)}: {found}")
        out.append(found[0])
    return out


def chunk_list(converted_dir, run, max_chunks):
    files = sorted(glob.glob(os.path.join(converted_dir, run, "chunks", "chunk_*.root")))
    if not files:
        raise SystemExit(f"ERROR: no decoded chunks in {os.path.join(converted_dir, run, 'chunks')}")
    return files[:max_chunks] if max_chunks > 0 else files


def plot(diag_root, png, th, run):
    """The fit as the reference script drew it: 2D histogram, slice medians, the line, the band and the switch."""
    import numpy as np
    import uproot
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    f = uproot.open(diag_root)
    counts, xe, ye = f["anchor_hist"].to_numpy()
    k = f["k"].member("fVal")
    c = f["c"].member("fVal")
    lo, hi, sw = (f[n].member("fVal") for n in ("fiducial_low", "fiducial_high", "switch"))
    gx, gy = f["slice_medians"].values()
    fx, fy = f["fit_points"].values()
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.pcolormesh(xe, ye, np.where(counts > 0, counts, np.nan).T, norm=LogNorm(), cmap="Greys", shading="flat")
    ax.axvspan(lo, hi, color="#1baf7a", alpha=0.10, lw=0, label=f"fiducial band {lo:.0f}–{hi:.0f} ADC")
    ax.axvline(sw, color="#eb6834", ls="--", lw=1.4, label=f"switch to the low gain, {sw:.0f} ADC")
    ax.plot(gx, gy, ".", color="#2a78d6", ms=4, label="median adc_low per adc_high slice")
    ax.plot(fx, fy, "o", color="#2a78d6", ms=4, mfc="none")
    xx = np.linspace(xe[0], xe[-1], 100)
    ax.plot(xx, k * xx + c, color="#1baf7a", lw=2, label=f"fit: k = {k:.4f}, c = {c:+.2f}   (1/k = {1 / k:.2f})")
    ax.set_xlim(-50, 2350)
    ax.set_ylim(-20, 300)
    ax.set_xlabel("adc_high − pedestal [ADC]")
    ax.set_ylabel("adc_low − pedestal [ADC]")
    ax.set_title(f"{run} (th{th}): LG→HG anchor, GainAnchorCalibrator")
    ax.legend(loc="upper left", fontsize=10)
    fig.tight_layout()
    fig.savefig(png, dpi=110)
    plt.close(fig)


def run_period(args, th):
    run = args.run or PERIOD_RUNS[th]
    out_dir = os.path.join(args.calib_dir, "anchor", f"th{th}")
    out = os.path.join(out_dir, f"gain_anchor_th{th}.txt")
    if os.path.exists(out) and not args.force:
        raise SystemExit(f"ERROR: {out} exists; pass --force to overwrite it, or another --calib-dir")
    os.makedirs(out_dir, exist_ok=True)
    ped_hg, ped_lg = pedestal_files(args.ped_dir or args.calib_dir, args.ped_th)
    chunks = chunk_list(args.converted_dir, run, args.max_chunks)
    tag = run.split("_")[-1]
    diag = os.path.join(out_dir, f"gain_anchor_fit_{tag}_th{th}.root") if args.diagnostics else ""
    env = {
        **os.environ,
        "ANCHOR_INPUT_FILES": ",".join(chunks),
        "ANCHOR_PEDESTAL_FILE": ped_hg,
        "ANCHOR_PEDESTAL_FILE_LOWGAIN": ped_lg,
        "ANCHOR_OUTPUT_FILE": out,
        "ANCHOR_DIAGNOSTICS_FILE": diag,
        "ANCHOR_MAX_ACQ": str(args.max_acq),
        "ANCHOR_FID_LO": str(args.fid_lo),
        "ANCHOR_FID_HI": str(args.fid_hi),
        "ANCHOR_MIN_HITBITS": str(args.min_hitbits),
        "ANCHOR_MAX_HITBITS": str(args.max_hitbits),
        "ANCHOR_SWITCH": str(args.switch),
        "ANCHOR_SOURCE": run,
        "ANCHOR_TH": th,
    }
    print(f"[anchor th{th}] {run}: {len(chunks)} chunks, pedestals th{args.ped_th} -> {out}", flush=True)
    res = subprocess.run(["k4run", _STEERING], env=env)
    if res.returncode != 0:
        raise SystemExit(f"ERROR: k4run failed for th{th} ({run})")
    vals = {ln.split()[0]: ln.split()[1] for ln in open(out) if ln.strip() and not ln.startswith("#")}
    print(f"[anchor th{th}] k = {vals['k']}, c = {vals['c']}  ({vals['n_hits']} fired channels, "
          f"{vals['n_slices']} slices)")
    if diag and not args.no_plot:
        png = os.path.splitext(diag)[0] + ".png"
        plot(diag, png, th, run)
        print(f"[anchor th{th}] plot: {png}")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    which = p.add_mutually_exclusive_group(required=True)
    which.add_argument("--th", choices=sorted(PERIOD_RUNS), help="One data-taking period")
    which.add_argument("--all-periods", action="store_true", help="th210, th220 and th230")
    p.add_argument("--run", default=None, help="Electron run to measure on (default: the period's, see above); "
                                                "only with --th")
    p.add_argument("--calib-dir", default=_DEFAULT_CALIB_DIR,
                   help=f"Calibration folder: anchor/th<N>/ is written here (default {_DEFAULT_CALIB_DIR})")
    p.add_argument("--ped-dir", default=None, help="Folder with pedestals/th<N>/ (default: --calib-dir)")
    p.add_argument("--ped-th", default="210",
                   help="Pedestal tables to subtract: those the reconstruction uses (default 210)")
    p.add_argument("--converted-dir", default=_DEFAULT_CONVERTED_DIR, help="Decoded chunks: <dir>/<run>/chunks/")
    p.add_argument("--max-chunks", type=int, default=60, help="First N chunks of the run (0 = all; default 60)")
    p.add_argument("--max-acq", type=int, default=8000, help="Acquisitions read (0 = all; default 8000)")
    p.add_argument("--fid-lo", type=float, default=150.0)
    p.add_argument("--fid-hi", type=float, default=1200.0)
    p.add_argument("--min-hitbits", type=int, default=0, help="Occupancy window: hit bits per chip-SCA (default 0)")
    p.add_argument("--max-hitbits", type=int, default=64, help="(default 64 = no cut)")
    p.add_argument("--switch", type=float, default=1500.0, help="Written as 'switch' (AdcSaturationThreshold)")
    p.add_argument("--no-diagnostics", dest="diagnostics", action="store_false",
                   help="Do not write the diagnostics ROOT file and plot")
    p.add_argument("--no-plot", action="store_true", help="Diagnostics ROOT file only, no PNG")
    p.add_argument("--force", action="store_true", help="Overwrite an existing anchor file")
    args = p.parse_args(argv)
    if args.run and args.all_periods:
        p.error("--run only with --th")
    for th in (sorted(PERIOD_RUNS) if args.all_periods else [args.th]):
        run_period(args, th)
    print("\n[Done]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
