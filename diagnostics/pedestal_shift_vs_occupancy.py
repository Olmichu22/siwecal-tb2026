"""Does a SKIROC chip's baseline move when many of its channels fire?

For every (acquisition, slab, chip, SCA) of the decoded chunks: the number of channels with a hit bit, and the
high-gain ADC of the channels WITHOUT one minus their pedestal (th210 table of MuonCalib_gaudi_fixed, same SCA).
With no common-mode effect the residual is ~0 at every occupancy. A residual that goes negative with the number
of fired channels is a baseline drop: in a busy chip every amplitude reads low, so hits lose energy and the
smallest fall below the trigger threshold -- which the simulation (no such effect) does not have.

The baseline is read from the MODE of the residual distribution (peak, +/-1.5 ADC centroid): the median and the
mean are pulled up by real signals below the trigger threshold, more so the busier the chip.

    python diagnostics/pedestal_shift_vs_occupancy.py <outdir> [--runs e286,e291] [--chunks 10-15]
"""
import argparse
import os

import numpy as np
import uproot

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
CONV = f"{B}/Data/rundata_converted_fixed"
PED = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "calibration/MuonCalib_gaudi_fixed/pedestals/th210/Pedestal_TB2026CERN_run_000th210_highgain.txt")
OCC_BINS = [0, 1, 3, 6, 10, 15, 20, 30, 40, 65]


def load_ped(path):
    """ped[layer, chip, channel, sca] (nan where not fitted)."""
    ped = np.full((15, 16, 64, 15), np.nan)
    for line in open(path):
        if line.startswith("#"):
            continue
        v = line.split()
        if len(v) < 3 + 45:
            continue
        l, c, ch = int(v[0]), int(v[1]), int(v[2])
        m = np.array(v[3:3 + 45], float).reshape(15, 3)
        ok = (m[:, 0] > 0) & (m[:, 2] > 0)
        ped[l, c, ch] = np.where(ok, m[:, 0], np.nan)
    return ped


def run_name(r):
    return f"TB2026CERN_eudaq_run_{int(r[1:]):06d}" if r.startswith("e") else f"TB2026CERN_run_{int(r):06d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--runs", default="e286,e291")
    ap.add_argument("--chunks", default="10-15")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ped = load_ped(PED)
    c0, c1 = (int(x) for x in a.chunks.split("-"))
    lines = []
    for r in a.runs.split(","):
        sums = {k: [] for k in range(len(OCC_BINS) - 1)}
        for ci in range(c0, c1 + 1):
            p = f"{CONV}/{run_name(r)}/chunks/chunk_{ci:04d}.root"
            if not os.path.exists(p):
                continue
            t = uproot.open(p)["siwecaldecoded"]
            d = t.arrays(["chipid", "nColumns", "adc_high", "hitbit_high"], library="np")
            adc, hb = d["adc_high"], d["hitbit_high"]           # [acq, slab, chip, sca, ch]
            chipid, ncol = d["chipid"], d["nColumns"]           # [acq, slab, chip]
            nacq = adc.shape[0]
            sca = np.arange(15)
            valid = sca[None, None, None, :] < ncol[..., None]  # [acq, slab, chip, sca]
            nbits = hb.sum(axis=-1)                              # [acq, slab, chip, sca]
            # pedestal per (acq, slab, chip, sca, ch) through the chip id of each slot
            cid = np.clip(chipid, 0, 15)
            P = ped[np.arange(15)[None, :, None], cid]          # [acq, slab, chip, ch, sca]
            P = np.moveaxis(P, -1, -2)                           # [acq, slab, chip, sca, ch]
            res = adc - P
            use = valid[..., None] & (hb == 0) & np.isfinite(P) & (chipid[..., None, None] >= 0) & (adc > 0)
            occ = np.broadcast_to(nbits[..., None], res.shape)
            for k, (lo, hi) in enumerate(zip(OCC_BINS[:-1], OCC_BINS[1:])):
                m = use & (occ >= lo) & (occ < hi)
                if m.any():
                    sums[k].append(res[m].astype(np.float32))
            print(f"[{r}] chunk {ci}: {nacq} acquisitions", flush=True)
        lines.append(f"== {r}: non-hit channels, adc_high - pedestal [ADC], by the number of hit bits in the same "
                     f"chip and SCA")
        for k, (lo, hi) in enumerate(zip(OCC_BINS[:-1], OCC_BINS[1:])):
            if not sums[k]:
                continue
            x = np.concatenate(sums[k])
            # the mode is the baseline itself; the median and the mean are pulled up by real sub-threshold signal
            h, e = np.histogram(x, bins=np.arange(-30.25, 30.26, 0.5))
            i = np.argmax(h)
            ctr = 0.5 * (e[1:] + e[:-1])
            w = slice(max(i - 3, 0), i + 4)
            mode = float(np.sum(ctr[w] * h[w]) / max(h[w].sum(), 1))
            lines.append(f"  hit bits {lo:2d}-{hi - 1:2d}: mode {mode:+6.2f}  median {np.median(x):+6.2f}  "
                         f"mean {x.mean():+6.2f}  (N = {len(x)})")
    txt = "\n".join(lines)
    print(txt)
    open(f"{a.outdir}/pedestal_shift_summary.txt", "w").write(txt + "\n")


if __name__ == "__main__":
    main()
