"""The LG->HG anchor against chip occupancy, for the common-mode study.

The event builder converts a hit above the high-gain switch (1500 ADC pedestal-subtracted) to high-gain units with
the anchor line adc_low = k * adc_high + c (k = 0.0924, c = 2.35 for th210), traced on quiet data. If the chip's
baseline drops differently in the two gains when the chip is busy, c moves with occupancy, and the energy of the
saturated hits of busy chips -- most of the saturated hits -- is biased by delta_c / k (x10.8).

Collects every fired channel (hitbit_high) with an unsaturated high gain: (adc_high - ped_hg, adc_low - ped_lg,
hit bits in the chip-SCA), th210 pedestals of MuonCalib_gaudi_fixed. Also the per-chip-SCA far baselines of both
gains are in cm_baseline_collect.py.

    python diagnostics/cm_anchor_collect.py <outdir> --run e286 [--chunks 10-15] [--max-acq 1500]
"""
import argparse
import os

import numpy as np
import uproot

from cm_baseline_collect import CONV, CAL, load_ped, run_name


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--run", required=True)
    ap.add_argument("--chunks", default="10-15")
    ap.add_argument("--max-acq", type=int, default=1500)
    a = ap.parse_args()
    ped_hg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_highgain.txt")
    ped_lg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_lowgain.txt")
    c0, c1 = (int(x) for x in a.chunks.split("-"))
    files = [f"{CONV}/{run_name(a.run)}/chunks/chunk_{i:04d}.root" for i in range(c0, c1 + 1)]
    files = [f for f in files if os.path.exists(f)] or [f"{CONV}/{run_name(a.run)}/chunks/chunk_0000.root"]
    out, n_acq = [], 0
    for path in files:
        for d in uproot.open(path)["siwecaldecoded"].iterate(
                ["chipid", "nColumns", "adc_high", "adc_low", "hitbit_high"], step_size=20, library="np"):
            if n_acq >= a.max_acq:
                break
            cid = d["chipid"]; cidc = np.clip(cid, 0, 15); L = np.arange(15)[None, :, None]
            rh = d["adc_high"].astype(np.float32) - ped_hg[L, cidc]
            rl = d["adc_low"].astype(np.float32) - ped_lg[L, cidc]
            hb = d["hitbit_high"]
            valid = (np.arange(15)[None, None, None, :] < d["nColumns"][..., None]) & (cid >= 0)[..., None]
            nbits = (hb * valid[..., None]).sum(-1)
            m = valid[..., None] & (hb == 1) & np.isfinite(rh) & np.isfinite(rl) & (rh > 20) & (rh < 1500)
            nb = np.broadcast_to(nbits[..., None], rh.shape)
            sl = np.broadcast_to(np.arange(15)[None, :, None, None, None], rh.shape)
            out.append(np.stack([rh[m], rl[m], nb[m], sl[m]], 1).astype(np.float32))
            n_acq += len(cid)
        print(f"[{a.run}] {os.path.basename(path)}: {n_acq} acquisitions", flush=True)
    arr = np.concatenate(out)
    np.save(f"{a.outdir}/anchor_{a.run}.npy", arr)
    print(f"[{a.run}] {len(arr)} unsaturated hits -> {a.outdir}/anchor_{a.run}.npy")


if __name__ == "__main__":
    main()
