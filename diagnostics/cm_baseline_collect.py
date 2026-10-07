"""Collect the chip baseline (ADC - pedestal of the channels without a hit bit) per chip-SCA, for the common-mode study.

For every (acquisition, slab, chip, SCA) of the decoded chunks (Data/rundata_converted_fixed, fixed SCA pairing):

  nbits   channels with hitbit_high
  q       charge of the fired channels, from the LOW gain (does not saturate), in high-gain ADC:
          sum (adc_low - ped_lg) / k over the fired channels, k = 0.0924 (th210 LG->HG anchor)
  base_*  median of adc - pedestal over the NON-fired channels, high gain, split by distance to the nearest fired
          channel of the same chip: near (adjacent pad, 1 pitch incl. diagonals) and far (>= 2 pitches away)
  lg_far  the same 'far' median in the low gain, LG ADC

and per-channel residual histograms (high gain, 0.25 ADC bins) by nbits bin and near/far, so that the baseline can be
read from the MODE (the medians are pulled up by real signal below the trigger threshold).

Pedestals: MuonCalib_gaudi_fixed th210, both gains. Slab 12 (COB) uses its own pad map.

    python diagnostics/cm_baseline_collect.py <outdir> --run e286 [--chunks 10-12] [--max-acq N]
"""
import argparse
import os

import numpy as np
import uproot

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
CONV = f"{B}/Data/rundata_converted_fixed"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAL = f"{REPO}/calibration/MuonCalib_gaudi_fixed/pedestals/th210"
MAP = f"{REPO}/mappings/fev10_rotate_chip_channel_x_y_mapping.txt"
MAP_COB = f"{REPO}/mappings/fev11_cob_good_rotate_chip_channel_x_y_mapping.txt"
K_ANCHOR = 0.0924
NB_EDGES = np.array([0, 1, 2, 3, 5, 7, 10, 13, 16, 20, 25, 30, 35, 40, 50, 65])
HBINS = np.arange(-40, 60.001, 0.25)


def load_ped(path):
    ped = np.full((15, 16, 64, 15), np.nan, np.float32)
    for line in open(path):
        if line.startswith("#"):
            continue
        v = line.split()
        if len(v) < 48:
            continue
        m = np.array(v[3:48], float).reshape(15, 3)
        ok = (m[:, 0] > 0) & (m[:, 2] > 0)
        ped[int(v[0]), int(v[1]), int(v[2])] = np.where(ok, m[:, 0], np.nan)
    return np.moveaxis(ped, -1, -2)   # [layer, chip, sca, channel]


def adjacency(path):
    """A1[chip, ch, ch'] adjacent pads (<= 1 pitch, diagonals included); A2 within 2 pitches."""
    xy = np.full((16, 64, 2), np.nan)
    for line in open(path):
        v = line.split()
        if not v or not v[0].isdigit():
            continue
        xy[int(v[0]), int(v[3])] = float(v[4]), float(v[5])
    d = np.max(np.abs(xy[:, :, None, :] - xy[:, None, :, :]), axis=-1)   # Chebyshev distance [chip, 64, 64]
    pitch = 5.53
    a1 = (d < 1.5 * pitch) & (d > 0)
    a2 = (d < 2.5 * pitch) & (d > 0)
    return a1.astype(np.float32), a2.astype(np.float32)


def run_name(r):
    return f"TB2026CERN_eudaq_run_{int(r[1:]):06d}" if r.startswith("e") else f"TB2026CERN_run_{int(r):06d}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--run", required=True)
    ap.add_argument("--chunks", default="10-12")
    ap.add_argument("--max-acq", type=int, default=None)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    ped_hg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_highgain.txt")
    ped_lg = load_ped(f"{CAL}/Pedestal_TB2026CERN_run_000th210_lowgain.txt")
    A1, A2 = adjacency(MAP)
    A1c, A2c = adjacency(MAP_COB)
    c0, c1 = (int(x) for x in a.chunks.split("-"))
    files = [f"{CONV}/{run_name(a.run)}/chunks/chunk_{i:04d}.root" for i in range(c0, c1 + 1)]
    files = [f for f in files if os.path.exists(f)]
    if not files:   # LCIO-decoded runs come as a single chunk
        files = sorted(f for f in os.listdir(f"{CONV}/{run_name(a.run)}/chunks"))[:1]
        files = [f"{CONV}/{run_name(a.run)}/chunks/{f}" for f in files]
    nbb = len(NB_EDGES) - 1
    hist = {k: np.zeros((nbb, len(HBINS) - 1)) for k in ("all", "near", "far", "lg_far")}
    recs = []
    n_acq = 0
    for path in files:
        for d in uproot.open(path)["siwecaldecoded"].iterate(
                ["chipid", "nColumns", "adc_high", "adc_low", "hitbit_high"], step_size=20, library="np"):
            if a.max_acq and n_acq >= a.max_acq:
                break
            cid = d["chipid"]                                         # [a, slab, chip]
            ok_chip = cid >= 0
            cidc = np.clip(cid, 0, 15)
            L = np.arange(15)[None, :, None]
            P = ped_hg[L, cidc]                                        # [a, slab, chip, sca, ch]
            PL = ped_lg[L, cidc]
            adc = d["adc_high"].astype(np.float32)
            adl = d["adc_low"].astype(np.float32)
            hb = d["hitbit_high"].astype(np.float32)
            valid = (np.arange(15)[None, None, None, :] < d["nColumns"][..., None]) & ok_chip[..., None]
            good = valid[..., None] & np.isfinite(P) & np.isfinite(PL) & (adc > 0)
            r = adc - P
            rl = adl - PL
            # neighbourhood of the fired channels, per chip: hits @ adjacency (slab 12 with the COB map)
            Aa1 = np.where((np.arange(15) == 12)[None, :, None, None, None], A1c[cidc], A1[cidc])   # [a,slab,chip,64,64]
            Aa2 = np.where((np.arange(15) == 12)[None, :, None, None, None], A2c[cidc], A2[cidc])
            hbg = hb * good
            near1 = np.einsum("...sc,...cd->...sd", hbg, Aa1) > 0
            near2 = np.einsum("...sc,...cd->...sd", hbg, Aa2) > 0
            nob = good & (hb == 0)
            m_near = nob & near1
            m_far = nob & ~near2
            nbits = hbg.sum(-1)                                        # [a, slab, chip, sca]
            q = np.where(hbg > 0, rl, 0).sum(-1) / K_ANCHOR
            nan = np.float32(np.nan)
            base_all = np.nanmedian(np.where(nob, r, nan), axis=-1)
            base_near = np.nanmedian(np.where(m_near, r, nan), axis=-1)
            base_far = np.nanmedian(np.where(m_far, r, nan), axis=-1)
            lg_far = np.nanmedian(np.where(m_far, rl, nan), axis=-1)
            n_far = m_far.sum(-1)
            sel = valid & (good.sum(-1) > 8)
            ia, isl, ich, isc = np.nonzero(sel)
            recs.append(np.stack([np.full(len(ia), n_acq) + ia, isl, cid[ia, isl, ich], isc, nbits[sel], q[sel],
                                  base_all[sel], base_near[sel], base_far[sel], lg_far[sel], n_far[sel]],
                                 axis=1).astype(np.float32))
            kb = np.clip(np.searchsorted(NB_EDGES, nbits, side="right") - 1, 0, nbb - 1)
            kbb = np.broadcast_to(kb[..., None], r.shape)
            for name, m, x in (("all", nob, r), ("near", m_near, r), ("far", m_far, r), ("lg_far", m_far, rl / K_ANCHOR)):
                idx = kbb[m]
                xv = x[m]
                h2, _, _ = np.histogram2d(idx, xv, bins=[np.arange(nbb + 1) - 0.5, HBINS])
                hist[name] += h2
            n_acq += len(cid)
        print(f"[{a.run}] {os.path.basename(path)}: {n_acq} acquisitions so far", flush=True)
    rec = np.concatenate(recs)
    np.savez_compressed(f"{a.outdir}/cm_{a.run}.npz", rec=rec, nb_edges=NB_EDGES, hbins=HBINS, n_acq=n_acq,
                        cols=np.array(["acq", "slab", "chip", "sca", "nbits", "q", "base_all", "base_near", "base_far",
                                       "lg_far", "n_far"]), **{f"h_{k}": v for k, v in hist.items()})
    print(f"[{a.run}] {n_acq} acquisitions, {len(rec)} chip-SCAs -> {a.outdir}/cm_{a.run}.npz")


if __name__ == "__main__":
    main()
