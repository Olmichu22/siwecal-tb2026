"""Where does the data lose hits against the simulation as the beam energy rises?

The digitised simulation (final_v5) has 1.03x the data's hits at 20 GeV and 1.10x at 74 GeV (th210 scan,
shower-like events). If that is a readout problem -- SCAs, the BCID window, the chip-noise veto -- it should leave
fingerprints the simulation cannot have. This collects, per run of the th210 scan, data (Reconstruction) and
digitised simulation side by side:

  1. hits per event per layer and per ring around the shower axis (sim / data): a deficit in the dense core
     points at readout or saturation, one in the halo at the threshold or trigger model;
  2. data hits against the acquisition: n_events_acq, and the event's rank in its acquisition (by BCID) --
     if a chip's 15 SCAs fill up, the late events of an acquisition lose hits;
  3. which SCA the hits are read from;
  4. hits per (event, slab, chip): a ceiling in data the simulation does not have.

Selection as in the report: unmasked hits above 0.5 MIP, shower-like events (nhit > 0.5 x p90).

    python diagnostics/hit_deficit_vs_energy.py <outdir> [--max-events 40000]
"""
import argparse
import os

import numpy as np
import uproot
import awkward as ak

B = "/eos/experiment/drdcalo/siw-ecal/TB2026-06"
RECO = f"{B}/Reconstruction"
SIM = f"{B}/Simulation/Processed/adc_vs_tb/final_v5/trees"
SCAN = [("7.5", "e293", "17_-50"), ("10", "e292", "19_-51"), ("20", "e291", "26_-52"), ("34", "e288", "27_-53"),
        ("52", "e287", "24_-51"), ("74", "e286", "22_-49")]
HIT_CUT = 0.5
N_LAYERS = 15
RINGS = np.array([0, 5, 10, 15, 20, 30, 40, 60, 200.])
DATA_BR = ["spill", "bcid", "n_events_acq", "hit_slab", "hit_chip", "hit_sca", "hit_energy", "hit_x", "hit_y",
           "hit_ismasked"]
SIM_BR = ["hit_slab", "hit_chip", "hit_sca", "hit_energy", "hit_x", "hit_y", "hit_ismasked"]


def run_path(r):
    name = f"TB2026CERN_eudaq_run_{int(r[1:]):06d}"
    return f"{RECO}/{name}/ecal_{name}.root"


def load(path, branches, max_events):
    a = uproot.open(path)["ecal"].arrays(branches, entry_stop=max_events)
    keep = (a["hit_ismasked"] == 0) & (a["hit_energy"] > HIT_CUT)
    h = {b: a[b][keep] for b in branches if b.startswith("hit_")}
    nhit = ak.to_numpy(ak.num(h["hit_slab"]))
    sel = nhit > 0.5 * np.percentile(nhit, 90)
    ev = {b: ak.to_numpy(a[b])[sel] for b in branches if not b.startswith("hit_")}
    return {b: v[sel] for b, v in h.items()}, ev, nhit[sel]


def per_layer(h):
    sl = ak.to_numpy(ak.flatten(h["hit_slab"])).astype(int)
    return np.bincount(sl, minlength=N_LAYERS)[:N_LAYERS] / len(h["hit_slab"])


def per_ring(h):
    e = h["hit_energy"]
    bx = ak.sum(h["hit_x"] * e, axis=1) / ak.sum(e, axis=1)
    by = ak.sum(h["hit_y"] * e, axis=1) / ak.sum(e, axis=1)
    r = np.sqrt((h["hit_x"] - bx) ** 2 + (h["hit_y"] - by) ** 2)
    return np.histogram(ak.to_numpy(ak.flatten(r)), bins=RINGS)[0] / len(e)


def chip_mult(h):
    """Hits per (event, slab, chip), flattened."""
    key = h["hit_slab"] * 100 + h["hit_chip"]
    out = []
    for k in key:
        _, c = np.unique(ak.to_numpy(k), return_counts=True)
        out.append(c)
    return np.concatenate(out)


def rank_in_acq(ev):
    """0 for the first event (lowest BCID) of its acquisition, 1 for the second, ..."""
    order = np.lexsort((ev["bcid"], ev["spill"]))
    rank = np.zeros(len(order), int)
    sp = ev["spill"][order]
    start = np.r_[0, np.flatnonzero(np.diff(sp)) + 1]
    for s, e in zip(start, np.r_[start[1:], len(sp)]):
        rank[order[s:e]] = np.arange(e - s)
    return rank


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("outdir")
    ap.add_argument("--max-events", type=int, default=40000)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    res = {}
    for E, run, pos in SCAN:
        hd, evd, nd = load(run_path(run), DATA_BR, a.max_events)
        hs, _, ns = load(f"{SIM}/ecal_th210_e{E}_{pos}.root", SIM_BR, None)
        r = {"E": float(E), "nhit_data": nd.mean(), "nhit_sim": ns.mean(),
             "layer_data": per_layer(hd), "layer_sim": per_layer(hs),
             "ring_data": per_ring(hd), "ring_sim": per_ring(hs)}
        # acquisition dependence (data only)
        nacq = evd["n_events_acq"]
        rank = rank_in_acq(evd)
        r["by_nacq"] = {lab: (nd[m].mean(), int(m.sum())) for lab, m in
                        (("1", nacq == 1), ("2", nacq == 2), ("3-4", (nacq >= 3) & (nacq <= 4)),
                         ("5-9", (nacq >= 5) & (nacq <= 9)), (">=10", nacq >= 10)) if m.sum() > 50}
        r["by_rank"] = {lab: (nd[m].mean(), int(m.sum())) for lab, m in
                        (("0", rank == 0), ("1", rank == 1), ("2-4", (rank >= 2) & (rank <= 4)),
                         ("5-9", (rank >= 5) & (rank <= 9)), (">=10", rank >= 10)) if m.sum() > 50}
        # SCA of the hits, and how deep chips go
        sca_d = ak.to_numpy(ak.flatten(hd["hit_sca"])).astype(int)
        sca_s = ak.to_numpy(ak.flatten(hs["hit_sca"])).astype(int)
        r["sca_data"] = np.bincount(sca_d, minlength=15)[:15] / len(sca_d)
        r["sca_sim"] = np.bincount(sca_s, minlength=15)[:15] / max(len(sca_s), 1)
        r["chip_data"] = np.bincount(chip_mult(hd), minlength=65)[:65]
        r["chip_sim"] = np.bincount(chip_mult(hs), minlength=65)[:65]
        res[E] = r
        print(f"[{E} GeV] data {len(nd)} ev, nhit {nd.mean():.1f}; sim {len(ns)} ev, nhit {ns.mean():.1f}; "
              f"sim/data {ns.mean() / nd.mean():.3f}", flush=True)
    np.save(f"{a.outdir}/hit_deficit.npy", res, allow_pickle=True)

    with open(f"{a.outdir}/hit_deficit_summary.txt", "w") as f:
        def out(s=""):
            print(s); f.write(s + "\n")
        out("sim/data hits per event, shower-like events, th210 scan")
        for E, r in res.items():
            out(f"  {E:>5} GeV  data {r['nhit_data']:7.1f}  sim {r['nhit_sim']:7.1f}  ratio {r['nhit_sim'] / r['nhit_data']:.3f}")
        out("\nsim/data per layer")
        out("  E     " + "".join(f"{l:>6}" for l in range(N_LAYERS)))
        for E, r in res.items():
            out(f"  {E:>5} " + "".join(f"{x:6.2f}" for x in r["layer_sim"] / r["layer_data"]))
        out("\nsim/data per ring around the energy barycentre [mm]: " + " ".join(f"{int(lo)}-{int(hi)}" for lo, hi in zip(RINGS[:-1], RINGS[1:])))
        for E, r in res.items():
            out(f"  {E:>5} " + "".join(f"{x:7.2f}" for x in r["ring_sim"] / r["ring_data"]) +
                "   | data hits/ev " + " ".join(f"{x:.1f}" for x in r["ring_data"]))
        out("\ndata nhit vs n_events_acq (mean, N)")
        for E, r in res.items():
            out(f"  {E:>5} " + "  ".join(f"{k}: {v[0]:.1f} ({v[1]})" for k, v in r["by_nacq"].items()))
        out("\ndata nhit vs rank of the event in its acquisition (mean, N)")
        for E, r in res.items():
            out(f"  {E:>5} " + "  ".join(f"{k}: {v[0]:.1f} ({v[1]})" for k, v in r["by_rank"].items()))
        out("\nfraction of hits read from SCA 0..14 (data | sim)")
        for E, r in res.items():
            out(f"  {E:>5} " + " ".join(f"{x:.3f}" for x in r["sca_data"]) + "  |  " +
                " ".join(f"{x:.2f}" for x in r["sca_sim"][:3]))
        out("\nhits per (event, slab, chip): fraction of chip-entries with >= 20, 30, 40, 50 hits (data | sim) and max")
        for E, r in res.items():
            def fr(c):
                c = np.asarray(c, float); t = c.sum()
                return " ".join(f"{c[k:].sum() / t:.4f}" for k in (20, 30, 40, 50)) + f"  max {np.flatnonzero(c)[-1]}"
            out(f"  {E:>5} " + fr(r["chip_data"]) + "  |  " + fr(r["chip_sim"]))


if __name__ == "__main__":
    main()
