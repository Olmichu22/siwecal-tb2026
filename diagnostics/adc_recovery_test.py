"""Can an ADC selection, corrected for the chip common mode, give back the hits the baseline drop cost?

Input: the same chunks rebuilt twice (diagnostics/out_cm/run_adcsel.sh): HitSelection=hitbit, and
HitSelection=adc with AdcHitThreshold=8 so that every candidate without a hit bit is stored. "Has a bit" is
hit_bit_window (a bit in any SCA of the window): in adc mode the channel is read from its largest-amplitude SCA, and
hit_bit only refers to that one.

Offline, per (event, slab, chip, SCA): charge Q = sum of hit_lg / k over its hits WITH a bit (k = 0.0924); common
mode CM = SLOPE * Q / 1000 (SLOPE = -0.8 ADC, diagnostics/cm_baseline_analysis.py); corrected amplitude
a' = hit_hg - CM. A hit is kept if it has a bit, or if a' > T (scan of T: the trigger threshold is 16.1 ADC at
th210). Then the 0.5 MIP cut on the corrected energy, unmasked hits, shower-like events (nhit > 0.5 x p90), as in
the report. Compared with the digitised simulation (final_v5) as it is.

    python diagnostics/adc_recovery_test.py <adcsel dir>
"""
import sys

import numpy as np
import uproot
import awkward as ak

SIM = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Simulation/Processed/adc_vs_tb/final_v5/trees"
K = 0.0924
SLOPE = -0.8
THRESHOLDS = [None, 30, 20, 16, 12]
GROUPS = [("291", "20", "26_-52"), ("286", "74", "22_-49")]
BR = ["hit_slab", "hit_chip", "hit_sca", "hit_hg", "hit_lg", "hit_energy", "hit_ismasked", "hit_bit_window"]


def sim_nhit(E, pos):
    a = uproot.open(f"{SIM}/ecal_th210_e{E}_{pos}.root")["ecal"].arrays(["hit_energy", "hit_ismasked"])
    n = ak.to_numpy(ak.num(a.hit_energy[(a.hit_ismasked == 0) & (a.hit_energy > 0.5)]))
    return n[n > 0.5 * np.percentile(n, 90)].mean()


def main():
    d = sys.argv[1]
    for run, E, pos in GROUPS:
        ns = sim_nhit(E, pos)
        hb = uproot.open(f"{d}/ecal_{run}_hitbit.root")["ecal"].arrays(["hit_energy", "hit_ismasked"])
        nh = ak.to_numpy(ak.num(hb.hit_energy[(hb.hit_ismasked == 0) & (hb.hit_energy > 0.5)]))
        nh_sel = nh[nh > 0.5 * np.percentile(nh, 90)].mean()
        a = uproot.open(f"{d}/ecal_{run}_adc.root")["ecal"].arrays(BR)
        print(f"== {E} GeV (eudaq {run}), {len(nh)} events; sim {ns:.1f} hits; hit-bit rebuild {nh_sel:.1f} "
              f"(sim/data {ns / nh_sel:.3f})")
        res = {}
        for T in THRESHOLDS:
            n_ev, rec_quiet, rec_busy, n_q, n_b = [], 0, 0, 0, 0
            for sl, ch, sc, hg, lg, e, m, bit in zip(a.hit_slab, a.hit_chip, a.hit_sca, a.hit_hg, a.hit_lg,
                                                     a.hit_energy, a.hit_ismasked, a.hit_bit_window):
                sl, ch, sc, hg, lg, e, m, bit = (ak.to_numpy(x) for x in (sl, ch, sc, hg, lg, e, m, bit))
                key = sl * 10000 + ch * 100 + sc
                u, inv = np.unique(key, return_inverse=True)
                q = np.bincount(inv, weights=np.where(bit == 1, lg / K, 0.0))[inv]
                nbit = np.bincount(inv, weights=(bit == 1).astype(float))[inv]
                if T is None:
                    cm = np.zeros_like(hg); keep = bit == 1
                else:
                    cm = SLOPE * q / 1000.0
                    keep = (bit == 1) | (hg - cm > T)
                scale = np.where((hg > 0) & (hg < 1500), (hg - cm) / np.maximum(hg, 1e-3), 1.0)
                ok = keep & (m == 0) & (e * scale > 0.5)
                n_ev.append(ok.sum())
                new = ok & (bit == 0)
                rec_quiet += (new & (nbit <= 3)).sum(); n_q += ((nbit >= 1) & (nbit <= 3)).sum() > 0
                rec_busy += (new & (nbit >= 15)).sum()
            n_ev = np.array(n_ev)
            sel = n_ev > 0.5 * np.percentile(n_ev, 90)
            res[T] = n_ev[sel].mean()
            lab = "hit bit only (adc rebuild)" if T is None else f"bit, or corrected ADC > {T} ADC"
            print(f"  {lab:<32} {res[T]:7.1f} hits   sim/data {ns / res[T]:.3f}   recovered per event: in chips "
                  f"<= 3 bits {rec_quiet / len(n_ev):5.1f}, >= 15 bits {rec_busy / len(n_ev):5.1f}")


if __name__ == "__main__":
    main()
