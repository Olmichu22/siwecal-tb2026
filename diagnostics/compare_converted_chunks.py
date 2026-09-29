#!/usr/bin/env python
"""
Compare two decodings of the same raw data, chunk by chunk (e.g. the legacy
Data/rundata_converted_gaudi against Data/rundata_converted_fixed, re-decoded
after the SCA column-order fix).

The fix only re-pairs a chip's SCA data blocks with its BCIDs, so per chunk:

  * the chunk list and the entry count must be identical;
  * the acqNumber sequence must be identical;
  * hitbit_high == hitbit_low must still hold everywhere (both gains are
    reversed together; a mismatch means one gain was re-paired and not the other).

Also counts the chips whose data really moved (nColumns >= 2 and some
adc_high differs), as a sanity check that the fixed decoding is not the
legacy one copied.

    python diagnostics/compare_converted_chunks.py --runs TB2026CERN_run_000012,... \
        [--legacy DIR] [--fixed DIR] [--max-chunks N] [--out report.txt]
"""
import argparse
import glob
import os
import sys

import numpy as np
import uproot

BASE = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data"


def compare_chunk(fa, fb):
    ta, tb = uproot.open(fa)["siwecaldecoded"], uproot.open(fb)["siwecaldecoded"]
    res = {"entries_a": ta.num_entries, "entries_b": tb.num_entries, "acq_equal": False,
           "hitbit_hl_mismatch": 0, "multi_col_chips": 0, "moved_chips": 0}
    if ta.num_entries != tb.num_entries:
        return res
    br = ["acqNumber", "nColumns", "adc_high"]
    a = ta.arrays(br, library="np")
    b = tb.arrays(br + ["hitbit_high", "hitbit_low"], library="np")
    res["acq_equal"] = bool(np.array_equal(a["acqNumber"], b["acqNumber"]))
    for e in range(len(b["acqNumber"])):
        res["hitbit_hl_mismatch"] += int(np.count_nonzero(b["hitbit_high"][e] != b["hitbit_low"][e]))
        multi = b["nColumns"][e] >= 2
        res["multi_col_chips"] += int(multi.sum())
        if multi.any():
            diff = (a["adc_high"][e] != b["adc_high"][e]).any(axis=(-1, -2)) if a["adc_high"][e].ndim == 4 \
                else (a["adc_high"][e] != b["adc_high"][e]).reshape(multi.shape + (-1,)).any(axis=-1)
            res["moved_chips"] += int((diff & multi).sum())
    return res


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--runs", required=True, help="Comma-separated run names, or 'all' (every run in --fixed)")
    p.add_argument("--legacy", default=os.path.join(BASE, "rundata_converted_gaudi"))
    p.add_argument("--fixed", default=os.path.join(BASE, "rundata_converted_fixed"))
    p.add_argument("--max-chunks", type=int, default=0,
                   help="Content check (acqNumber, hit bits) on at most N chunks per run, spread over the run "
                        "(0 = all). Chunk list and entry counts are always checked on every chunk.")
    p.add_argument("--out", default=None)
    args = p.parse_args(argv)

    runs = sorted(os.listdir(args.fixed)) if args.runs == "all" else args.runs.split(",")
    lines, bad = [], 0
    for run in runs:
        la = sorted(glob.glob(os.path.join(args.legacy, run, "chunks", "chunk_*.root")))
        lb = sorted(glob.glob(os.path.join(args.fixed, run, "chunks", "chunk_*.root")))
        na, nb = [os.path.basename(f) for f in la], [os.path.basename(f) for f in lb]
        if not lb:
            lines.append(f"{run}: NO fixed chunks"); bad += 1; continue
        if not la:
            lines.append(f"{run}: {len(lb)} fixed chunks, no legacy decoding to compare"); continue
        if na != nb:
            lines.append(f"{run}: chunk lists differ (legacy {len(na)}, fixed {len(nb)})"); bad += 1; continue
        ent_bad = [n for fa, fb, n in zip(la, lb, na)
                   if uproot.open(fa)["siwecaldecoded"].num_entries != uproot.open(fb)["siwecaldecoded"].num_entries]
        idx = range(len(la)) if not args.max_chunks or len(la) <= args.max_chunks else \
            np.unique(np.linspace(0, len(la) - 1, args.max_chunks).astype(int))
        tot = {"acq_bad": 0, "hitbit_hl_mismatch": 0, "multi_col_chips": 0, "moved_chips": 0}
        for i in idx:
            r = compare_chunk(la[i], lb[i])
            tot["acq_bad"] += int(not r["acq_equal"])
            for k in ("hitbit_hl_mismatch", "multi_col_chips", "moved_chips"):
                tot[k] += r[k]
        ok = not ent_bad and tot["acq_bad"] == 0 and tot["hitbit_hl_mismatch"] == 0
        bad += int(not ok)
        lines.append(f"{run}: {'OK ' if ok else 'BAD'} chunks {len(lb)}, entry mismatches {len(ent_bad)}, "
                     f"acqNumber mismatches {tot['acq_bad']}/{len(idx)}, hitbit_high!=low {tot['hitbit_hl_mismatch']}, "
                     f"multi-column chips {tot['multi_col_chips']}, re-paired {tot['moved_chips']}")
        print(lines[-1], flush=True)
    lines.append(f"\n{len(runs)} run(s), {bad} with a problem")
    print(lines[-1])
    if args.out:
        with open(args.out, "w") as f:
            f.write("\n".join(lines) + "\n")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
