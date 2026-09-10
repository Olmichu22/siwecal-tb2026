"""
Jitter vs retriggering: la poblacion de hueco=1, ?es la misma particula vista
con skew entre slabs, o es el SKIROC redisparando en el BCID siguiente?

Distingue mirando `badbcid` (que el clustering NO usa) de los SCA cuyo bcid
tiene un predecesor valido en bcid-1 dentro de la misma acquisition, frente a
los SCA aislados. Si los "con predecesor" estan sistematicamente flageados,
es retrigger y fusionarlos mete ruido en el evento; si no, es skew real.
"""
import glob
import sys

import numpy as np
import uproot

SKIP_BCID_START = 50
DROP_BCIDS = (0, 901)

chunk_dir = ("/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi/"
             "TB2026CERN_run_000044/chunks")
n_chunks = int(sys.argv[1]) if len(sys.argv) > 1 else 20
files = sorted(glob.glob(f"{chunk_dir}/chunk_*.root"))[:n_chunks]

from collections import Counter
bad_with_pred = Counter()   # SCA cuyo bcid tiene predecesor valido en bcid-1
bad_isolated = Counter()    # SCA sin predecesor en bcid-1
nhits_with_pred = []
nhits_isolated = []
# ?el slab del "segundo" BCID es el mismo que el del primero? (retrigger = mismo chip)
same_chip_pairs = 0
diff_chip_pairs = 0

for path in files:
    t = uproot.open(path)["siwecaldecoded"]
    bcid = t["bcid"].array(library="np")
    nhits = t["nhits"].array(library="np")
    chipid = t["chipid"].array(library="np")
    badbcid = t["badbcid"].array(library="np")

    valid = (nhits > 0) & (bcid >= SKIP_BCID_START) & (chipid[..., None] >= 0)
    for d in DROP_BCIDS:
        valid &= bcid != d

    for i in range(bcid.shape[0]):
        v = valid[i]
        if not v.any():
            continue
        b = bcid[i][v]
        bad = badbcid[i][v]
        nh = nhits[i][v]
        slab_idx, chip_idx, _ = np.nonzero(v)
        present = set(b.tolist())
        has_pred = np.array([(x - 1) in present for x in b])

        for flag, n in zip(bad[has_pred], nh[has_pred]):
            bad_with_pred[int(flag)] += 1
        nhits_with_pred.extend(nh[has_pred].tolist())
        for flag, n in zip(bad[~has_pred], nh[~has_pred]):
            bad_isolated[int(flag)] += 1
        nhits_isolated.extend(nh[~has_pred].tolist())

        # para los que tienen predecesor: ?estaba el predecesor en el mismo chip?
        for k in np.flatnonzero(has_pred):
            prev_same = ((b == b[k] - 1) & (slab_idx == slab_idx[k])
                         & (chip_idx == chip_idx[k])).any()
            if prev_same:
                same_chip_pairs += 1
            else:
                diff_chip_pairs += 1

n_pred = sum(bad_with_pred.values())
n_iso = sum(bad_isolated.values())
print(f"SCA validos con predecesor en bcid-1 : {n_pred:,}")
print(f"SCA validos aislados                 : {n_iso:,}\n")

print("badbcid   con-predecesor        aislados")
print("-" * 46)
for flag in sorted(set(bad_with_pred) | set(bad_isolated)):
    a, c = bad_with_pred[flag], bad_isolated[flag]
    print(f"{flag:7d} {a:10,} ({a/max(n_pred,1):6.1%}) {c:10,} ({c/max(n_iso,1):6.1%})")

print(f"\nnhits medio  con-predecesor: {np.mean(nhits_with_pred):.2f}   "
      f"aislados: {np.mean(nhits_isolated):.2f}")
print(f"\nEl predecesor en bcid-1 estaba en:")
print(f"  el MISMO chip  : {same_chip_pairs:,}  <- firma de retrigger")
print(f"  OTRO chip/slab : {diff_chip_pairs:,}  <- firma de skew entre slabs")
