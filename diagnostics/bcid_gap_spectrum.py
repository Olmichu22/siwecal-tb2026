"""
Pasos 1 y 2 del procedimiento: distribucion de huecos entre BCIDs ocupados dentro
de una misma acquisition, PRE-clustering, con y sin enmascarado de retriggers.

Replica exactamente EventBuilder.h::collectValidBcids + mergeIntoWindows:
  * itera todos los SCA (0..14) de cada (slab, chip) con chipid >= 0
  * descarta nhits == 0  (MaxHitsPerSca esta deshabilitado: 1e18)
  * descarta bcid < 0, bcid < skipBcidStart (50), bcid en dropBcids {0, 901}
  * junta los BCID unicos de TODA la acquisition (no por chip)
  * MergeDelta cierra ventana cuando el hueco al siguiente es >= MergeDelta

Usa el `bcid` CRUDO, no `corrected_bcid` (nota 1 de la cabecera de EventBuilder.h).

Sobre esa base corre TRES modos de enmascarado en una sola pasada (misma
estadistica exacta para los tres, y una sola lectura de los ficheros):

  none      lo que hace el builder hoy: nada. Linea base del paso 1.
  follower  el corte de la referencia (bcid_handling.py:38-41,69): dentro de la
            memoria de cada chip, un SCA cuyo bcid dista <= DROP_RETRIGGER_DELTA
            del SCA anterior es retrigger -> se enmascara SOLO EL SEGUIDOR; el
            SCA lider se conserva.
  chain     enmascara todo SCA con badbcid == 3. NO es un candidato a
            implementar: SlbFrameDecoder.h::tagBadBcid marca con 3 la cadena
            entera, lider incluido, asi que este modo borra tambien la senal
            fisica que provoco el retrigger. Esta aqui para MEDIR esa
            diferencia: (chain - follower) es lo que costaria la version
            ingenua.

El orden de las mascaras replica _get_bcids: primero nhits==0 -> bad_value,
despues el corte de retrigger, y solo despues skip_bcid_start / drop_values.

Modelo nulo: para cada acquisition, N BCIDs uniformes en el rango permitido
(N = el de ese modo). Reproduce los huecos accidentales de particulas
independientes; el exceso del observado sobre el nulo a huecos pequenos es el
jitter intra-particula.
"""
import glob
import sys

import numpy as np
import uproot

SKIP_BCID_START = 50
DROP_BCIDS = (0, 901)
BCID_OVERFLOW = 4096
BAD_VALUE = -999
DROP_RETRIGGER_DELTA = 2   # config_run.cfg:48
BADBCID_RETRIGGER = 3      # SlbFrameDecoder.h:418-426
MAX_GAP = 40               # eje del histograma
MERGE_SCAN = range(1, 9)
MODES = ("none", "follower", "chain")

chunk_dir = sys.argv[1] if len(sys.argv) > 1 else (
    "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi/"
    "TB2026CERN_run_000044/chunks")
n_chunks = int(sys.argv[2]) if len(sys.argv) > 2 else 20

files = sorted(glob.glob(f"{chunk_dir}/chunk_*.root"))[:n_chunks]
print(f"[in] {len(files)} chunk(s) de {chunk_dir}")
print(f"[in] modos: {', '.join(MODES)}  (delta retrigger = {DROP_RETRIGGER_DELTA})\n")

rng = np.random.default_rng(12345)
# Valores de BCID permitidos por las mascaras, para el modelo nulo.
allowed = np.array([b for b in range(SKIP_BCID_START, BCID_OVERFLOW)
                    if b not in DROP_BCIDS], dtype=np.int64)

gaps_obs = {m: np.zeros(MAX_GAP + 2, dtype=np.int64) for m in MODES}
gaps_null = {m: np.zeros(MAX_GAP + 2, dtype=np.int64) for m in MODES}
n_unique_all = {m: [] for m in MODES}
n_acq_empty = {m: 0 for m in MODES}
n_sca_valid = {m: 0 for m in MODES}     # SCA que sobreviven a todas las mascaras
n_hits_valid = {m: 0 for m in MODES}    # suma de nhits sobre esos SCA
n_windows = {m: {d: 0 for d in MERGE_SCAN} for m in MODES}
max_span = {m: {d: 0 for d in MERGE_SCAN} for m in MODES}
span_sum = {m: {d: 0 for d in MERGE_SCAN} for m in MODES}
n_acq = 0


def accumulate(gaps, hist):
    if gaps.size == 0:
        return
    clipped = np.minimum(gaps, MAX_GAP + 1)
    hist += np.bincount(clipped, minlength=MAX_GAP + 2)[:MAX_GAP + 2]


def follower_mask(bcid, occupied):
    """Port de BCIDHandler._is_retrigger (bcid_handling.py:38-41).

    Se aplica sobre la matriz de bcid POR CHIP con los SCA vacios ya puestos a
    bad_value, comparando slots de SCA consecutivos (no bcids ordenados). Un
    bad_value en cualquiera de los dos extremos da un delta fuera de rango, asi
    que nunca cuenta como retrigger. Marca solo la columna [1:], el seguidor.
    """
    m = np.where(occupied, bcid, BAD_VALUE)
    delta = m[..., 1:] - m[..., :-1]
    mask = np.zeros(bcid.shape, dtype=bool)
    mask[..., 1:] = (delta > 0) & (delta <= DROP_RETRIGGER_DELTA)
    return mask


for path in files:
    tree = uproot.open(path)["siwecaldecoded"]
    bcid = tree["bcid"].array(library="np")        # (n, 15, 16, 15)
    nhits = tree["nhits"].array(library="np")
    chipid = tree["chipid"].array(library="np")    # (n, 15, 16)
    badbcid = tree["badbcid"].array(library="np")

    # Orden de _get_bcids: SCA ocupados -> corte de retrigger -> cortes de bcid.
    occupied = (nhits > 0) & (chipid[..., None] >= 0)
    retrigger = {
        "none": np.zeros(bcid.shape, dtype=bool),
        "follower": follower_mask(bcid, occupied),
        "chain": badbcid == BADBCID_RETRIGGER,
    }
    keep_bcid = bcid >= SKIP_BCID_START
    for drop in DROP_BCIDS:
        keep_bcid &= bcid != drop

    n_acq += bcid.shape[0]
    for mode in MODES:
        valid = occupied & keep_bcid & ~retrigger[mode]
        n_sca_valid[mode] += int(valid.sum())
        n_hits_valid[mode] += int(nhits[valid].sum())

        for i in range(bcid.shape[0]):
            vals = np.unique(bcid[i][valid[i]])
            if vals.size == 0:
                n_acq_empty[mode] += 1
                continue
            n_unique_all[mode].append(vals.size)
            gaps = np.diff(vals)
            accumulate(gaps, gaps_obs[mode])

            # ventanas por MergeDelta: un hueco >= delta cierra ventana
            for d in MERGE_SCAN:
                breaks = np.flatnonzero(gaps >= d)
                n_windows[mode][d] += breaks.size + 1
                edges = np.concatenate(([0], breaks + 1, [vals.size]))
                spans = vals[edges[1:] - 1] - vals[edges[:-1]]
                max_span[mode][d] = max(max_span[mode][d], int(spans.max()))
                span_sum[mode][d] += int(spans.sum())

            # modelo nulo: mismos N, posiciones uniformes
            draw = np.unique(rng.choice(allowed, size=vals.size, replace=False))
            accumulate(np.diff(draw), gaps_null[mode])

n_unique_all = {m: np.asarray(v) for m, v in n_unique_all.items()}
tot_obs = {m: gaps_obs[m].sum() for m in MODES}
tot_null = {m: gaps_null[m].sum() for m in MODES}
n_ok = {m: n_acq - n_acq_empty[m] for m in MODES}

print(f"[acq] {n_acq:,} acquisitions\n")
print("[mask] efecto de cada mascara sobre la entrada del clustering")
print(" modo         acq vacias      SCA validos          hits (sum nhits)   BCID/acq")
print(" " + "-" * 80)
base_sca, base_hits = n_sca_valid["none"], n_hits_valid["none"]
for m in MODES:
    d_sca = 100.0 * (n_sca_valid[m] - base_sca) / base_sca
    d_hit = 100.0 * (n_hits_valid[m] - base_hits) / base_hits
    print(f" {m:10s} {n_acq_empty[m]:10,} {n_sca_valid[m]:14,} ({d_sca:+6.1f}%) "
          f"{n_hits_valid[m]:14,} ({d_hit:+6.1f}%) {n_unique_all[m].mean():9.2f}")
print("\n Los hits son la suma de nhits de los SCA supervivientes: la caida es COTA")
print(" SUPERIOR de la perdida real, porque _best_sca_per_channel ya deduplica.")
print(f" chain - follower = "
      f"{100.0 * (n_hits_valid['chain'] - n_hits_valid['follower']) / base_hits:+.1f}% "
      "de los hits <- lo que costaria enmascarar tambien el lider.\n")

print("[gap] espectro de huecos: observado y exceso sobre el nulo")
print(" hueco" + "".join(f"{m:>26s}" for m in MODES))
print("      " + "".join(f"{'obs':>11s}{'frac':>7s}{'exc/obs':>8s}" for _ in MODES))
print(" " + "-" * (5 + 26 * len(MODES)))
for g in range(1, min(13, MAX_GAP + 1)):
    row = f" {g:4d}"
    for m in MODES:
        o, nul = gaps_obs[m][g], gaps_null[m][g]
        exc = o - nul * (tot_obs[m] / tot_null[m]) if tot_null[m] else 0.0
        row += f"{o:11,}{o / tot_obs[m]:7.3f}{(exc / o if o else 0.0):8.2f}"
    print(row)
print(" huecos totales: " + ", ".join(f"{m}={tot_obs[m]:,}" for m in MODES))

print("\n[scan] ventanas por acquisition (pre-MinSlabsHit) y span medio")
print(" delta" + "".join(f"{m:>22s}" for m in MODES))
print("      " + "".join(f"{'vent/acq':>11s}{'span':>11s}" for _ in MODES))
print(" " + "-" * (5 + 22 * len(MODES)))
for d in MERGE_SCAN:
    row = f" {d:4d}"
    for m in MODES:
        w = n_windows[m][d]
        row += f"{w / n_ok[m]:11.2f}{span_sum[m][d] / w:11.2f}"
    print(row)
