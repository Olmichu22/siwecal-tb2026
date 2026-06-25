# k4SiWEcalReco

Gaudi/k4FWCore stage that turns the per-event reconstructed `ecal` tree (written
by `siwecal_eventbuilder`) into the high-level reconstruction product, computing
the particle-discrimination shower variables in C++.

It is also **the single place where the cut/cleaned event collections are
generated**: it applies an optional hit-level MIP cut, then the event selection
(the same `CutSet` as `siwecal_validation`), and writes the cut-passing events in
either of two formats — an **EDM4hep** PID file (`ecal_<run>.edm4hep.root`) or a
plain **valtree** TTree (`ecal_<run>.valtree.root`). `siwecal_validation` no
longer generates any tree; it only reads these outputs and makes plots.

```
ecal_<run>.root  ──►  EcalToEDM4hep  ──►  CalorimeterHitCollection (+ EventHeader)
   (1 row/event)      (hit-MIP cut)               │
                                                  ▼
                                          EcalPidTransformer   (shower variables)
                                                  │
                                ClusterCollection (1 Cluster/event)
                                                  │
                      run_pid_batch.py: event cuts (CutSet) + format select
                            ├──►  ecal_<run>.edm4hep.root   (filtered EDM4hep)
                            └──►  ecal_<run>.valtree.root   (plain valcache-schema TTree)
```

The Gaudi chain itself is a strict **1→1** transform (the `ecal` tree is already
one row per physics event — no BCID fan-out): it writes a *full* EDM4hep file to
a temporary path. `run_pid_batch.py` then applies the event cuts and writes the
requested final, cut-passing output(s), deleting the temporary file.

## Components

- **`EcalToEDM4hep`** (`src/components/EcalToEDM4hep.cpp`) — reads the (non-podio)
  `ecal` tree with `TTreeReader` and emits, per event:
  - `CalorimeterHitCollection ECalHits` — `energy = hit_energy`,
    `position = (hit_x, hit_y, hit_z)`, `cellID` encoding `slab/chip/channel/sca`,
    and `type = slab` so the layer is readable natively (no bitfield decode).
  - `EventHeaderCollection EventHeader` — run, event (+ bcid in `timeStamp`,
    spill in `weight`).
  - Parallel `UserDataCollection`s (same order as the hits) for the per-hit
    quantities `CalorimeterHit` has no native field for:
    `ECalHitChip/Chan/Sca` (int) and `ECalHitHG/LG` (float).

  Channels flagged `hit_ismasked` in the `ecal` tree (no MIP calibration) are
  dropped from `ECalHits` and every parallel collection at once, so the output
  carries only calibrated hits and `EcalPidTransformer` recomputes the shower
  variables on the filtered set. The branch is optional: pre-mask `ecal` files
  keep all hits.

  **Hit-level MIP cut** (`HitMipCut` property): hits with `energy < threshold`
  are dropped the same index-aligned way, so every shower variable is recomputed
  on the cleaned hit set. Disabled by a negative value (the default when run
  directly); `run_pid_batch.py` sets it to `0.5` in physics mode and disables it
  in `--validation` mode.
- **`EcalPidTransformer`** (`src/components/EcalPidTransformer.cpp`) — one input
  `CalorimeterHitCollection` → one `Cluster`. The physics lives in
  `include/k4SiWEcalReco/EcalShowerVars.h` (a C++ port of
  `siwecal_validation/metrics.py`, the parity oracle). All derived variables go
  into `Cluster::shapeParameters`; their names are published as the frame-level
  metadata parameter `ECalPid_shapeParameterNames`.

### shapeParameters layout (canonical, see `EcalShowerVars.h`)
```
[ scalarNames() ]                                              21 base scalars
[ hits_per_layer[15] | energy_per_layer[15] | weighte_per_layer[15] ]
[ mip05_<scalarNames()> ] [ mip1_<scalarNames()> ]            MIP-cut variants
```
The `mip05_/mip1_` variant blocks are computed **only in `--validation` mode**
(they feed the `event_viewer`'s interactive threshold slider). In the default
physics mode the hits are already cleaned by the `0.5` MIP hit cut, so the blocks
are omitted and the layout is just the 21 base scalars + the three per-layer
profiles. Readers (`PidFileReader`) auto-detect which layout a file uses.

## Build & run (under key4hep)

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
cmake -S k4SiWEcalReco -B k4SiWEcalReco/build -DCMAKE_BUILD_TYPE=Release
cmake --build k4SiWEcalReco/build -j4

export LD_LIBRARY_PATH=$PWD/k4SiWEcalReco/build:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/k4SiWEcalReco/build/genConfDir:$PWD:$PYTHONPATH
```

### Batch driver (recommended): `run_pid_batch.py`

Resolves inputs exactly like `siwecal_validation` (`--run/--file/--all/--point/--cfg`),
runs the Gaudi stage per `ecal_<run>.root`, then applies the cuts and writes the
final output(s):

```bash
# physics mode (default): hit-MIP cut >= 0.5, energy>0, no mip05_/mip1_ blocks
python k4SiWEcalReco/run_pid_batch.py --run TB2026CERN_run_000007 --outdir /tmp/pid

# both formats + an event cut (cuts use the same flags as siwecal_validation)
python k4SiWEcalReco/run_pid_batch.py --all --format both --nhit-min 20 --outdir /tmp/pid

# validation/visualizer mode: no hit cut, compute the mip05_/mip1_ slider blocks
python k4SiWEcalReco/run_pid_batch.py --run TB2026CERN_run_000007 --validation
```

Key flags:

| Flag | Meaning |
|---|---|
| `--format {edm4hep,valtree,both}` | output format(s); default `edm4hep` |
| `--validation` | visualizer mode: compute `mip05_/mip1_` blocks and disable the hit cut |
| `--hit-mip-cut FLOAT` | override the per-hit MIP cut (default `0.5`; off with `--validation`) |
| `--<var>-min` / `--<var>-max`, `--is-shower`/`--no-shower` | event selection cuts (one pair per `CutSet` variable, identical to `siwecal_validation`) |

Cut policy: all cuts are off by default **except total per-event energy > 0**,
which is always enforced. In `--all`/`--point` a per-energy `cuts:` block in the
YAML overrides the CLI cuts (the two-level `general.merge(per_energy)` scheme).
Outputs are named `ecal_<label>.edm4hep.root` and `ecal_<label>.valtree.root`,
next to each input or under `--outdir` / `settings.yml` `pid_dir`.

### Single file via `k4run` (low level)

```bash
# ECAL_HIT_MIP_CUT (<0 disables) and ECAL_MIP_THRESHOLDS ("" = no variant blocks)
# control the mode; ECAL_PID_OUT receives the *unfiltered* EDM4hep file.
ECAL_FILE=/path/ecal_<run>.root ECAL_PID_OUT=out.root \
ECAL_HIT_MIP_CUT=0.5 ECAL_MIP_THRESHOLDS= \
    k4run k4SiWEcalReco/options/run_pid.py
```
(`k4run` alone does not apply the event-level cuts — use `run_pid_batch.py` for
the cut-passing, correctly-named outputs.)

## Verification (parity vs `metrics.py`)

Validated on `TB2026CERN_run_000007` (6342 events) against a freshly recomputed
`metrics.py` oracle: **all 21 scalar variables agree to float precision**,
NaN-consistent. Worst cases: `energy`/`weighte` ≤5e-4 absolute (~1e-7 relative —
inherent to summing float32 hit energies), every other variable ≤8e-6.

Matching the oracle to this level required carrying the weight/barycenter/cumulant
math in **double** (not float): `W/X0`, the Molière containment fraction and the
energy-ordered cumulant are all 1-ULP-sensitive at the 90%-containment boundary,
where a float32 truncation can pick a neighbouring hit (a ~0.3 mm jump on rare
events). See `EcalShowerVars.h`.

> Note: a *stale* `*.valcache.root` checked in next to the data was produced with
> an older pad-map sign convention (negated x/y). The EDM4hep output matches the
> **current** `ecal` tree; regenerate caches before comparing.
