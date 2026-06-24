# k4SiWEcalReco

Gaudi/k4FWCore stage that turns the per-event reconstructed `ecal` tree (written
by `siwecal_eventbuilder`) into **EDM4hep**, computing the particle-discrimination
shower variables in C++.

```
ecal_<run>.root  ──►  EcalToEDM4hep  ──►  CalorimeterHitCollection (+ EventHeader)
   (1 row/event)        (1→1 source)              │
                                                  ▼
                                          EcalPidTransformer
                                                  │
                                                  ▼
                              ClusterCollection (1 Cluster/event), the shower
                              variables flattened into Cluster::shapeParameters
```

Because the `ecal` tree is already one row per physics event, the whole chain is
a strict **1→1** transform — no BCID fan-out, no custom event-loop machinery. The
source sets `EvtMax = GetEntries()` and serves one event per call.

## Components

- **`EcalToEDM4hep`** (`src/components/EcalToEDM4hep.cpp`) — reads the (non-podio)
  `ecal` tree with `TTreeReader` and emits, per event, a `CalorimeterHitCollection`
  (`energy = hit_energy`, `position = (hit_x, hit_y, hit_z)`, `cellID` encoding
  `slab/chip/channel/sca`) and an `EventHeaderCollection` (run, event, bcid, spill).
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

## Build & run (under key4hep)

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
cmake -S k4SiWEcalReco -B k4SiWEcalReco/build -DCMAKE_BUILD_TYPE=Release
cmake --build k4SiWEcalReco/build -j4

export LD_LIBRARY_PATH=$PWD/k4SiWEcalReco/build:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/k4SiWEcalReco/build/genConfDir:$PYTHONPATH
ECAL_FILE=/path/ecal_<run>.root ECAL_PID_OUT=out.root \
    k4run k4SiWEcalReco/options/run_pid.py
```

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
