# SiW-ECAL TB2026CERN Analysis Suite

End-to-end analysis tools for the **SiW-ECAL** (Silicon-Tungsten Electromagnetic
Calorimeter) prototype tested at the CERN **TB2026CERN** beam test. The suite
takes the decoded raw data, builds physics events, validates them, and lets you
inspect them interactively — three packages that share one geometry, one
calibration and one configuration file.

```
 RAW2ROOT (.root, siwecaldecoded tree)        <-- produced upstream (not in this repo)
        │
        ▼
 ┌─────────────────────────┐   reconstructed events
 │  siwecal_eventbuilder   │   ─────────────────────►  ecal_<run>.root  (ecal tree)
 └─────────────────────────┘
        │
        ▼
 ┌─────────────────────────┐   per-event metrics + cuts (C++), one of two formats:
 │      gaudi_source       │   ─────────────────────►  ecal_<run>.edm4hep.root
 │  shower vars + selection│                           ecal_<run>.valtree.root
 └─────────────────────────┘
        │                                   │
        ▼                                   ▼
 ┌─────────────────────────┐        ┌─────────────────────────┐
 │   siwecal_validation    │        │      event_viewer       │
 │   validation plots only │        │   interactive Dash app  │
 │   (reads the metrics)   │        │   (3-D / 2-D / dists)    │
 └─────────────────────────┘        └─────────────────────────┘
```

The per-event metrics **and the cut/cleaned event collections** are produced by
`gaudi_source` (in C++, fast). `siwecal_validation` and `event_viewer` only read
that output — neither recomputes metrics nor writes any tree.

## Packages

| Package | What it does | Docs |
|---|---|---|
| [`siwecal_eventbuilder`](siwecal_eventbuilder/README.md) | Turns the decoded raw tree into reconstructed `ecal` events (BCID clustering, calibration, pad/geometry mapping). | builder README |
| [`siwecal_validation`](siwecal_validation/README.md) | Validation plots only: reads the per-event metrics from the `gaudi_source` output (EDM4hep / valtree), applies cuts, fits the energy peak, writes PNGs + results. Generates no trees. | validation README |
| [`event_viewer`](event_viewer/README.md) | Plotly Dash app to browse events one by one and explore file-level distributions with dynamic cuts and clustering. | viewer README |
| [`event_display`](event_display/README.md) | Standalone ROOT TEve 3-D single-event display. Runs directly under key4hep (no virtualenv needed) — just `source setup.sh` and launch. | display README |
| [`gaudi_source`](gaudi_source/README.md) | Gaudi/k4FWCore stage (C++) that computes the particle-ID shower variables (a parity-validated port of `siwecal_validation.metrics`), applies the hit-level + event-level cuts, and writes the cut-passing events as **EDM4hep** (`ecal_<run>.edm4hep.root`) or a plain **valtree** TTree (`ecal_<run>.valtree.root`). | gaudi README |
| [`gaudi_jobs`](gaudi_jobs/) | Driver de batch (`run_pid_batch.py`) y ejemplos de jobs concretos (e.g. `run000013/`). | — |
| `siwecal_common` | Shared infrastructure: `paths`, the single source of truth for every filesystem location (driven by `settings.yml`). | — |

### How they depend on each other

```
event_viewer  ──imports──▶  siwecal_eventbuilder.{geometry, pad_map}
              └─imports──▶  siwecal_validation.metrics
gaudi_source  ──imports──▶  siwecal_validation.{selection, event_data, vars_cache}
siwecal_eventbuilder ─┐
siwecal_validation  ──┼──▶  siwecal_common.paths
gaudi_source        ──┘
```

`metrics.py` is the parity oracle: `k4SiWEcalReco` is a C++ port of it, and its
batch driver reuses the validation's `CutSet` and tree schema so the cut logic
and the output branches never drift. `siwecal_validation` then reads the metrics
back from the `k4SiWEcalReco` output (it no longer touches the raw `ecal` tree).
All modules resolve paths through `siwecal_common.paths`, so **no module
hard-codes an absolute `/eos` path**. `settings.yml` is the default data root; a
`data_reference*.yml` may set `main_path` to override it for a whole run list (and
a per-entry absolute `path:` wins over both), so the run lists stay the single
convenient place to point the suite at your data.

### Reconstruction stage (Gaudi): `gaudi_source` + `gaudi_jobs`

`gaudi_source` is the **Gaudi/k4FWCore** stage that computes the particle-ID
shower variables in C++ (a parity-validated port of `siwecal_validation.metrics`)
and produces the **cut-passing event collections** consumed downstream. Per run
it applies an optional hit-level MIP cut and the event selection (the same
`CutSet` as `siwecal_validation`), then writes one of two formats:

- **EDM4hep** PID file `ecal_<run>.edm4hep.root` — one `Cluster` per event
  (variables in `shapeParameters`) + `CalorimeterHit` collections.
- **valtree** `ecal_<run>.valtree.root` — the same per-event variables in a plain
  TTree (the valcache schema).

Two modes: the default **physics** mode bakes in a `0.5` MIP per-hit cut, while
`--validation` mode keeps raw hits and also computes the `mip05_/mip1_` variant
blocks that feed the viewer's interactive threshold slider. All cuts are off by
default except total per-event energy > 0 (always enforced).

```bash
# Build it once (./install.sh does this for you) and load the env:
./install.sh --no-viewer      # or the manual cmake build, see k4SiWEcalReco/README.md
source setup.sh               # puts the build on LD_LIBRARY_PATH/PYTHONPATH for k4run

# batch driver, same input options as siwecal_validation (--run/--file/--all/--point/--cfg)
python gaudi_jobs/run_pid_batch.py --run TB2026CERN_run_000007
python gaudi_jobs/run_pid_batch.py --all --format both --nhit-min 20
python gaudi_jobs/run_pid_batch.py --run TB2026CERN_run_000007 --validation   # viewer mode

# concrete job example (steering file for run000013, physics mode)
k4run gaudi_jobs/run000013/steer_run000013.py

# validation / viewer then consume the output directly
python -m siwecal_validation --run TB2026CERN_run_000007
```

> ⚠️ **Masked channels.** Two classes of channels are flagged `hit_ismasked = 1`
> in the `ecal` tree and excluded from all downstream stages (EDM4hep / valtree,
> metrics, plots):
> - **No MIP calibration** (`mpv = 0` in the MIP file).
> - **No valid pedestal** (all 15 SCA means are `-nan` in the pedestal file — the
>   calibration tool's sentinel for channels with insufficient statistics).
>
> Channels with only *partial* NaN coverage (some SCAs finite, some `-nan`) are
> **not** masked; the first valid SCA mean is substituted for the missing ones.
> See [`siwecal_eventbuilder/README.md`](siwecal_eventbuilder/README.md) for the
> full discussion and open questions.
>
> If you open the viewer directly on a **raw** `ecal_<run>.root` (instead of the
> `gaudi_source` output — the EDM4hep or valtree file), masked hits *will* show
> up — the `ecal` tree keeps them as a raw record.

See [`gaudi_source/README.md`](gaudi_source/README.md) for details.

## Tree schemas

Three on-disk products carry the per-event data, in order of the pipeline:

1. the reconstructed **`ecal`** tree (`ecal_<run>.root`) — raw events, *no metrics*;
2. the **EDM4hep** PID file (`ecal_<run>.edm4hep.root`) — cut-passing events with metrics;
3. the **valtree** (`ecal_<run>.valtree.root`) — the same content as a plain TTree.

In both `k4SiWEcalReco` outputs the masked channels are already dropped, so there
is no `hit_ismasked` branch and the per-event counters/sums are recomputed from
the surviving hits.

### `ecal` tree — reconstructed events (no metrics)

One entry per BCID window. Per-hit branches are variable-length, counted by
`nhit_chan` (also named `nhit` downstream).

| Var          | Description                                            |
| ------------ | ------------------------------------------------------ |
| run          | run number                                             |
| event        | event number. `spill_index * 1000 + event_index`       |
| spill        | spill number                                           |
| bcid         | bcid id                                                |
| nhit_slab    | Number of distinct slabs contributing to this window.  |
| nhit_chip    | Number of distinct ``(slab, chip)`` pairs with hits .  |
| sum_hg       | Sum of pedestal-subtracted high-gain ADC over all hits |
| sum_energy   | Sum of calibrated energy over all hits                 |
| hit_slab     | physical slab slot (ib, 0..14)                         |
| hit_chip     | hardware chip ID (chipid)                              |
| hit_chan     | pixel index (ipix, 0..63)                              |
| hit_sca      | SCA memory cell (isca, 0..14)                          |
| hit_hg       | high-gain ADC minus pedestal                           |
| hit_lg       | low-gain ADC minus pedestal                            |
| hit_energy   | high-gain signal in MIP units                          |
| hit_x        | transverse position in mm                              |
| hit_y        | transverse position in mm                              |
| hit_z        | longitudinal position in mm                            |
| hit_ismasked | masked in ped or MIP file                              |

### EDM4hep PID file (`ecal_<run>.edm4hep.root`)

One podio `events` frame per cut-passing event. Three collections plus five
parallel per-hit `UserDataCollection`s. The legacy column each maps to (read back
by `siwecal_common.edm4hep_pid`) is shown in parentheses.

| Collection / field                        | Description                                                      |
| ----------------------------------------- | --------------------------------------------------------------- |
| `EventHeader.runNumber`                   | run number (→ `run`)                                            |
| `EventHeader.eventNumber`                 | event number (→ `event`)                                        |
| `EventHeader.weight`                      | spill number (→ `spill`)                                        |
| `EventHeader.timeStamp`                   | bcid id (→ `bcid`)                                              |
| `ECalHits` (`CalorimeterHit`) `.type`     | physical slab slot 0..14 (→ `hit_slab`)                         |
| `ECalHits.position` (x, y, z)             | hit position in mm (→ `hit_x` / `hit_y` / `hit_z`)             |
| `ECalHits.energy`                         | calibrated hit signal in MIP units (→ `hit_energy`)            |
| `ECalHitChip` (UserData)                  | hardware chip ID (→ `hit_chip`)                                |
| `ECalHitChan` (UserData)                  | pixel index 0..63 (→ `hit_chan`)                               |
| `ECalHitSca` (UserData)                   | SCA memory cell 0..14 (→ `hit_sca`)                            |
| `ECalHitHG` (UserData)                    | high-gain ADC minus pedestal (→ `hit_hg`)                      |
| `ECalHitLG` (UserData)                    | low-gain ADC minus pedestal (→ `hit_lg`)                       |
| `ECalPid` (`Cluster`) `.shapeParameters`  | the per-event metrics, one cluster/event (see table below)      |

The `shapeParameters` floats are stored in a fixed order (names listed once in
the `metadata` frame parameter `ECalPid_shapeParameterNames`):

| shapeParameter      | Description                                                       |
| ------------------- | ---------------------------------------------------------------- |
| nhit                | number of (unmasked) hit channels in the event                   |
| zbary               | energy-weighted mean slab index (layer units)                    |
| energy              | summed hit energy (MIP units)                                    |
| mip_likeness        | inverse-hit-rate MIP-likeness score (Σ 1/hits_layer / n_layers)  |
| weighte             | tungsten-weighted energy Σ E·W/X0                                |
| bar_x               | energy-weighted transverse barycenter x (mm)                     |
| bar_y               | energy-weighted transverse barycenter y (mm)                     |
| bar_r               | radial barycenter √(bar_x² + bar_y²) (mm)                        |
| moliere             | 90% transverse containment (Molière) radius in mm (showers; 0 else) |
| transverse_rms      | energy-weighted RMS hit radius about the shower axis (mm)        |
| is_shower           | longitudinal shower flag (0/1): rising EM-like edge with a peak  |
| shower_start        | first shower layer (NaN if not a shower)                         |
| shower_max          | peak (maximum) layer of the profile                              |
| shower_end          | last shower layer                                                |
| shower_start_10     | first layer above 10% of the peak                               |
| shower_end_10       | last layer above 10% of the peak                                |
| shower_length       | number of layers above the shower threshold                     |
| first_layer         | lowest layer with a hit                                         |
| last_layer          | highest layer with a hit                                        |
| n_layers_hit        | number of layers with at least one hit                          |
| e_over_nhit         | energy / nhit (hit energy density)                              |
| hits_per_layer_0..14    | per-layer hit count (one float per layer)                   |
| energy_per_layer_0..14  | per-layer Σ E (one float per layer)                         |
| weighte_per_layer_0..14 | per-layer tungsten-weighted energy (one float per layer)    |

In `--validation` mode the same 21 scalars are appended twice more, recomputed
after a per-hit cut, under the prefixes `mip05_` (`hit_energy ≥ 0.5`) and `mip1_`
(`hit_energy ≥ 1.0`) — e.g. `mip05_moliere`, `mip1_is_shower`. These feed the
viewer's interactive threshold slider; physics mode omits them.

### valtree (`ecal_<run>.valtree.root`)

A plain `ecal` TTree (the valcache schema) carrying exactly the EDM4hep content
flattened into branches: the event identifiers, the recomputed counters/sums, the
per-hit arrays, and the per-event metrics.

| Var                     | Description                                                       |
| ----------------------- | ---------------------------------------------------------------- |
| run, event, spill, bcid | event identifiers (as in the `ecal` tree)                        |
| nhit_chan               | number of surviving hit channels (counts the `hit_*` arrays)     |
| nhit_slab               | number of distinct slabs with hits                              |
| nhit_chip               | number of distinct `(slab, chip)` pairs with hits               |
| sum_hg                  | Σ high-gain ADC (pedestal-subtracted) over the surviving hits    |
| sum_energy              | Σ calibrated energy (MIP units) over the surviving hits          |
| hit_slab/chip/chan/sca  | per-hit integer arrays, as in the `ecal` tree (masked removed)   |
| hit_hg/lg               | per-hit high/low-gain ADC minus pedestal                        |
| hit_energy              | per-hit signal in MIP units                                     |
| hit_x/y/z               | per-hit position in mm                                          |
| nhit … e_over_nhit      | the 21 per-event scalar metrics — same as the EDM4hep `shapeParameters` table |
| hits/energy/weighte_per_layer | the three `[n_layers]` per-layer profiles                  |
| mip05_*, mip1_*         | per-hit-cut scalar variants (validation mode only, as above)     |

## The software stack

- **[key4hep](https://key4hep.github.io/)** from CVMFS provides the scientific
  Python stack used everywhere: `numpy`, `scipy`, `pandas`, `scikit-learn`,
  `matplotlib`, `pyyaml`, `uproot`, `awkward` and **ROOT** (PyROOT). The event
  builder and `gaudi_source` write ROOT trees with PyROOT; `siwecal_validation`
  and the viewer read with `uproot` (no PyROOT needed there).
- **dash + plotly** (the only packages not in key4hep) power the `event_viewer`
  web UI. They live in a local `--system-site-packages` virtualenv,
  `.venv-viewer` (not versioned — recreate it per machine, see below).
- **`gaudi_source`** (the Gaudi/k4FWCore plugin) and its jobs driver
  (`gaudi_jobs/run_pid_batch.py`) also require key4hep. Build with CMake once
  and export `LD_LIBRARY_PATH` / `PYTHONPATH` as shown above.

## Quick start

```bash
# 1. Clone and enter the repo
git clone <this-repo> siwecal-tb2026 && cd siwecal-tb2026

# 2. Tell the suite where your data lives
cp settings.example.yml settings.yml
$EDITOR settings.yml          # set data_dir to your run/event directory

# 3. Install the stack once (key4hep + .venv-viewer + k4SiWEcalReco build)
./install.sh

# 4. Load the environment in every new shell
source setup.sh
```

### Two scripts: `install.sh` (once) vs `setup.sh` (every shell)

The stack has three pieces; `install.sh` performs the one-time, heavy setup of
all of them and `setup.sh` is the lightweight per-shell environment loader.

| Piece | What it is | `install.sh` does | `setup.sh` does |
|---|---|---|---|
| **key4hep** | the cvmfs scientific stack (numpy/scipy/uproot/ROOT…) | sources it | sources it |
| **`.venv-viewer`** | the viewer's `dash`/`plotly` (not in key4hep, not committed) | creates it + `pip install` | activates it if present |
| **`k4SiWEcalReco`** | the compiled Gaudi plugin `k4run` loads | builds it (cmake) | puts the build on `LD_LIBRARY_PATH`/`PYTHONPATH` |

```bash
./install.sh                 # all three pieces (default)
./install.sh --no-k4         # skip the C++ build (e.g. viewer-only machine)
./install.sh --no-viewer     # skip the dash/plotly venv
./install.sh --force-venv    # recreate .venv-viewer from scratch
./install.sh -j 8            # build with 8 parallel jobs
KEY4HEP_RELEASE=2026-02-01 ./install.sh   # pin a different key4hep release
```

The default key4hep release lives in **`.key4hep-release`** — the single source
of truth read by both `install.sh` and `setup.sh`. Edit that file to move the
whole suite to a new release; the `KEY4HEP_RELEASE` env var overrides it per run.

`.venv-viewer` is **not** committed (it bakes in absolute CVMFS paths), and the
`k4SiWEcalReco/build/` tree is machine-specific — both are recreated by
`install.sh` per machine. After installing, `source setup.sh` wires everything
up: it sources key4hep, puts the repo + the `k4SiWEcalReco` build on
`PYTHONPATH`, and activates `.venv-viewer` automatically.

## Configuration: `settings.yml`

A single file at the repo root, read by all packages through
`siwecal_common.paths`. Copy `settings.example.yml` and edit:

```yaml
data_dir:                                   # heavy run/event ROOT files (NOT in repo)
  - /eos/home-o/oarquero/TB2026CERN         #   single path or a LIST of search roots
calib_dir:    ./calibration                 # pedestal / MIP files (vendored)
geometry_dir: ./mappings                    # pad maps, slab z, tungsten map
configs_dir:  ./configs/data                # data_reference*.yml run lists
output_dir:   ./validation_output           # validation plots / results
pid_dir:      null                          # gaudi_source outputs (null = next to input)
cache_dir:    null                          # legacy *.valcache.root (null = next to input)
```

Notes:
- `data_dir` may be a **list of search roots**: inputs are looked up in each, in
  order, and the first match wins. It is the default data root.
- A `data_reference*.yml` may set `main_path` to point a whole run list at a
  specific root; when present it **takes precedence** over `settings.yml`'s
  `data_dir`. A per-entry absolute `path:` wins over both, and files not found
  under `main_path` fall back to the `data_dir` roots.
- `pid_dir` is where `k4SiWEcalReco` writes (and `siwecal_validation` /
  `event_viewer` look for) `ecal_<run>.edm4hep.root` / `ecal_<run>.valtree.root`
  — set it when your data directory is read-only.
- `cache_dir` is the legacy `*.valcache.root` location (still read by the viewer
  if present); `siwecal_validation` no longer writes caches.
- Point to a different settings file with `export SIWECAL_SETTINGS=/path/to.yml`.

## Repository layout

```
siwecal-tb2026/
├── settings.yml / settings.example.yml   shared configuration
├── install.sh                           one-time installer (key4hep venv + k4 build)
├── setup.sh                              per-shell environment loader
├── requirements.txt                      dash + plotly (rest from key4hep)
├── siwecal_common/                       shared path resolution
├── siwecal_eventbuilder/                 event building (+ README)
├── siwecal_validation/                   validation + metrics (+ README)
├── event_viewer/                         Dash viewer (+ README)
├── event_display/                        ROOT TEve 3-D event display (native key4hep, no extra deps)
├── gaudi_source/                         Gaudi/k4FWCore C++ source + CMake + steering file genérico
├── gaudi_jobs/                           driver de batch (run_pid_batch.py) + ejemplos de jobs
│   └── run000013/                        ejemplo concreto para TB2026CERN_run_000013
├── mappings/                             pad maps, slab_z_positions.yml, Tungsten_thickness.yml
├── configs/data/                         data_reference*.yml run lists
└── calibration/                          vendored dummy pedestal/MIP files
```

Heavy per-run data (`ecal_<run>.root`, the `gaudi_source` outputs
`ecal_<run>.edm4hep.root` / `ecal_<run>.valtree.root`, raw inputs) live
**outside** the repo, under `data_dir`.

## End-to-end example

```bash
source setup.sh

# 1. Build events for one run (or a whole energy point with --energy / --all)
python -m siwecal_eventbuilder --run TB2026CERN_run_000013 --th 220

# 2. Reconstruct: shower variables + cuts -> ecal_<run>.edm4hep.root (needs key4hep,
#    see the gaudi_source README; add --validation for the viewer's slider blocks)
python gaudi_jobs/run_pid_batch.py --run TB2026CERN_run_000007

# 3. Validate it: plots + fits with the cuts you want (reads the step-2 output)
python -m siwecal_validation --run TB2026CERN_run_000007 --nhit-min 20

# 4. Inspect events interactively (open via SSH tunnel, see viewer README)
python -m event_viewer

# 4b. Or use the standalone TEve 3-D display (native key4hep, requires X11/ssh -Y)
cd event_display/
./launch.sh ../data/TB2026CERN_run_000007/ecal_TB2026CERN_run_000007.valtree.root
```

See each package's README for the full option list and design notes.
