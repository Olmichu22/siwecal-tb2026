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
 │      k4SiWEcalReco      │   ─────────────────────►  ecal_<run>.edm4hep.root
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
`k4SiWEcalReco` (in C++, fast). `siwecal_validation` and `event_viewer` only read
that output — neither recomputes metrics nor writes any tree.

## Packages

| Package | What it does | Docs |
|---|---|---|
| [`siwecal_eventbuilder`](siwecal_eventbuilder/README.md) | Turns the decoded raw tree into reconstructed `ecal` events (BCID clustering, calibration, pad/geometry mapping). | builder README |
| [`siwecal_validation`](siwecal_validation/README.md) | Validation plots only: reads the per-event metrics from the `k4SiWEcalReco` output (EDM4hep / valtree), applies cuts, fits the energy peak, writes PNGs + results. Generates no trees. | validation README |
| [`event_viewer`](event_viewer/README.md) | Plotly Dash app to browse events one by one and explore file-level distributions with dynamic cuts and clustering. | viewer README |
| [`event_display`](event_display/README.md) | Standalone ROOT TEve 3-D single-event display. Runs directly under key4hep (no virtualenv needed) — just `source setup.sh` and launch. | display README |
| [`k4SiWEcalReco`](k4SiWEcalReco/README.md) | Gaudi/k4FWCore stage (C++) that computes the particle-ID shower variables (a parity-validated port of `siwecal_validation.metrics`), applies the hit-level + event-level cuts, and writes the cut-passing events as **EDM4hep** (`ecal_<run>.edm4hep.root`) or a plain **valtree** TTree (`ecal_<run>.valtree.root`). | k4 reco README |
| `siwecal_common` | Shared infrastructure: `paths`, the single source of truth for every filesystem location (driven by `settings.yml`). | — |

### How they depend on each other

```
event_viewer  ──imports──▶  siwecal_eventbuilder.{geometry, pad_map}
              └─imports──▶  siwecal_validation.metrics
k4SiWEcalReco ──imports──▶  siwecal_validation.{selection, event_data, vars_cache}
siwecal_eventbuilder ─┐
siwecal_validation  ──┼──▶  siwecal_common.paths
k4SiWEcalReco       ──┘
```

`metrics.py` is the parity oracle: `k4SiWEcalReco` is a C++ port of it, and its
batch driver reuses the validation's `CutSet` and tree schema so the cut logic
and the output branches never drift. `siwecal_validation` then reads the metrics
back from the `k4SiWEcalReco` output (it no longer touches the raw `ecal` tree).
All modules resolve paths through `siwecal_common.paths`, so **no module
hard-codes an absolute `/eos` path** — change `settings.yml` and everything follows.

### Reconstruction stage (Gaudi): `k4SiWEcalReco`

`k4SiWEcalReco` is the **Gaudi/k4FWCore** stage that computes the particle-ID
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
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
cmake -S k4SiWEcalReco -B k4SiWEcalReco/build && cmake --build k4SiWEcalReco/build -j4
export LD_LIBRARY_PATH=$PWD/k4SiWEcalReco/build:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/k4SiWEcalReco/build/genConfDir:$PWD:$PYTHONPATH

# batch driver, same input options as siwecal_validation (--run/--file/--all/--point/--cfg)
python k4SiWEcalReco/run_pid_batch.py --run TB2026CERN_run_000007 --outdir /tmp/pid
python k4SiWEcalReco/run_pid_batch.py --all --format both --nhit-min 20 --outdir /tmp/pid
python k4SiWEcalReco/run_pid_batch.py --run TB2026CERN_run_000007 --validation   # viewer mode

# validation / viewer then consume the output directly (found next to the input)
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
> `k4SiWEcalReco` output — the EDM4hep or valtree file), masked hits *will* show
> up — the `ecal` tree keeps them as a raw record.

See [`k4SiWEcalReco/README.md`](k4SiWEcalReco/README.md) for details.

## The software stack

- **[key4hep](https://key4hep.github.io/)** from CVMFS provides the scientific
  Python stack used everywhere: `numpy`, `scipy`, `pandas`, `scikit-learn`,
  `matplotlib`, `pyyaml`, `uproot`, `awkward` and **ROOT** (PyROOT). The event
  builder and `k4SiWEcalReco` write ROOT trees with PyROOT; `siwecal_validation`
  and the viewer read with `uproot` (no PyROOT needed there).
- **dash + plotly** (the only packages not in key4hep) power the `event_viewer`
  web UI. They live in a local `--system-site-packages` virtualenv,
  `.venv-viewer` (not versioned — recreate it per machine, see below).

## Quick start

```bash
# 1. Clone and enter the repo
git clone <this-repo> siwecal-tb2026 && cd siwecal-tb2026

# 2. Tell the suite where your data lives
cp settings.example.yml settings.yml
$EDITOR settings.yml          # set data_dir to your run/event directory

# 3. Load the environment (key4hep + repo on PYTHONPATH + .venv-viewer if present)
source setup.sh
```

### First-time setup of the viewer virtualenv

`.venv-viewer` is **not** committed (147 MB, with absolute CVMFS paths baked in).
Create it once per machine, on top of key4hep:

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
python -m venv --system-site-packages .venv-viewer
source .venv-viewer/bin/activate
pip install -r requirements.txt        # dash + plotly
```

After that, `source setup.sh` activates it automatically.

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
pid_dir:      null                          # k4SiWEcalReco outputs (null = next to input)
cache_dir:    null                          # legacy *.valcache.root (null = next to input)
```

Notes:
- `data_dir` may be a **list of search roots**: inputs are looked up in each, in
  order, and the first match wins. Absolute paths written inside the
  `data_reference*.yml` files always take precedence.
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
├── setup.sh                              environment loader
├── requirements.txt                      dash + plotly (rest from key4hep)
├── siwecal_common/                       shared path resolution
├── siwecal_eventbuilder/                 event building (+ README)
├── siwecal_validation/                   validation + metrics (+ README)
├── event_viewer/                         Dash viewer (+ README)
├── event_display/                        ROOT TEve 3-D event display (native key4hep, no extra deps)
├── mappings/                             pad maps, slab_z_positions.yml, Tungsten_thickness.yml
├── configs/data/                         data_reference*.yml run lists
└── calibration/                          vendored dummy pedestal/MIP files
```

Heavy per-run data (`ecal_<run>.root`, the `k4SiWEcalReco` outputs
`ecal_<run>.edm4hep.root` / `ecal_<run>.valtree.root`, raw inputs) live
**outside** the repo, under `data_dir`.

## End-to-end example

```bash
source setup.sh

# 1. Build events for one run (or a whole energy point with --energy / --all)
python -m siwecal_eventbuilder --run TB2026CERN_run_000013 --th 220

# 2. Reconstruct: shower variables + cuts -> ecal_<run>.edm4hep.root (needs key4hep,
#    see the k4SiWEcalReco README; add --validation for the viewer's slider blocks)
python k4SiWEcalReco/run_pid_batch.py --run TB2026CERN_run_000007

# 3. Validate it: plots + fits with the cuts you want (reads the step-2 output)
python -m siwecal_validation --run TB2026CERN_run_000007 --nhit-min 20

# 4. Inspect events interactively (open via SSH tunnel, see viewer README)
python -m event_viewer

# 4b. Or use the standalone TEve 3-D display (native key4hep, requires X11/ssh -Y)
cd event_display/
./launch.sh ../data/TB2026CERN_run_000007/ecal_TB2026CERN_run_000007.valtree.root
```

See each package's README for the full option list and design notes.
