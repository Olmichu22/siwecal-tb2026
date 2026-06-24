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
        │                                   │
        ▼                                   ▼
 ┌─────────────────────────┐        ┌─────────────────────────┐
 │   siwecal_validation    │        │      event_viewer       │
 │  plots + particle-ID    │        │   interactive Dash app  │
 │  metrics, *.valcache    │        │   (3-D / 2-D / dists)    │
 └─────────────────────────┘        └─────────────────────────┘
```

## Packages

| Package | What it does | Docs |
|---|---|---|
| [`siwecal_eventbuilder`](siwecal_eventbuilder/README.md) | Turns the decoded raw tree into reconstructed `ecal` events (BCID clustering, calibration, pad/geometry mapping). | builder README |
| [`siwecal_validation`](siwecal_validation/README.md) | Validation plots + particle-discrimination metrics over the `ecal` tree; caches derived variables in `*.valcache.root`. | validation README |
| [`event_viewer`](event_viewer/README.md) | Plotly Dash app to browse events one by one and explore file-level distributions with dynamic cuts and clustering. | viewer README |
| [`event_display`](event_display/README.md) | Standalone ROOT TEve 3-D single-event display. Runs directly under key4hep (no virtualenv needed) — just `source setup.sh` and launch. | display README |
| [`k4SiWEcalReco`](k4SiWEcalReco/README.md) | Gaudi/k4FWCore stage (C++) that turns the `ecal` tree into **EDM4hep** (`CalorimeterHit` + a per-event `Cluster` carrying the particle-ID shower variables in `shapeParameters`). A 1→1 port of `siwecal_validation.metrics`, parity-validated. | k4 reco README |
| `siwecal_common` | Shared infrastructure: `paths`, the single source of truth for every filesystem location (driven by `settings.yml`). | — |

### How they depend on each other

```
event_viewer  ──imports──▶  siwecal_eventbuilder.{geometry, pad_map}
              └─imports──▶  siwecal_validation.metrics
siwecal_eventbuilder ─┐
siwecal_validation  ──┴──▶  siwecal_common.paths
```

`siwecal_eventbuilder` and `siwecal_validation` are independent of each other;
the viewer reuses geometry from the builder and metrics from the validation. All
three resolve paths through `siwecal_common.paths`, so **no module hard-codes an
absolute `/eos` path** — change `settings.yml` and everything follows.

## The software stack

- **[key4hep](https://key4hep.github.io/)** from CVMFS provides the scientific
  Python stack used everywhere: `numpy`, `scipy`, `pandas`, `scikit-learn`,
  `matplotlib`, `pyyaml`, `uproot`, `awkward` and **ROOT** (PyROOT). The event
  builder and validation write/read ROOT trees with PyROOT; the validation and
  viewer read with `uproot` (no PyROOT needed there).
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
geometry_dir: ./geometry                    # pad maps, slab z, tungsten map
configs_dir:  ./configs/data                # data_reference*.yml run lists
output_dir:   ./validation_output           # validation plots / results
cache_dir:    null                          # *.valcache.root location (null = next to input)
```

Notes:
- `data_dir` may be a **list of search roots**: inputs are looked up in each, in
  order, and the first match wins. Absolute paths written inside the
  `data_reference*.yml` files always take precedence.
- `cache_dir` redirects the validation's `*.valcache.root` caches to one
  directory — set it when your data directory is read-only.
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
├── geometry/                             pad maps, slab_z_positions.yml, Tungsten_thickness.yml
├── configs/data/                         data_reference*.yml run lists
└── calibration/                          vendored dummy pedestal/MIP files
```

Heavy per-run data (`ecal_<run>.root`, `*.valcache.root`, raw inputs) live
**outside** the repo, under `data_dir`.

## End-to-end example

```bash
source setup.sh

# 1. Build events for one run (or a whole energy point with --energy / --all)
python -m siwecal_eventbuilder --run TB2026CERN_run_000007

# 2. Validate it: plots + metrics, caching derived variables
python -m siwecal_validation --run TB2026CERN_run_000007

# 3. Inspect events interactively (open via SSH tunnel, see viewer README)
python -m event_viewer

# 3b. Or use the standalone TEve 3-D display (native key4hep, requires X11/ssh -Y)
cd event_display/
./launch.sh ../data/TB2026CERN_run_000007/ecal_TB2026CERN_run_000007.valcache.root
```

See each package's README for the full option list and design notes.
