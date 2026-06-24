# siwecal_validation

Validation plots and **particle-identification metrics** for the reconstructed
`ecal` tree produced by [`siwecal_eventbuilder`](../siwecal_eventbuilder/README.md).
For one file or a whole energy scan it computes per-event variables, applies
selection cuts, fits the energy peak, and writes PNG plots plus a results table.
It is the object-oriented successor of `plot_val_plots.py`.

## How it works

For each input file the pipeline (one class per concern) does:

1. **Read + derive** (`event_data.py`, `metrics.py`) — load the `ecal` tree with
   `uproot` and compute every per-event variable: basic ones (`nhit`, `zbary`,
   `energy`, `hits_per_layer`, `mip_likeness`) and the particle-discrimination
   metrics (shower start/length, tungsten-weighted energy, Molière radius,
   transverse RMS, ...). Events with no hits or non-positive energy are skipped.
2. **Cache** (`vars_cache.py`) — the full original tree plus all derived
   variables is written to a `<stem>.valcache.root` so re-runs are instant. The
   cache stores a fingerprint of the config; a config change invalidates it
   automatically.
3. **Select** (`selection.py`) — apply a `CutSet` (general cuts from the CLI, or
   per-energy cuts from the YAML, which override the general ones).
4. **Plot + fit** (`plots.py`, `fits.py`) — one `Plotter` subclass per plot type
   (energy histogram with a Gaussian fit, MIP-likeness, nHit-vs-⟨Z⟩, ...),
   produced individually and/or as a combined grid.
5. **Results** (`results.py`, `output.py`) — a structured output tree under
   `output_dir`, plus `results.csv` / `results.txt` (signal rate, fitted μ/σ,
   summary metrics). In `--all`/`--point` mode it also makes cross-energy summary
   plots (energy calibration, resolution vs E).

`PlotConfig` (`config.py`) holds the tunables (layer count, fit windows,
tungsten map, shower thresholds, ...). Tungsten thickness and slab z default to
values matching `geometry/Tungsten_thickness.yml` / `slab_z_positions.yml`.

## Usage

```bash
source setup.sh                      # from the repo root

# Single run (resolved to <data_dir>/<run>/ecal_<run>.root, searching all roots)
python -m siwecal_validation --run TB2026CERN_run_000007

# A specific file, with cuts
python -m siwecal_validation --file data/ecal_P2_7_5GeV_runs_000023_000024.root \
    --nhit-min 20 --energy-min 30

# Whole scan from the data_reference YAML, or just one beam point
python -m siwecal_validation --all
python -m siwecal_validation --point 2
```

### Main options

| Flag | Meaning |
|---|---|
| `--file FILE` | a specific `ecal`/`valcache` ROOT file |
| `--run RUN` | a single run (resolved against the `settings.yml` data roots) |
| `--all` / `--point N` | every entry of the YAML `event_data` map / only beam point P*N* |
| `--cfg FILE` | `data_reference` YAML (default: `configs/data/data_reference.yml`) |
| `--out DIR` | output base (default: `settings.yml` `output_dir`) |
| `--grid-only` / `--no-grid` | only the combined grid PNG / only the individual PNGs |
| `--no-plots` | produce no plots at all, only the ROOT file(s) (valcache / `--save-tree`) — for the `event_viewer` |
| `--<var>-min` / `--<var>-max` | general selection cuts (one pair per `CutSet` variable) |
| `--is-shower` / `--no-shower` | keep only shower / non-shower events |
| `--create-tree` | force (re)generate the metrics cache (after a config change) |
| `--save-tree` | also write a cut-applied augmented tree under the output dir |
| `--cache-dir DIR` | put the `*.valcache.root` caches here instead of next to the input — **use this when the data directory is read-only** |

Run `python -m siwecal_validation --help` for the complete list.

## Caching: `*.valcache.root`

The first run over a file computes the derived variables and writes a
`<stem>.valcache.root`; later runs read it back (the fast path) unless the config
changed or `--create-tree` is given. By default the cache sits next to the input;
set `cache_dir` in `settings.yml` (or `--cache-dir`) to collect all caches in one
writable directory. The same caches are what `event_viewer` reads to show metrics.

## Stack

Python with `uproot` (reading) + **PyROOT** (writing the valcache/augmented
trees), `numpy`, `scipy` (Gaussian fits), `matplotlib` (PNGs, styled via
`MLPConfig/newams.mplstyle`), `pyyaml`. Provided by key4hep — see the top-level
[README](../README.md). Run the unit tests with `pytest siwecal_validation/tests`.
