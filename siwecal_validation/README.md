# siwecal_validation

Validation plots and **particle-identification metrics** for the SiW-ECAL.
For one file or a whole energy scan it reads the per-event variables, applies
selection cuts, fits the energy peak, and writes PNG plots plus a results table.
It is the object-oriented successor of `plot_val_plots.py`.

The per-event metrics are **not computed here** — they are read from the
[`k4SiWEcalReco`](../k4SiWEcalReco/README.md) output next to each input
(`ecal_<run>.edm4hep.root`, or its `ecal_<run>.valtree.root` tree if the EDM4hep
file is absent). This module is plots-only; generate the metrics first with
`k4SiWEcalReco/run_pid_batch.py`.

## How it works

For each input file the pipeline (one class per concern) does:

1. **Read** (`event_data.py`) — load every per-event variable from the
   `k4SiWEcalReco` output: basic ones (`nhit`, `zbary`, `energy`,
   `hits_per_layer`, `mip_likeness`) and the particle-discrimination metrics
   (shower start/length, tungsten-weighted energy, Molière radius, transverse
   RMS, ...). Metrics come straight from the EDM4hep Cluster (or the valtree's
   derived branches); nothing is recomputed. Events with no hits or non-positive
   energy were already dropped upstream.
2. **Select** (`selection.py`) — apply a `CutSet` (general cuts from the CLI, or
   per-energy cuts from the YAML, which override the general ones).
3. **Plot + fit** (`plots.py`, `fits.py`) — one `Plotter` subclass per plot type
   (energy histogram with a Gaussian fit, MIP-likeness, nHit-vs-⟨Z⟩, ...),
   produced individually and/or as a combined grid.
4. **Results** (`results.py`, `output.py`) — a structured output tree under
   `output_dir`, plus `results.csv` / `results.txt` (signal rate, fitted μ/σ,
   summary metrics). In `--all`/`--point` mode it also makes cross-energy summary
   plots (energy calibration, resolution vs E).

`PlotConfig` (`config.py`) holds the tunables (layer count, fit windows,
tungsten map, shower thresholds, ...). Tungsten thickness and slab z default to
values matching `mappings/Tungsten_thickness.yml` / `slab_z_positions.yml`.

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
| `--file FILE` | a specific event-builder `ecal_<run>.root` (its metrics file is found next to it) |
| `--run RUN` | a single run (resolved against the `settings.yml` data roots) |
| `--all` / `--point N` | every entry of the YAML `event_data` map / only beam point P*N* |
| `--cfg FILE` | `data_reference` YAML (default: `configs/data/data_reference.yml`) |
| `--out DIR` | output base (default: `settings.yml` `output_dir`) |
| `--grid-only` / `--no-grid` | only the combined grid PNG / only the individual PNGs |
| `--<var>-min` / `--<var>-max` | general selection cuts (one pair per `CutSet` variable) |
| `--is-shower` / `--no-shower` | keep only shower / non-shower events |

Run `python -m siwecal_validation --help` for the complete list.

## Metrics source

The per-event metrics are read from the `k4SiWEcalReco` output located next to
each input (or in `settings.yml` `pid_dir`): the EDM4hep PID file
`ecal_<run>.edm4hep.root` first, else its `ecal_<run>.valtree.root` tree (the
same derived variables in a plain TTree). If neither exists the run errors,
asking you to generate one with `k4SiWEcalReco/run_pid_batch.py`. This module
never writes those files.

## Stack

Python with `uproot` (reading), `numpy`, `scipy` (Gaussian fits), `matplotlib`
(PNGs, styled via `MLPConfig/newams.mplstyle`), `pyyaml`. Provided by key4hep —
see the top-level [README](../README.md). Run the unit tests with
`pytest siwecal_validation/tests`.
