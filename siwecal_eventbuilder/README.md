# siwecal_eventbuilder

Reconstructs physics **events** for the SiW-ECAL prototype from the decoded raw
data. It reads the converted `siwecaldecoded` ROOT tree, groups channel hits into
events by their bunch-crossing time, calibrates each hit into an energy, attaches
a transverse `(x, y)` and longitudinal `z` position, and writes a flat `ecal`
tree — one entry per reconstructed event — ready for analysis.

This is the object-oriented successor of the original `build_events.py`; its
defaults reproduce the reference event builder (3191 events on run 7 / 74 GeV).

## How it works

The raw tree stores, per acquisition, every triggered SCA cell. Event building
turns that into events in a few well-separated stages:

1. **BCID clustering** (`bcid_clustering.py`) — each SCA carries a *bunch-crossing
   ID*. Nearby BCIDs are merged into time windows; noisy starts, known-artefact
   BCIDs and windows spanning too few slabs are dropped.
2. **Hit collection** (`hit_collector.py`) — for every channel inside a window,
   keep the ones that fired and turn their ADC into a `Hit`.
3. **Calibration** (`calibration.py`) — subtract the per-channel **pedestal** and
   divide by the **MIP** scale to get a physical energy. Pedestal/MIP can be read
   from files or computed on the fly from dedicated runs.
4. **Geometry & pad map** (`geometry.py`, `pad_map.py`) — map each hit's
   `(slab, chip, channel)` to a position: `(x, y)` from the pad map text files,
   `z` from `slab_z_positions.yml`.
5. **Assembly & output** (`event_builder.py`, `pipeline.py`, `root_io.py`) — a
   parallel driver splits the run across worker processes, each builds events and
   writes a temporary ROOT file; the parts are merged with `hadd`. `root_io.py`
   is the only module that talks to PyROOT.

Configuration precedence (lowest → highest):

```
dataclass defaults  →  config.yml  →  data_reference.yml  →  command-line flags
```

`BuilderConfig` (`config.py`) holds every threshold and cut as a default; an
optional `config.yml` overrides only the keys you list (sections `builder`,
`geometry`, `paths`, `calibration`, `mapping`, parsed by `settings.py`); the CLI
wins last. Filesystem locations default to the shared `settings.yml` via
`siwecal_common.paths`.

## Usage

```bash
source setup.sh                      # from the repo root

# One run -> <data_dir>/<run>/ecal_<run>.root
python -m siwecal_eventbuilder --run TB2026CERN_run_000007

# A whole energy point (group of runs from data_reference.yml) into one file
python -m siwecal_eventbuilder --energy P1_20GeV --outdir /tmp/ecal_20gev

# Every entry of a data-reference file
python -m siwecal_eventbuilder --all --data-reference configs/data/data_reference_base.yml
```

### Main options

| Flag | Meaning |
|---|---|
| `--run RUN` / `--energy LABEL` / `--all` | what to build: a single run, an energy point, or everything in the data-reference file |
| `--data-reference PATH` | run-list YAML (energy → runs). Default: `configs/data/data_reference_base.yml` |
| `--config PATH` | optional `config.yml` of tunables/paths overrides |
| `--output FILE` / `--outdir DIR` | output file (single job) or output directory |
| `--workers N` | parallel worker processes (default 5) |
| `--max-entries N` | cap input entries (quick tests) |
| `--calib-dir DIR` / `--ped-file` / `--mip-file` | calibration source files (default: `calibration/` from `settings.yml`) |
| `--compute-calib` / `--ped-run` / `--mip-run` | compute pedestal/MIP from dedicated runs instead of files |
| `--no-calib` | skip calibration (raw ADC) |
| `--no-mapping` | skip the pad `(x, y)` mapping |
| `--exclude LABEL...` / `--ref-out PATH` | (with `--all`) skip labels / write an event-data reference |

Run `python -m siwecal_eventbuilder --help` for the complete list.

## Inputs the package needs (resolved via `settings.yml`)

- **Calibration**: `calibration/dummy_pedestal_15_highgain.txt`,
  `dummy_mip_map_15_highgain.txt` (vendored; override with `--calib-dir`).
- **Pad maps**: `geometry/fev10_rotate_chip_channel_x_y_mapping.txt` (default
  wafer) and `geometry/fev11_cob_good_rotate_chip_channel_x_y_mapping.txt`
  (FEV11 chip-on-board, slab 12 override). These are the *rotated* maps that put
  every layer in its true detector orientation.
- **Slab z**: `geometry/slab_z_positions.yml` (live source of truth for `hit_z`).
- **Run lists**: `configs/data/data_reference*.yml`.

## Stack

Python + **PyROOT** (read/write the ROOT trees), `numpy`, `pyyaml`,
`multiprocessing` for the parallel driver, `hadd` (ROOT) to merge parts. Provided
by key4hep — see the top-level [README](../README.md).

## Output tree (`ecal`)

One entry per reconstructed event with per-event scalars (event id, BCID, number
of hits, total energy, ...) and per-hit arrays (`hit_slab`, `hit_chip`,
`hit_chan`, `hit_x`, `hit_y`, `hit_z`, `hit_energy`, ...). `hit_energy` can be
negative when an ADC sits below pedestal — expected for the current dummy
calibration, not a physics signal.
