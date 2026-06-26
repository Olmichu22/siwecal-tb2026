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
`siwecal_common.paths`. A documented template lives at the repo root —
`cp config.example.yml config.yml` and uncomment only what you need.

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
- **Pad maps**: `mappings/fev10_rotate_chip_channel_x_y_mapping.txt` (default
  wafer) and `mappings/fev11_cob_good_rotate_chip_channel_x_y_mapping.txt`
  (FEV11 chip-on-board, slab 12 override). These are the *rotated* maps that put
  every layer in its true detector orientation.
- **Slab z**: `mappings/slab_z_positions.yml` (live source of truth for `hit_z`).
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

Channels with no MIP calibration (`mpv = 0`) cannot be turned into an energy.
They are kept in the `ecal` tree as a raw record (`hit_energy = 0`) but tagged
with the `hit_ismasked` branch (`1` = masked). The reconstruction drops them:
`k4SiWEcalReco`'s `EcalToEDM4hep` omits masked channels from its output, so they
never reach the shower variables or the `siwecal_validation` plots.

## Calibration file format and NaN handling

### Muon calibration tables (`MuonCalib_it2*`)

`--th N` loads the cumulative MuonCalib tables for threshold `th<N>` from
**`calibration/MuonCalib_it2_corrected/`** — this `_corrected` set is the one the
builder uses (see `MUON_CALIB_DIR` in `cli.py`). It is the corrected variant of
the original `MuonCalib_it2/`: in `_corrected` the pedestal sentinels have been
resolved — finite means are kept, missing SCAs are **imputed** from the first
valid SCA in the row, and channels whose 15 SCAs are **all** sentinel are masked
out entirely. The original `MuonCalib_it2/` is kept alongside for reference/
provenance; the layout (`mips/th<N>/`, `pedestals/th<N>/`) is identical, only the
pedestal values differ. The same sentinel logic is also applied at read time by
`_read_pedestal_file` (below), so a raw set still works — `_corrected` just bakes
it into the files.

### Pedestal file format

```
layer  chip  channel   mean0 err0 wid0   mean1 err1 wid1  ...  mean14 err14 wid14
```

45 columns per row: `(slab_id, chip_id, channel)` followed by 15 × 3 values
(mean, error, width) for each SCA (Switched Capacitor Array cell, 0–14).

The calibration tool writes **`-nan`** as a sentinel for SCAs that had
insufficient statistics (too few triggers on that SCA index). This is normal:
muon runs may not cycle through all 15 SCAs evenly.

### How sentinel entries are treated (current approach — under discussion)

The calibration tool marks uncalibrated SCAs with two interchangeable sentinels:
`-nan` (insufficient statistics) and `0` (no data). Both are treated identically.

`_read_pedestal_file` applies the following logic per row:

| Condition | Action |
|-----------|--------|
| Some SCAs are sentinel (`0` or `-nan`), at least one is a valid non-zero finite value | Sentinel SCAs are replaced by the **first valid mean** in the same row. |
| **All** 15 SCA means are sentinels | The whole channel is added to the **masked set** — its hits get `hit_energy = 0` and `hit_ismasked = 1`, and are excluded from all downstream metrics. |

In `MuonCalib_it2_corrected/run_000142`: 395 channels have all-zero SCAs (→ masked), 11
have a mix of NaN and valid SCAs (→ substitution), none have zeros mixed with
valid values.

**Rationale**: the pedestal of a channel varies only slightly between SCAs
(electronic baseline drift is small). A channel that has at least one calibrated
SCA can safely use that value for the uncalibrated SCAs. A channel with no
calibrated SCA at all has no reference and is treated the same way as a MIP-masked
channel.

### Open questions for team discussion

The following choices were made pragmatically; they should be reviewed with the
calibration team before being treated as the reference strategy:

1. **First-valid substitution vs. averaging**: when multiple SCAs have finite
   means, using the *first* valid one rather than the *average* or a
   *nearest-neighbour* interpolation was chosen for simplicity. If the pedestal
   drift between SCAs is measurable and systematic, an average or a weighted
   interpolation may be more accurate.

2. **Masking all-NaN channels vs. falling back to `pedestal_fallback`**: the
   previous behaviour was to use `pedestal_fallback = 250 ADC` (a global
   constant) when a key was missing from the map, which applied silently to *any*
   uncalibrated channel — including all-NaN ones. The new behaviour masks them
   explicitly. The physical question is: should an uncalibrated channel produce
   hits at all (with a rough global fallback), or should it be suppressed
   entirely? If suppressing is too aggressive and hides real signal, one could
   revert all-NaN channels to the fallback instead.

3. **MIP file**: currently channels with `mpv = 0` are masked (separate
   mechanism in `_read_mip_file`). `-nan` MPV values are *not* present in the
   MIP file (the tool only writes `0` there) so no equivalent fix is needed
   there yet.

4. **Coverage in `MuonCalib_it2_corrected`** (`run_000142` file): calibration coverage is
   nearly complete — most chips on every layer have finite pedestal means.
   Known gaps: chips 0–3 on layer 0, and chip 9 on layers 6 and 13. Those
   channels use `pedestal_fallback` (250 ADC) and the global median MPV. Whether
   250 ADC is a good estimate for those specific chips should be verified against
   raw ADC distributions.
