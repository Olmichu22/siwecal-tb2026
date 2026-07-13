# gaudi_source

Gaudi/k4FWCore stage covering the **full reconstruction chain**, from raw SLB
DAQ binary files to high-level physics-discrimination variables, entirely in
C++/ROOT:

```
run/*_raw.bin_NNNN  (READ-ONLY DAQ data)
        │
        ▼  EcalRawDecoder            (raw2root: binary -> siwecaldecoded)
siwecaldecoded.root  { TTree "siwecaldecoded" }, 1 row/acquisition cycle
        │
        ├──▼  PedestalMipCalibrator, Mode=Pedestal  ──►  Pedestal_<label>_<gain>gain.txt
        └──▼  PedestalMipCalibrator, Mode=Mip        ──►  MIP_pedestalsubmode1_<label>_<gain>gain.txt
               (both modes read the SAME data-taking run(s) -- see below)
        │
        ▼  EcalEventBuilder          (BCID clustering + calibration + geometry)
ecal_<run>.root  { TTree "ecal" }, 1 row/physics event
        │
        ▼  EcalToEDM4hep             (hit-MIP cut)
CalorimeterHitCollection (+ EventHeader)
        │
        ▼  EcalPidTransformer        (shower variables)
ClusterCollection (1 Cluster/event)
        │
        ▼  run_pid_batch.py: event cuts (CutSet) + format select
        ├──►  ecal_<run>.edm4hep.root   (filtered EDM4hep)
        └──►  ecal_<run>.valtree.root   (plain valcache-schema TTree)
```

Everything above the `EcalToEDM4hep` line used to be the pure-Python
`siwecal_eventbuilder` package (still present, kept as a cross-check
reference for now); it has been ported to Gaudi/C++ so the whole chain runs
as `k4run` processes. Pedestal and MIP calibration are computed from the
**same** data-taking run(s) — the reference tool's single event loop
(`DecodedSLBAnalysis.cc::NSlabsAnalysis`) fills the pedestal histograms from
samples where `hitbit==0` and the MIP histograms from samples where
`hitbit==1` in one pass over one input file, and the deployed
`MuonCalib_it2/{pedestals,mips}/th<N>/` area always pairs a `Pedestal_*` and
a `MIP_*` file for the same run number. Calibration is still a separate
upstream step from the physics run being converted, just not a *different*
data-taking run from each other.

**SAFETY**: the raw DAQ directory (`.../Data/rundata/`) is **READ-ONLY**.
Every driver in `gaudi_jobs/` writes decoded/reconstructed output to a
separate `--converted-dir` and never deletes anything under `.../Data/`.

It is also **the single place where the cut/cleaned event collections are
generated**: `EcalToEDM4hep` applies an optional hit-level MIP cut, then the
event selection (the same `CutSet` as `siwecal_validation`), and writes the
cut-passing events in either of two formats — an **EDM4hep** PID file
(`ecal_<run>.edm4hep.root`) or a plain **valtree** TTree
(`ecal_<run>.valtree.root`). `siwecal_validation` no longer generates any
tree; it only reads these outputs and makes plots.

The `EcalToEDM4hep`/`EcalPidTransformer` chain itself is a strict **1→1**
transform (the `ecal` tree is already one row per physics event — no BCID
fan-out): it writes a *full* EDM4hep file to a temporary path.
`run_pid_batch.py` then applies the event cuts and writes the requested
final, cut-passing output(s), deleting the temporary file.

## Components

### Raw2root, calibration, event building

- **`EcalRawDecoder`** (`src/components/EcalRawDecoder.cpp`,
  `include/k4SiWEcalReco/SlbFrameDecoder.h`) — decodes the binary SLB DAQ
  chunk files of one run (`<run>_raw.bin`, `_raw.bin_0001`, ...) into the
  `siwecaldecoded` ROOT tree: Gray-code ADC/BCID decoding, per-(slab,chip,sca)
  cycle buffering with overflow-jump retry, and the bad-BCID
  (retrigger/empty-event) tagging state machine — a faithful port of the
  external reference tool `SLBraw2ROOT.cc`. Also reads `Run_Settings.txt`
  (via `include/k4SiWEcalReco/RunSettings.h`, property `RunSettingsFile`) and
  writes six extra run-constant branches, repeated on every entry (same
  convention as `threshold_dac`/`run` further down the chain):
  `acqWindowMs`, `delayBetweenCycleMs`, `thresholdDac`, `holdDelay`,
  `fsPeakTime`, `gainSelectionThreshold`. The four per-chip fields are
  aggregated by **mode** (most frequent value across all `## ChipIndex:`
  lines) so a handful of miscalibrated/test chips don't skew the result;
  `acqWindowMs`/`delayBetweenCycleMs` are read once (they appear once per
  file). Verified on real TB2026CERN data (see git history for the
  parity/segfault-fix notes).

- **`PedestalMipCalibrator`** (`src/components/PedestalMipCalibrator.cpp`,
  `include/k4SiWEcalReco/PedestalMipCalib.h`) — computes pedestal (Gaussian
  fit via `TSpectrum` peak search) and MIP (Landau⊗Gauss "langau" fit)
  calibration tables from a `siwecaldecoded` tree, in the exact text format
  `siwecal_eventbuilder/calibration.py::from_files` already reads (`Mode`
  property: `Pedestal` or `Mip`). Ported from
  `SiWECAL-TB-analysis/SLBperformance/{DecodedSLBAnalysis.cc,TBchecks/analysis.h}`
  -- verified to be the method that actually produced this repo's deployed
  `calibration/MuonCalib_it2*` files (header lines match verbatim).

  **Fit-robustness policy** (ported from `calibration/clean_pedestals.py`'s
  already-validated technique directly into the generation step, instead of
  a separate cleanup pass):
  - *Pedestal*: a fitted SCA is "usable" only if the Gaussian fit converged
    **and** its mean is finite, positive, and at or below
    `PedestalMaxMeanAdc` (default `300`, matching `clean_pedestals.py`'s
    `MAX_MEAN`). Every unusable SCA (NaN/Inf, out-of-range, or a failed fit)
    is filled with the **plain average of the channel's other usable SCAs**
    (flagged with the existing `error=-10` sentinel). If a channel has **no**
    usable SCA at all, the whole row is written as 15× `0 -10 0` — an
    all-zero-means row, which `calibration.py::_read_pedestal_file` already
    treats as **masked**.
  - *MIP*: if the Langau fit fails (`ndf<=0`), instead of masking
    immediately the combined histogram's **peak bin** is checked as a rough
    MPV estimate; if positive and at or below `MipFallbackMaxAdc` (default
    `300`; lower it for low gain, e.g. `~30`) it is used, flagged with a new
    `empv=-2` sentinel (fit failed, histogram-peak fallback used) so it's
    distinguishable from a genuine fit (`empv>=0`), a fit that didn't look
    like a real MIP (`empv=-1`), a fit that failed with no usable fallback
    (`empv=-5`), or too-few-entries (`empv=-10`).
  - Cross-masking between the pedestal and MIP files (a channel masked in
    *either* is masked overall) already happens at **read time** in
    `include/k4SiWEcalReco/CalibrationTables.h::fromFiles` (union of both
    masked-channel sets) — the two generated `.txt` files themselves are
    not rewritten against each other.

  **Cross-check plots** (`DiagnosticsFile` property, empty by default): an
  optional ROOT file with 2D `layer*20+chip` vs. `sca*100+channel` (pedestal)
  or `layer*20+chip` vs. `channel` (MIP) maps of the written value **and** a
  per-cell fit-status code (masked / genuine fit / fallback-used / etc. --
  see the sentinel scheme above), plus 1D value distributions
  (`pedestal_mean_dist`, `mip_mpv_dist`, ...). Makes the new robustness
  policy's decisions directly inspectable in a `TBrowser` or a quick PyROOT
  script, without needing to parse the `.txt` tables by hand. `--diagnostics`
  on `run_calibration_batch.py` enables it, writing
  `<output>.diagnostics.root` next to each calibration table.

- **`EcalEventBuilder`** (`src/components/EcalEventBuilder.cpp`,
  `include/k4SiWEcalReco/{EventBuilder,CalibrationTables,PadMapGeometry}.h`)
  — reconstructs physics events from a `siwecaldecoded` tree: BCID
  clustering across the 15 slabs (with 12-bit BCID-overflow unwrapping,
  independent of the raw2root converter's own `corrected_bcid` branch, which
  has known per-frame-reset quirks — see `EventBuilder.h`), per-channel
  best-SCA hit selection, pedestal/MIP calibration application (reads the
  `CalibrationTables.h`-format text files above), and pad-map + slab-z/X0
  geometry lookup. Writes the `ecal` tree in the **exact** schema
  `siwecal_eventbuilder/root_io.py::EcalWriter` already produces, so
  `EcalToEDM4hep` reads it unchanged. `threshold_dac` is read directly from
  the `siwecaldecoded` tree's `thresholdDac` branch (written by
  `EcalRawDecoder` above), not re-parsed from `Run_Settings.txt`.
  **Verified against the Python `siwecal_eventbuilder` reference**: exact
  event-count match, exact `(spill,bcid)` key-set match, and exact
  hit-by-hit field match (position, energy, mask) on sampled events with up
  to ~1300 hits each.

All three are plain `Gaudi::Algorithm` (not `k4FWCore::Producer`/
`Transformer`): each is a batch reduction over an entire run whose output
row count isn't known ahead of time, and none of them produce an
EDM4hep/podio collection — they read/write plain ROOT trees or ASCII tables.
All work happens in `initialize()`; `execute()` is a no-op.

### PID / EDM4hep

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
cmake -S gaudi_source -B gaudi_source/build -DCMAKE_BUILD_TYPE=Release
cmake --build gaudi_source/build -j4

export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
```

### Calibration batch driver: `gaudi_jobs/calibration/run_calibration_batch.py`

Chains `EcalRawDecoder` → `PedestalMipCalibrator` (two `k4run` invocations, so
a run can be re-calibrated without re-decoding its raw binary). Pedestal and
MIP are computed from the **same** data-taking run(s) — the reference tool's
event loop (`DecodedSLBAnalysis::NSlabsAnalysis`) fills pedestal histograms
from `hitbit==0` samples and MIP histograms from `hitbit==1` samples in one
pass over the same input file, and the deployed `MuonCalib_it2/` area always
pairs a `Pedestal_*` and a `MIP_*` file for the same run number. So there is
a single `--runs`, used as input for **both** stages, and statistics from
every run given are pooled together into one combined pedestal calibration
and one combined MIP calibration (maximizes usable statistics per ASU;
`PedestalMipCalibrator` already accumulates over every `InputFiles` entry
into shared histograms for both modes — see `forEachSelectedChannel` in
`PedestalMipCalibrator.cpp`).

`--runs` takes a comma-separated list of run **folder paths** (bare run
names are also accepted, resolved against `--raw-base`, default
`.../Data/rundata/`). Before decoding anything, every run's
`Run_Settings.txt` is read to check they all share the same `ThresholdDAC`
— that shared value fixes the output `th<N>` subdirectory automatically
(printed as it's checked); mixing thresholds without an explicit `--th`
override aborts with a run→threshold table instead of guessing:

```bash
python gaudi_jobs/calibration/run_calibration_batch.py \
    --runs TB2026CERN_run_000060,TB2026CERN_run_000061,TB2026CERN_run_000062

# explicit folder paths, pedestal-only
python gaudi_jobs/calibration/run_calibration_batch.py \
    --runs /eos/.../rundata/TB2026CERN_run_000060 --pedestal-only
```

Output layout mirrors `MuonCalib_it2/{pedestals,mips}/th<N>/` under the
output base (default `settings.yml calib_dir()/MuonCalib_gaudi`). A single
run produces `Pedestal_<run>_<gain>gain.txt`/`MIP_pedestalsubmode1_<run>_
<gain>gain.txt` (matching the individual per-run files already deployed);
pooling more than one run produces
`Pedestal_TB2026CERN_run_000th<N>_<gain>gain.txt`/`MIP_pedestalsubmode1_
TB2026CERN_run_000th<N>_<gain>gain.txt` — the same `run_000th<N>` literal
already used by the deployed cumulative MIP files and expected by
`siwecal_eventbuilder.cli.resolve_muon_calib_files(th)`.

Use **`gaudi_jobs/calibration/list_run_thresholds.py`** to see which runs share a
threshold without opening each `Run_Settings.txt` by hand — it groups every
run under `--raw-base` (or an explicit `--runs` list) by `ThresholdDAC` and
prints a ready-to-paste `--runs` value per group (`--output FILE.txt` also
exports the same report to a text file):

```bash
python gaudi_jobs/calibration/list_run_thresholds.py
python gaudi_jobs/calibration/list_run_thresholds.py --runs TB2026CERN_run_000060,TB2026CERN_run_000004 --output th_map.txt
```

Fit-robustness tuning (`--max-nhit`/`--nslabs-hit` are existing event-selection
flags; the new outlier-handling ceilings are set via `k4run` properties on
`PedestalMipCalibrator` if the defaults documented above need overriding for a
specific run — `PedestalMaxMeanAdc`, `MipFallbackMaxAdc`).

**SAFETY**: run folders under `--raw-base` (default `.../Data/rundata/`) are
read-only; decoded `siwecaldecoded.root` files always go to
`--converted-dir` (default `.../Data/rundata_converted_gaudi/`), enforced by
an `_assert_not_under()` guard that refuses to write inside a run's raw
folder. `list_run_thresholds.py` is read-only (only opens `Run_Settings.txt`
files, never touches raw binaries or calls `k4run`).

### Full pipeline: raw2root → event building → PID

`gaudi_jobs/run_full_pipeline_batch.py` (generic) and
`gaudi_jobs/run000013/run_full_pipeline_000013.py` (concrete example: 52 GeV
electrons, run `TB2026CERN_run_000013`) run the **entire** chain from raw
binary to EDM4hep. Both are plain Python scripts (run with `python`, not
`k4run`) that orchestrate **two** `k4run` processes:

1. `options/run_raw2root.py`, **one process per raw chunk**, into
   `<converted-dir>/<run>/chunks/` — followed by a health check that the
   acquisitions add up (`gaudi_jobs/decode_chunks.py`). This used to be a single
   `run_raw2root_and_eventbuilder.py` process decoding the whole run at once;
   that silently dropped ~75% of the acquisitions on some runs, so the steering
   file is gone and `run_raw2root.py` warns if handed more than one chunk.
2. `options/run_event_builder.py` (`EcalEventBuilder`), chaining those chunks via
   `EVBLD_INPUT_FILES`.
3. `options/run_pid.py` (`TopAlg=[EcalToEDM4hep, EcalPidTransformer]`) — a
   **separate** process, because these are `k4FWCore::Producer`/
   `Transformer` components that need `EvtMax` (the `ecal` tree's entry
   count) fixed *before* the process starts, and that count is only known
   after stage 1 has run. The driver script reads it from the just-written
   `ecal_<run>.root` in between the two `k4run` calls.

```bash
python gaudi_jobs/run000013/run_full_pipeline_000013.py

python gaudi_jobs/run_full_pipeline_batch.py --run TB2026CERN_run_000013 --th 230
```

Calibration files are auto-resolved via
`siwecal_eventbuilder.cli.resolve_muon_calib_files(th)` (the same lookup the
Python pipeline already uses) unless `--pedestal-file`/`--mip-file` are given
explicitly.

### Batch driver (recommended): `gaudi_jobs/run_pid_batch.py`

Resolves inputs exactly like `siwecal_validation` (`--run/--file/--all/--point/--cfg`),
runs the Gaudi stage per `ecal_<run>.root`, then applies the cuts and writes the
final output(s):

```bash
# physics mode (default): hit-MIP cut >= 0.5, energy>0, no mip05_/mip1_ blocks
python gaudi_jobs/run_pid_batch.py --run TB2026CERN_run_000007

# both formats + an event cut (cuts use the same flags as siwecal_validation)
python gaudi_jobs/run_pid_batch.py --all --format both --nhit-min 20

# validation/visualizer mode: no hit cut, compute the mip05_/mip1_ slider blocks
python gaudi_jobs/run_pid_batch.py --run TB2026CERN_run_000007 --validation
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
Outputs are named `ecal_<label>.edm4hep.root` and `ecal_<label>.valtree.root`
under `settings.yml` `pid_dir` (or `--outdir`).

### Ejemplo concreto: `gaudi_jobs/run000013/`

Two example scripts for TB2026CERN_run_000013, covering different starting points:

- `steer_run000013.py` — **PID-only**, `k4run` steering file, assumes
  `ecal_TB2026CERN_run_000013.root` already exists (produced beforehand, by
  either `siwecal_eventbuilder` or `EcalEventBuilder`):
  ```bash
  k4run gaudi_jobs/run000013/steer_run000013.py
  ```
- `run_full_pipeline_000013.py` — the **entire chain from raw binary**
  (see "Full pipeline" above), a plain Python script:
  ```bash
  python gaudi_jobs/run000013/run_full_pipeline_000013.py
  ```

### Single file via `k4run` (low level)

Each Gaudi stage also has a low-level, env-var-driven `k4run` steering file
under `gaudi_source/options/` for running it standalone (no batch-driver cut
application or output-file naming):

| Stage | Steering file | Key env vars |
|---|---|---|
| raw2root (ONE chunk per process) | `run_raw2root.py` | `RAW_FILES`, `RAW2ROOT_OUT`, `RAW2ROOT_RUN_SETTINGS_FILE` |
| pedestal/MIP calibration | `run_pedestal_mip.py` | `CALIB_INPUT_FILES`, `CALIB_MODE`, `CALIB_GAIN`, `CALIB_OUTPUT_PEDESTAL_FILE`/`CALIB_OUTPUT_MIP_FILE` |
| event building | `run_event_builder.py` | `EVBLD_INPUT`, `EVBLD_OUTPUT`, `EVBLD_PEDESTAL_FILE`, `EVBLD_MIP_FILE`, `EVBLD_PADMAP_DEFAULT`, `EVBLD_SLAB_Z_FILE` |
| PID/EDM4hep | `run_pid.py` | `ECAL_FILE`, `ECAL_PID_OUT`, `ECAL_HIT_MIP_CUT` (`<0` disables), `ECAL_MIP_THRESHOLDS` (`""` = no variant blocks) |

```bash
ECAL_FILE=/path/ecal_<run>.root ECAL_PID_OUT=out.root \
ECAL_HIT_MIP_CUT=0.5 ECAL_MIP_THRESHOLDS= \
    k4run gaudi_source/options/run_pid.py
```
(`k4run` alone does not apply the event-level cuts — use `gaudi_jobs/run_pid_batch.py` for
the cut-passing, correctly-named outputs.)

## Verification

### `EcalEventBuilder` vs. `siwecal_eventbuilder` (Python reference)

Run on the same `siwecaldecoded.root`, same calibration/pad-map/geometry
inputs: **exact event-count match** (3019/3019), **exact `(spill,bcid)`
key-set match** (0 in either direction), and **exact hit-by-hit field match**
(slab/chip/channel/sca/hg/lg/energy/x/y/z/X0/ismasked, all rounded to 3
decimals) on sampled events with 347–1335 hits each. Independently confirmed:
`hit_z` follows `mappings/slab_z_positions.yml`'s negative-z convention, and
`hit_X0` matches the expected cumulative-W/X0 calculation exactly for every
sampled slab.

### PID/EDM4hep parity vs `metrics.py`

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
