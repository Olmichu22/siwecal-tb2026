# `gaudi_jobs/calibration/`

Tools to generate pedestal/MIP calibration tables (`Pedestal_*.txt`,
`MIP_pedestalsubmode1_*.txt`) from raw SiW-ECAL test-beam data. Two ways to
run the same underlying `PedestalMipCalibrator` Gaudi component:

- **`run_calibration_batch.py`** — sequential, single process. Good for a
  handful of runs or a quick check.
- **`condor/`** — HTCondor DAGMan pipeline (fill per raw chunk → merge per
  run → merge per threshold → fit). Good for pooling a whole threshold's
  worth of runs, including very large ones (some individual runs here are
  1000+ raw chunk files / tens of GB), without one long fragile sequential
  process or risking a Condor slot's memory/walltime limits.

Rule of thumb: few runs / iterating on parameters → `run_calibration_batch.py`.
A full threshold group, especially with huge runs in it → `condor/`.

## Background: pedestal and MIP always share the same input runs

The reference tool's event loop
(`SiWECAL-TB-analysis/SLBperformance/DecodedSLBAnalysis.cc::NSlabsAnalysis`)
fills the pedestal histograms from samples where `hitbit==0` and the MIP
histograms from samples where `hitbit==1`, in **one pass over one input
file** — there is no separate "pedestal run" vs. "MIP run". Both tools here
take a single `--runs` list used as input for both calibrations, and pool
statistics from every run given into one combined pedestal table and one
combined MIP table (maximizes usable statistics per ASU).

## Finding which runs share a threshold

Calibration tables are organised per `ThresholdDAC` (`pedestals/th<N>/`,
`mips/th<N>/`) — mixing runs from different thresholds into one calibration
doesn't make physical sense, so both `run_calibration_batch.py` and the
Condor pipeline require every given run to agree on one `ThresholdDAC`
(read from each run's `Run_Settings.txt`, checked *before* touching any raw
data) and abort with a run→threshold table if they don't.

To find that grouping without opening each `Run_Settings.txt` by hand:

```bash
# every run under the default raw-data area, grouped by threshold
python gaudi_jobs/calibration/list_run_thresholds.py

# just a specific list of runs/folders
python gaudi_jobs/calibration/list_run_thresholds.py --runs TB2026CERN_run_000060,TB2026CERN_run_000004

# same as the first, but written to a .txt report too
gaudi_jobs/calibration/list_all_run_thresholds.sh
```

`check_selected_run_thresholds.sh` is a small checked-in wrapper recording
a *specific* run selection (edit its `RUNS` array) rather than retyping a
long `--runs` list each time — see the file for the current selection.

Both tools print, for each threshold, a ready-to-paste `--runs` value
(full folder paths, comma-separated) for the next two sections.

**Naming quirk**: raw chunk files come in two shapes depending on trigger
mode — plain runs use `<run>_raw.bin`/`<run>_raw.bin_NNNN`, but
`..._eudaq_run_...` runs use `<run>.bin`/`<run>.bin_NNNN` (no `_raw`
infix). `calib_run_utils.raw_chunk_files()` tries the plain pattern first
and falls back to the eudaq one automatically — this is why plain and
eudaq run names can be mixed freely in one `--runs` list.

## `run_calibration_batch.py` (sequential)

```bash
python gaudi_jobs/calibration/run_calibration_batch.py \
    --runs TB2026CERN_run_000060,TB2026CERN_run_000061,TB2026CERN_run_000062

# explicit folder paths, --th override, pedestal-only, cross-check plots
python gaudi_jobs/calibration/run_calibration_batch.py \
    --runs /eos/.../rundata/TB2026CERN_run_000060 --th 220 --pedestal-only --diagnostics
```

Chains `EcalRawDecoder` (raw2root, skipped if the decoded file already
exists) → `PedestalMipCalibrator` `Mode=Pedestal` and/or `Mode=Mip`, both
reading the same pooled `--runs` file list. See `gaudi_source/README.md`'s
"Calibration batch driver" section for the full flag reference (output
naming convention, fit-robustness tuning, `--diagnostics`).

## `condor/` (fill → merge → merge → fit)

Two **separate** entry points, run as two separate manual steps — Convert
is the *only* stage that ever opens a file under the read-only raw data
area; Calibration only ever reads already-decoded files from EOS scratch:

1. **Convert** (`condor/generate_convert_jobs.py`): one independent Condor
   job per raw chunk file, running `EcalRawDecoder` on just that chunk.
   No dependencies between jobs, so no DAG is needed — a single
   multi-queue `.sub` submitted directly.
   ```bash
   python gaudi_jobs/calibration/condor/generate_convert_jobs.py \
       --runs <paths from list_run_thresholds.py> --out-dir gaudi_jobs/calibration/condor/generated/th<N>
   condor_submit gaudi_jobs/calibration/condor/generated/th<N>/convert.sub
   ```
   Wait for all jobs to finish (`condor_q`) before moving on. Both this and
   the next step derive the decoded/histogram scratch layout from
   `--fill-scratch-dir` (default: a shared constant, see `condor_common.py`)
   — leave it at its default (or pass the same value to both) so Calibration
   finds what Convert wrote.

2. **Calibration** (`condor/generate_dag.py`): validates that Convert's
   decoded output exists for every requested run, then writes an HTCondor
   DAGMan pipeline: Fill (per chunk, histograms only, zero-suppressed) →
   merge per run (`hadd`) → merge per threshold (`hadd`) → Fit (the exact
   same fit/robustness/diagnostics code `run_calibration_batch.py` uses,
   just reading pre-built histograms instead of raw trees).
   ```bash
   python gaudi_jobs/calibration/condor/generate_dag.py \
       --runs <same paths> --dag-dir gaudi_jobs/calibration/condor/generated/th<N>
   condor_submit_dag gaudi_jobs/calibration/condor/generated/th<N>/calibration_th<N>.dag
   ```
   Can be re-run (regenerate + resubmit) as many times as needed to tune
   fit parameters (`--pedestal-max-mean-adc`, `--min-mip-integral`, ...)
   **without re-running Convert** — decoding only has to happen once.

Both generators only **write** files (`.sub`/`.sh`/`.dag`, file lists for
`hadd`, a `logs/` dir) and print the exact `condor_submit`/
`condor_submit_dag` command to run — they never submit anything
themselves. Submitting is a deliberate, explicit step: Convert alone can
mean thousands of jobs for a large threshold group.

Generated artifacts land in `condor/generated/` (gitignored — they're
per-invocation output, not code) and intermediate histograms/decoded files
in a dedicated EOS scratch area, `.../Data/calib_fill_scratch/{decoded,hist}/`.

## Safety rules (apply to every tool in this directory)

- `/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata/` is
  **READ-ONLY**. Nothing here writes or deletes anything under it.
- **Never delete anything under `/Data/`** (including the Convert/fill
  scratch area) — cleanup is the user's call, not automated by any script
  here.
- Decoded/intermediate output always goes to a separate location
  (`rundata_converted_gaudi/` for `run_calibration_batch.py`,
  `calib_fill_scratch/` for the Condor pipeline), never mixed with the raw
  data or with each other.
