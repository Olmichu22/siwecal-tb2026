# angular_scan pipeline — runs 179-239 (th210)

Event-builds + validates (EDM4hep/PID) the angular-scan run range, from
chunks **already decoded** into `rundata_converted_gaudi/<RUN>/chunks/`
(see `gaudi_jobs/condor/generate_convert_eudaq_dag.py` and
`gaudi_jobs/decode_lcio_runs.py` for how those chunks were produced). This
pipeline does **no CONVERT step** — it only reads chunks that must already
exist.

Same one-file-per-run convention as `gaudi_jobs/run000166/` and
`gaudi_jobs/run000091/`: every run gets its own `steer_<RUN>.py` (event
builder) and `validate_<RUN>.py` (EDM4hep + PID), copied rather than
parameterized, on purpose — so any single run's calibration or options can be
tweaked without touching the others. All 55 runs here are **th210**, so every
steering file resolves the same `calibration/MuonCalib_gaudi/*/th210/` tables.

## Runs

179-239 inclusive (61 run numbers), minus 6 excluded:

- **197, 217, 218, 219, 220**: no raw or converted data exists for these run
  numbers at all (no `rundata/TB2026CERN_..._run_0002NN` folder under either
  naming) -- a gap in the run numbering, not a filter.
- **225**: a `chunks/chunk_0000.root` exists but holds 0 entries -- this is
  one of the 5 runs `gaudi_jobs/condor/classify_eudaq_routes.py` could not
  resolve during the eudaq conversion (`route: "failed"` in
  `gaudi_jobs/condor/generated/convert_eudaq/route_map.json`: 0 entries via
  the raw decoder and no `Data/eudaq/ROC_run_000225_tp.slcio` to fall back
  to). Its input simply doesn't exist yet; re-run
  `classify_eudaq_routes.py --runs TB2026CERN_eudaq_run_000225` if the raw or
  LCIO source data ever turns up, then regenerate its steer/validate pair
  from another run's pair as a template.

That leaves **55** run pairs (`steer_*.py` / `validate_*.py`), all
`TB2026CERN_eudaq_run_0001NN`.

## Run one

```bash
source setup.sh
export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH

k4run gaudi_jobs/angular_scan/steer_TB2026CERN_eudaq_run_000179.py
k4run gaudi_jobs/angular_scan/validate_TB2026CERN_eudaq_run_000179.py
```

Output goes to `rundata_converted_gaudi/<RUN>/events/` (next to that run's
`chunks/`, not the shared `Reconstruction/` area):
`ecal_<RUN>.root` (event builder) then `ecal_<RUN>.edm4hep.root` (PID/EDM4hep).

## Run all 56

```bash
source setup.sh
export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH

bash gaudi_jobs/angular_scan/run_all.sh
```

Sequential, local, idempotent (skips a run whose `ecal_<RUN>.edm4hep.root`
already exists unless `--force`). For farm-scale parallelism instead, use
`gaudi_jobs/condor/generate_reco_dag.py --runs <comma-separated run names>
--out-dir ... ` (its CONVERT stage will just verify the chunks are already
there and fall straight through to RECO/VALIDATE — see that script's
`--convert-only` note).

## Adapt to a different run / threshold

Copy one `steer_*.py`/`validate_*.py` pair, change `_RUN` (and `_TH` if the
new run isn't th210 — check its calibration tables exist under
`calibration/MuonCalib_gaudi/{pedestals,mips,anchor}/th<N>/` first).
