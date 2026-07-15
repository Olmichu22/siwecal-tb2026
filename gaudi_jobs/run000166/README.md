# Example pipeline — run000166 (th210, electrons)

Event-builds `TB2026CERN_eudaq_run_000166` from its already-converted chunks into
an `ecal` tree, using this project's th210 calibration
(`calibration/MuonCalib_gaudi/{pedestals,mips,anchor}/th210/`: pedestal + MIP +
LG→HG anchor line, all th210's own).

## Run it

```bash
source setup.sh
export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH

k4run gaudi_jobs/run000166/steer_run000166.py
```

Output: `ecal_TB2026CERN_eudaq_run_000166.root` in this folder (override with
`EVBLD_OUTPUT=/path/out.root`).

Needs the chunks to exist first. If they don't, convert once:

```bash
python gaudi_jobs/decode_lcio_runs.py --runs 166      # EUDAQ/LCIO run
# (raw SL-native runs instead: python gaudi_jobs/run_convert_batch.py --runs <run>)
```

## Adapt to multiple files

- **Many chunks of the same run:** nothing to do — the steer globs
  `chunks/chunk_*.root`, so it already reads every chunk of the run at once.
- **A different run / threshold:** copy this steer, and change `_RUN` and `_TH`
  at the top (one copy per run, on purpose — see run000091 for the th220 twin).
  Everything else (calibration, anchor, switch) is resolved from `_TH`.
