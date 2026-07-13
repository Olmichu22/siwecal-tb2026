# `gaudi_jobs/condor/` — raw → EDM4hep on the batch farm

[`generate_reco_dag.py`](generate_reco_dag.py) builds one DAGMan pipeline that
takes any number of runs, of any type (physics, muon, calibration), from raw
binary to EDM4hep:

```
CONVERT  one job per raw chunk, all in parallel   →  <RUN>/chunks/chunk_NNNN.root
   │                                                  (decoded, siwecaldecoded)
   ▼
RECO     one job per run                          →  <RUN>/ecal_<RUN>.root
         event building (chains the chunks)          <RUN>/ecal_<RUN>.edm4hep.root
         + PID/EDM4hep
```

Per run, the **output is still exactly one `ecal` file and one EDM4hep file** —
the thing analyses actually open. What disappears is the intermediate merged
`siwecaldecoded`.

```bash
python gaudi_jobs/condor/generate_reco_dag.py \
    --runs TB2026CERN_run_000012,TB2026CERN_run_000013 \
    --out-dir gaudi_jobs/condor/generated/physics

condor_submit_dag gaudi_jobs/condor/generated/physics/reco.dag
```

Calibration tables are resolved **per run** from that run's own `ThresholdDAC`
(read from its `Run_Settings.txt`), so runs taken at different thresholds can go
in the same submission. Low-gain tables are picked up when they exist; when they
don't, the run is still reconstructed, just without saturation recovery (the job
log says so).

## Why not the serial drivers

`gaudi_jobs/run_full_pipeline_*.py` decode every raw chunk of a run inside **one**
k4run process and write a single multi-GB `siwecaldecoded` file, which the event
builder then reads back. That file is pure overhead, and it costs:

* **Space.** It duplicates the decoded data on EOS — ~10 GB for one physics run,
  ~390 GB across the calibration runs — for something no analysis ever opens.
* **Time.** Decoding is serial. One physics run takes ~20 min of wall clock;
  300 parallel chunk jobs do it in a couple of minutes.
* **Reliability.** Closing a 10 GB file over XRootD times out often enough to be
  a real operational problem (it stalled `run_000013` three times in one night),
  and a failed `Close()` leaves a 0-byte file that fails the next stage.
* **A hard ceiling.** Past ~100 GB it is not even possible: ROOT's
  `TTree::fgMaxTreeSize` rolls the tree over into a second file, so **no single
  file can hold it**. `TB2026CERN_run_000004`'s 1369 chunks come to ~176 GB —
  `hadd` dies on it, and `TFileMerger` ran for 5 h without producing a byte.

The fix was to stop needing the merge: `EcalEventBuilder` takes `InputFiles` (a
list or a glob) and chains them with a `TChain`. Chunking was always a
parallelisation detail of CONVERT, never a property of the data — the calibrator
already read chunks this way.

The serial drivers still work and are still the simplest way to run one run
interactively; they are just not how you process a campaign.

## Layout on EOS

Two areas, split by *can this be recreated?*

```
Data/rundata_converted_gaudi/<RUN>/          # data — cannot be recreated
└── chunks/chunk_0000.root …                 decoded (CONVERT). Keep: the
                                             calibration Fill reads these, and
                                             so does every re-reco.

Reconstruction/<RUN>/                        # products — regenerate at will
├── ecal_<RUN>.root                          event-built (RECO)
├── ecal_<RUN>.edm4hep.root                  PID/EDM4hep (RECO)
└── energy_dist/ weighte/ moliere/ …         validation plots for this run
```

Everything below the event builder is a *product*: it is rebuilt whenever the
calibration or the geometry changes, so it lives outside `Data/` and the whole
`Reconstruction/` tree can be wiped without risking a byte of data. The
validation labels each sample by run name and defaults its output base to
`reconstruction_dir` (`settings.yml`), which is why a run's plots land beside the
very events they were made from. Its cross-run results tables
(`results_id0N.csv/.txt`) stay at the top of `Reconstruction/`, since they span
runs.

Both locations come from `siwecal_common.paths` (`reconstruction_dir()` /
`DEFAULT_CONVERTED_DIR`) — no script hard-codes them a second time.

## Restarting, retries, logs

* **Idempotent.** `convert.sh` skips a chunk whose output already exists, so a
  resubmission only does the missing work.
* **Rescue.** A failed DAG leaves a `*.rescue00N` file. Resume it with a plain
  `condor_submit_dag reco.dag` — **without `-f`**, which throws the rescue away
  and restarts from scratch.
* **OOM self-heals.** A job held on `HoldReasonCode 34` (cgroup memory) is
  released automatically with more memory on each retry.
* **Logs.** Each node's `.out`/`.err` is deleted when it succeeds, to protect the
  AFS `logs/` directory's entry limit (see
  [`../calibration/condor/README.md`](../calibration/condor/README.md) for why —
  it bit us). Failed nodes keep theirs. Pass `--keep-job-logs` to keep them all
  on a small test run.
