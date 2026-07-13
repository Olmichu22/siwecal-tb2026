# Calibration Condor stage — logs, testing runs, and memory

This folder generates and drives the **CALIBRATION** DAGMan pipeline
(`fill` per decoded chunk → `merge` per run → `merge` per threshold → `fit`).
The generator is [`generate_dag.py`](generate_dag.py); see its module docstring
and [`../README.md`](../README.md) for the full two-stage design.

This note covers two operational things that bit us in July 2026:
per-job logs filling up AFS, and merge jobs running out of memory.

## Condor per-job logs and the AFS directory limit

Each DAG node writes three files into `<dag-dir>/logs/` (on AFS):

| file | written by | read by DAGMan? |
|------|------------|-----------------|
| `<JOB>.log` | Condor (job events) | **yes** — leave it alone while the node is active |
| `<JOB>.out` | Condor (job stdout) | no |
| `<JOB>.err` | Condor (job stderr) | no |

A calibration run has **thousands** of `fill` nodes (e.g. run_000085 has 2100+
chunks). That many files in one directory hits an **AFS limit**: an AFS
directory object maxes out around ~2 MB / ~64k entries, and once it is full,
creating a new **long-named** file fails with:

```
Transfer output files failure ... writing to file .../logs/fill_..._0195.out: (errno 27) File too large
```

This is **not** a disk-quota problem (`fs lq` can show plenty free) — it is the
directory itself being full. The job then goes on **HOLD** and the DAG stalls.

### Automatic fix: `.out`/`.err` deleted on success (default)

By default the generator attaches a DAGMan `SCRIPT POST` to every
`fill`/`merge_group`/`merge_run`/`fit` node
([`cleanup_job_logs.sh`](cleanup_job_logs.sh)) that **deletes that node's
`.out` and `.err` as soon as the node succeeds**. So the `logs/` dir only ever
holds the `.log` files plus the `.out/.err` of *currently-running* and
*failed* nodes — orders of magnitude fewer entries.

Guarantees:
- **Failed nodes keep their `.out/.err`** (the POST only deletes on success) —
  so you can still debug failures.
- The `.log` files are never touched.
- `RETRY` semantics are unchanged: the POST propagates the node's real exit
  code, so a failed node still retries and is never mis-marked as succeeded.
- The final `merge_threshold_th<N>` node keeps its own POST
  ([`check_merged_histograms.sh`](check_merged_histograms.sh)); it is a single
  node, so its logs are irrelevant to the directory limit.

### Testing runs: keep all logs with `--keep-job-logs`

When you are debugging a **small** run and want every job's stdout/stderr
around afterwards, generate the DAG with `--keep-job-logs`. That omits the
cleanup POST entirely, so **all `.out/.err` are kept**:

```bash
python gaudi_jobs/calibration/condor/generate_dag.py \
    --runs TB2026CERN_run_000063 \
    --dag-dir gaudi_jobs/calibration/condor/generated/th230_test \
    --keep-job-logs
condor_submit_dag gaudi_jobs/calibration/condor/generated/th230_test/calibration_th230.dag
```

A small test run has few nodes, so it is nowhere near the directory limit —
keeping all logs is safe there.

> This is a **generation-time** choice: it is baked into the `.dag` when it is
> written. To switch a run between keeping and deleting logs, regenerate the
> DAG. (Small tests get their own `--dag-dir` anyway, so this is natural.)

### If a `logs/` dir is already full

Delete the `.out`/`.err` of completed jobs (never read by DAGMan) to free
directory slots, then release the held jobs:

```bash
find <dag-dir>/logs -maxdepth 1 -name '*.out' -delete
find <dag-dir>/logs -maxdepth 1 -name '*.err' -delete
condor_release <user>
```

Only delete `.log` files of nodes that are **not** currently queued/running.

## Merge job memory (OOM → HOLD)

Large `hadd` merges were observed using **~11.8 GB** while requesting only
4 GB. On the batch farm that trips the cgroup memory limit: the job is put on
**HOLD** with `HoldReasonCode == 34`. A held job is **not** retried by DAGMan's
`RETRY`, so the DAG stalls.

The generated `.sub` files handle this automatically (see `_write_subs` in
[`generate_dag.py`](generate_dag.py)):

- **request_memory grows with each attempt**:
  `request_memory = ifthenelse(NumJobStarts <= 0, base, base * (NumJobStarts + 1))`
  — `base` on the first start, `2×base, 3×base, …` on each subsequent one.
- **OOM holds are auto-released**:
  `periodic_release = (JobStatus == 5) && (HoldReasonCode == 34) && (NumJobStarts < 5)`
  — an OOM-held job is re-run (up to 5 starts) and comes back with more memory,
  self-healing without manual intervention.
- The **merge base default is 8 GB** (`--merge-request-memory 8000`), so
  attempt 2 asks for 16 GB — comfortably above the ~11.8 GB observed.

Tune the base per job type with `--fill-request-memory` / `--merge-request-memory`
/ `--fit-request-memory` if a dataset needs more.
