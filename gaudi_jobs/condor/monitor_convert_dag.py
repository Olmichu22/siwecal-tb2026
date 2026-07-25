#!/usr/bin/env python3
"""
Watchdog for a running CONVERT DAG (generate_convert_eudaq_dag.py output).

One invocation = one check. It is meant to be called repeatedly (by hand, by a
loop, or by a scheduled agent) while a large conversion runs for hours:

  1. reports DAG progress (Done / Queued / Failed, from the dagman.out status
     line) and the queue state of this DAG's jobs;
  2. sweeps the AFS ``logs/`` directory -- deletes the per-job ``.log`` of jobs
     that are no longer in the queue (and ONLY the ``.log``: a surviving
     ``.out``/``.err`` belongs to a FAILED node, since the DAG's POST script
     removes them on success, and it is the only record of the failure).
     DAGMan does NOT read the per-job ``.log``: it takes its events from
     ``<dag>.nodes.log`` (verified in the
     dagman.out: "Default node log file is ... .nodes.log"), the per-job ``log =``
     copy is informational only. This matters because the AFS logs/ directory
     hits a hard directory-entry limit long before the volume quota does, and
     new jobs then fail to write their log at all ("errno 27 File too large");
  3. releases held jobs -- cgroup-memory holds (code 34) are already released by
     the submit file's ``periodic_release``, but a hold caused by a full logs/
     directory only clears once step 2 has run, so it is retried here;
  4. counts the decoded chunks actually on disk per run, so progress is measured
     on the output, not only on Condor's bookkeeping.

Read-only unless ``--sweep-logs``/``--release`` are given (both default ON;
use ``--no-sweep-logs``/``--no-release`` for a pure status check).

Usage::

    python gaudi_jobs/condor/monitor_convert_dag.py \\
        --dag-dir gaudi_jobs/condor/generated/convert_all_new
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", ".."))

_STATUS_RE = re.compile(
    r"Done\s+Pre\s+Queued\s+Post\s+Ready\s+Un-Ready\s+Failed", re.I)
_TS_RE = re.compile(r"^\s*\d\d/\d\d/\d\d\s+\d\d:\d\d:\d\d\s*")


def _run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True).stdout


def dag_progress(dag_dir):
    """Parse the last 'Done Pre Queued Post Ready Un-Ready Failed' table from
    the dagman.out. Returns a dict, or None if the DAG has not written one yet."""
    outs = [f for f in os.listdir(dag_dir) if f.endswith(".dag.dagman.out")]
    if not outs:
        return None, None
    path = os.path.join(dag_dir, sorted(outs)[0])
    with open(path, errors="replace") as fh:
        lines = fh.readlines()
    for i in range(len(lines) - 1, -1, -1):
        if _STATUS_RE.search(lines[i]):
            # The counts sit two lines below the header (a '===' rule between).
            # Every dagman.out line starts with 'MM/DD/YY HH:MM:SS ' -- strip it
            # before pulling numbers out, or the date digits become counts.
            for j in range(i + 1, min(i + 4, len(lines))):
                body = _TS_RE.sub("", lines[j])
                nums = re.findall(r"\b\d+\b", body)
                if len(nums) >= 7:
                    keys = ["done", "pre", "queued", "post", "ready", "unready", "failed", "futile"]
                    return dict(zip(keys, (int(n) for n in nums))), path
    return None, path


def queue_state(dag_dir):
    """condor_q rows belonging to this DAG's submit directory."""
    out = _run(["condor_q", "-json", "-attributes",
                "ClusterId,ProcId,JobStatus,HoldReasonCode,HoldReason,Iwd,Cmd,UserLog"])
    try:
        rows = json.loads(out) if out.strip() else []
    except json.JSONDecodeError:
        return []
    return [r for r in rows if os.path.abspath(r.get("Iwd", "")) == os.path.abspath(dag_dir)]


def sweep_logs(dag_dir, live_ids, min_age_s=900):
    """Delete per-job .log/.out/.err whose job is gone from the queue.

    ``live_ids`` is the set of log BASENAMES belonging to jobs still in the
    queue, taken from each job ad's own ``UserLog`` attribute -- so the
    protection is per PROCESS, not per run. That distinction is the whole point:
    a 20-job run keeps one live .log while the other 19 are finished, and
    protecting the run as a whole would hold those 19 (plus their .out/.err) on
    AFS for the many hours the run takes to drain.

    Still conservative on age: a file younger than ``min_age_s`` is left alone,
    since a just-submitted job may not be in the condor_q snapshot yet. The
    DAG's own *.dag.* files are never touched -- DAGMan is reading those.
    """
    log_dir = os.path.join(dag_dir, "logs")
    if not os.path.isdir(log_dir):
        return 0, 0
    now = time.time()
    removed = freed = 0
    for name in os.listdir(log_dir):
        # ONLY the per-job .log. The .out/.err are already deleted by the DAG's
        # POST script the moment a node SUCCEEDS, so any that survive belong to
        # a node that FAILED -- they are the only record of why, and sweeping
        # them (this script did, once) destroys the evidence needed to fix the
        # run. Leave them for a human.
        if not name.endswith(".log"):
            continue
        if name in live_ids:
            continue
        path = os.path.join(log_dir, name)
        try:
            st = os.stat(path)
        except OSError:
            continue
        if now - st.st_mtime < min_age_s:
            continue
        try:
            os.remove(path)
            removed += 1
            freed += st.st_size
        except OSError:
            pass
    return removed, freed


def decoded_chunks(converted_dir, runs):
    """(run -> n decoded chunk files on disk) for the runs of this DAG."""
    counts = {}
    for run in runs:
        cdir = os.path.join(converted_dir, run, "chunks")
        try:
            counts[run] = sum(1 for f in os.listdir(cdir) if f.endswith(".root"))
        except OSError:
            counts[run] = 0
    return counts


def expected_chunks(dag_dir):
    """(run -> n chunks the DAG intends to write), read from its grouplists."""
    groups = os.path.join(dag_dir, "groups")
    expected = {}
    if os.path.isdir(groups):
        for name in sorted(os.listdir(groups)):
            run = name.rsplit("_g", 1)[0]
            with open(os.path.join(groups, name)) as fh:
                expected[run] = expected.get(run, 0) + sum(1 for ln in fh if ln.strip())
    for name in os.listdir(dag_dir):
        if name.startswith("convert_lcio_") and name.endswith(".sub"):
            expected.setdefault(name[len("convert_lcio_"):-len(".sub")], 1)
    return expected


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dag-dir", required=True, help="Submit dir of the running CONVERT DAG.")
    p.add_argument("--converted-dir", default=None,
                   help="Where the decoded chunks land (default: read from the DAG's groups/).")
    p.add_argument("--sweep-logs", dest="sweep", action="store_true", default=True)
    p.add_argument("--no-sweep-logs", dest="sweep", action="store_false")
    p.add_argument("--min-log-age", type=int, default=600,
                   help="Only sweep log files older than this many seconds (default 600). The age guard "
                        "only covers jobs not yet visible in the condor_q snapshot; live jobs are "
                        "protected exactly, by their UserLog.")
    p.add_argument("--max-log-entries", type=int, default=2000,
                   help="If logs/ still holds more than this many entries after the normal sweep, sweep "
                        "again ignoring --min-log-age (live logs stay protected). Guards the AFS "
                        "directory-entry limit, which bites long before the volume quota does.")
    p.add_argument("--release", dest="release", action="store_true", default=True)
    p.add_argument("--no-release", dest="release", action="store_false")
    p.add_argument("--verbose-runs", action="store_true",
                   help="List every run's decoded/expected chunk count, not just the incomplete ones.")
    args = p.parse_args(argv)

    dag_dir = os.path.abspath(args.dag_dir)
    print(f"=== {time.strftime('%Y-%m-%d %H:%M:%S')}  {dag_dir}")

    prog, dagman_out = dag_progress(dag_dir)
    if prog:
        print(f"DAG: done={prog['done']}/{sum(prog.values()) - prog['pre'] - prog['post']} "
              f"queued={prog['queued']} ready={prog['ready']} "
              f"unready={prog['unready']} failed={prog['failed']} futile={prog.get('futile', 0)}")
    else:
        print("DAG: no status line yet in the dagman.out")

    rows = queue_state(dag_dir)
    status_names = {1: "idle", 2: "running", 3: "removing", 4: "completed", 5: "HELD", 6: "xfer", 7: "suspended"}
    tally = {}
    held = []
    for r in rows:
        tally[status_names.get(r.get("JobStatus"), "?")] = tally.get(status_names.get(r.get("JobStatus"), "?"), 0) + 1
        if r.get("JobStatus") == 5:
            held.append(r)
    print("queue: " + (", ".join(f"{k}={v}" for k, v in sorted(tally.items())) if tally else "no jobs"))
    for r in held[:10]:
        print(f"  HELD {r['ClusterId']}.{r['ProcId']} code={r.get('HoldReasonCode')} "
              f"{str(r.get('HoldReason'))[:120]}")

    # Exactly which log files are live: each job ad carries its own UserLog
    # path, so a finished process's log is swept even while sibling processes
    # of the same run are still running.
    live_logs = {os.path.basename(r["UserLog"]) for r in rows if r.get("UserLog")}
    if rows and not live_logs:
        # No UserLog in the ads (older schedd): fall back to protecting the
        # whole run rather than risk deleting a live job's log.
        for r in rows:
            m = re.search(r"(TB2026CERN_[a-z_]*run_\d+)", json.dumps(r))
            if m:
                live_logs.update(f"convert_{m.group(1)}_{p}.log" for p in range(200))
                live_logs.add(f"convert_lcio_{m.group(1)}.log")

    if args.sweep:
        log_dir = os.path.join(dag_dir, "logs")
        n, freed = sweep_logs(dag_dir, live_logs, args.min_log_age)
        left = len(os.listdir(log_dir)) if os.path.isdir(log_dir) else 0
        if left > args.max_log_entries:
            n2, freed2 = sweep_logs(dag_dir, live_logs, 0)
            n, freed = n + n2, freed + freed2
            left = len(os.listdir(log_dir))
            print(f"logs: over {args.max_log_entries} entries -- swept again ignoring age")
        print(f"logs: swept {n} file(s), {freed/1e6:.1f} MB freed; {left} entries left in logs/")

    if args.release and held:
        rel = _run(["condor_release", "-constraint",
                    f'JobStatus == 5 && Iwd == "{dag_dir}"'])
        print("release: " + rel.strip().replace("\n", " "))

    exp = expected_chunks(dag_dir)
    conv = args.converted_dir
    if conv is None:
        conv = _converted_dir_from_dag(dag_dir)
    if conv:
        got = decoded_chunks(conv, exp)
        total_exp = sum(exp.values())
        total_got = sum(min(got[r], exp[r]) for r in exp)
        pct = 100.0 * total_got / total_exp if total_exp else 0.0
        done_runs = sum(1 for r in exp if got[r] >= exp[r])
        print(f"chunks: {total_got}/{total_exp} ({pct:.1f}%) on disk under {conv}")
        print(f"runs complete: {done_runs}/{len(exp)}")
        for r in sorted(exp):
            if args.verbose_runs or got[r] < exp[r]:
                print(f"  {r}: {got[r]}/{exp[r]}")
    return 0


def _converted_dir_from_dag(dag_dir):
    """Recover the output base from any grouplist row ('<raw> <out_decoded>')."""
    groups = os.path.join(dag_dir, "groups")
    if os.path.isdir(groups):
        for name in sorted(os.listdir(groups)):
            with open(os.path.join(groups, name)) as fh:
                for line in fh:
                    parts = line.split()
                    if len(parts) == 2:
                        # <base>/<run>/chunks/chunk_NNNN.root -> <base>
                        return os.path.dirname(os.path.dirname(os.path.dirname(parts[1])))
    return None


if __name__ == "__main__":
    sys.exit(main())
