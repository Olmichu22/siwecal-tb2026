#!/usr/bin/env python
"""
List (or export) which ThresholdDAC each run was taken at, grouped so a
coherent ``--runs`` value for ``run_calibration_batch.py`` can be copy-pasted
straight out of the output -- no need to open each run's ``Run_Settings.txt``
by hand.

Two ways to pick the run list:

    # explicit list (same folder/name syntax as run_calibration_batch.py --runs)
    python gaudi_jobs/calibration/list_run_thresholds.py --runs r1,r2,...

    # auto-discover: every immediate subdirectory of --raw-base (default:
    # calib_run_utils.DEFAULT_RAW_BASE) is treated as a run
    python gaudi_jobs/calibration/list_run_thresholds.py
    python gaudi_jobs/calibration/list_run_thresholds.py --raw-base /some/other/rundata

Add ``--output FILE.txt`` to also write the same report to a text file.

Read-only: only opens Run_Settings.txt files, never touches raw binaries or
calls k4run.
"""
import argparse
import sys
from collections import defaultdict

from calib_run_utils import DEFAULT_RAW_BASE, discover_run_folders, parse_run_folder_list, run_threshold


def _build_report(run_folders):
    by_threshold = defaultdict(list)
    for run_name, raw_dir in run_folders:
        th = run_threshold(raw_dir)
        by_threshold[th].append((run_name, raw_dir))

    lines = []
    known_thresholds = sorted(th for th in by_threshold if th != -1)
    for th in known_thresholds:
        runs = by_threshold[th]
        lines.append(f"=== ThresholdDAC={th} ({len(runs)} run{'s' if len(runs) != 1 else ''}) ===")
        for run_name, raw_dir in runs:
            lines.append(f"  {run_name}  {raw_dir}")
        lines.append("--runs value:")
        lines.append(",".join(raw_dir for _, raw_dir in runs))
        lines.append("")

    if -1 in by_threshold:
        runs = by_threshold[-1]
        lines.append(f"=== ThresholdDAC unknown ({len(runs)} run{'s' if len(runs) != 1 else ''}, "
                      "Run_Settings.txt missing/unreadable) ===")
        for run_name, raw_dir in runs:
            lines.append(f"  {run_name}  {raw_dir}")
        lines.append("")

    return "\n".join(lines).rstrip("\n") + "\n"


def main(argv=None):
    p = argparse.ArgumentParser(description="List/export run<->ThresholdDAC equivalence, grouped by threshold.")
    p.add_argument("--runs", default=None,
                   help="Comma-separated run folder path(s) or bare run name(s) to inspect "
                        "(same syntax as run_calibration_batch.py --runs). "
                        "If omitted, auto-discovers every run under --raw-base.")
    p.add_argument("--raw-base", default=DEFAULT_RAW_BASE,
                   help=f"Base directory to auto-discover runs from when --runs is omitted. Default: {DEFAULT_RAW_BASE}")
    p.add_argument("--output", default=None, help="Also write the report to this text file.")
    args = p.parse_args(argv)

    if args.runs:
        run_folders = parse_run_folder_list(args.runs, default_base=args.raw_base)
    else:
        run_folders = discover_run_folders(args.raw_base)

    if not run_folders:
        raise SystemExit("ERROR: no runs to inspect")

    report = _build_report(run_folders)
    print(report, end="")

    if args.output:
        with open(args.output, "w") as fh:
            fh.write(report)
        print(f"\n[written] {args.output}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
