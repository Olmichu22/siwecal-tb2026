#!/usr/bin/env bash
# List ThresholdDAC for EVERY run under the raw data base directory (no
# fixed run list -- see check_selected_run_thresholds.sh for that case) and
# save the grouped report to a text file. Read-only: wraps
# list_run_thresholds.py, which only opens each run's Run_Settings.txt,
# never touches raw binaries or calls k4run.
#
# Usage:
#   gaudi_jobs/calibration/list_all_run_thresholds.sh                        # writes gaudi_jobs/calibration/run_thresholds_all.txt
#   gaudi_jobs/calibration/list_all_run_thresholds.sh --output other.txt      # custom output path
#   gaudi_jobs/calibration/list_all_run_thresholds.sh --raw-base /some/dir --output other.txt
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

DEFAULT_OUTPUT="${SCRIPT_DIR}/run_thresholds_all.txt"

# --output given here first, then "$@": argparse keeps the LAST occurrence
# of a store-type option, so any --output/--raw-base the caller passes in
# "$@" transparently overrides this default without extra bash logic.
exec python3 "${SCRIPT_DIR}/list_run_thresholds.py" --output "${DEFAULT_OUTPUT}" "$@"
