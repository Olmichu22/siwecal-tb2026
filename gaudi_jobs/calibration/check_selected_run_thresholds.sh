#!/usr/bin/env bash
# Threshold check for a specific run selection, kept here as a record of
# which runs were inspected/run together (rather than passing the list on
# the command line each time). Edit the RUNS array below to change the
# selection. Wraps list_run_thresholds.py -- read-only, only opens each
# run's Run_Settings.txt, never touches raw binaries or calls k4run.
#
# Usage:
#   gaudi_jobs/calibration/check_selected_run_thresholds.sh              # print report
#   gaudi_jobs/calibration/check_selected_run_thresholds.sh --output FILE.txt
set -euo pipefail

# Run numbers for this selection: 60-64, 85-89, 142, 143, 146-151.
RUNS=(
  60 61 62 63 64
  85 86 87 88 89
  142 143
  146 147 148 149 150 151
)

# Must match calib_run_utils.DEFAULT_RAW_BASE. Some run numbers exist under
# the plain "TB2026CERN_run_NNNNNN" name and others under the eudaq-trigger
# "TB2026CERN_eudaq_run_NNNNNN" name (the prefix is just the trigger mode,
# not a different dataset) -- resolve whichever one actually exists on disk
# instead of assuming a fixed prefix.
RAW_BASE="/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata"

resolve_run_name() {
  local n="$1" plain eudaq
  plain=$(printf "TB2026CERN_run_%06d" "${n}")
  eudaq=$(printf "TB2026CERN_eudaq_run_%06d" "${n}")
  if [ -d "${RAW_BASE}/${plain}" ]; then
    echo "${plain}"
  elif [ -d "${RAW_BASE}/${eudaq}" ]; then
    echo "${eudaq}"
  else
    echo "ERROR: no run directory found for run number ${n} (tried ${plain} and ${eudaq})" >&2
    exit 1
  fi
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH:-}"

runs_csv=""
for n in "${RUNS[@]}"; do
  run_name=$(resolve_run_name "${n}")
  runs_csv="${runs_csv}${run_name},"
done
runs_csv="${runs_csv%,}"

exec python3 "${SCRIPT_DIR}/list_run_thresholds.py" --runs "${runs_csv}" "$@"
