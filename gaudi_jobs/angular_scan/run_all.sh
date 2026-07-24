#!/bin/bash
# Run every steer_*.py + validate_*.py pair in this folder, sequentially.
#
# Idempotent: a run whose ecal_<RUN>.edm4hep.root already exists is skipped,
# unless --force. Assumes setup.sh + LD_LIBRARY_PATH/PYTHONPATH are already
# sourced (see README.md).
#
# Usage:
#   bash gaudi_jobs/angular_scan/run_all.sh [--force]
set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FORCE=0
[ "$1" = "--force" ] && FORCE=1

CONVERTED="${EVBLD_CONVERTED_DIR:-/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi}"

N=0
for steer in "$HERE"/steer_*.py; do
  RUN="$(basename "$steer" .py | sed 's/^steer_//')"
  VALIDATE="$HERE/validate_${RUN}.py"
  EDM4HEP="$CONVERTED/$RUN/events/ecal_${RUN}.edm4hep.root"

  if [ "$FORCE" != "1" ] && [ -s "$EDM4HEP" ]; then
    echo "[angular_scan] $RUN: already done -> $EDM4HEP (skip; use --force to redo)"
    continue
  fi

  echo "[angular_scan] $RUN: event builder"
  k4run "$steer"
  echo "[angular_scan] $RUN: validation (EDM4hep/PID)"
  k4run "$VALIDATE"
  N=$((N+1))
done

echo "[angular_scan] $N run(s) processed."
