#!/bin/bash
set -eo pipefail
# steer_validate.sh <run_dir> <run_tag>
#
# Condor executable for the per-run standalone pipelines (gaudi_jobs/run0000NN/):
# runs the k4run event-build steer, then the k4run EDM4hep/PID validate, in one
# job. Same key4hep + LD_LIBRARY_PATH/PYTHONPATH setup as reco.sh, so it behaves
# identically to running the steer/validate by hand -- just on a worker node.
#
#   <run_dir>  absolute path to the run's pipeline folder (has steer_<tag>.py /
#              validate_<tag>.py), e.g. .../gaudi_jobs/run000166
#   <run_tag>  the run tag, e.g. run000166
RUN_DIR="$1"; RUN_TAG="$2"

source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
REPO=/afs/cern.ch/user/m/marquezh/public/siwecal-tb2026
export LD_LIBRARY_PATH="$REPO/gaudi_source/build:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$REPO/gaudi_source/build/genConfDir:$REPO:${PYTHONPATH:-}"

echo "[pipeline] $RUN_TAG: event build (steer_$RUN_TAG.py)"
k4run "$RUN_DIR/steer_$RUN_TAG.py"

echo "[pipeline] $RUN_TAG: validate / EDM4hep+PID (validate_$RUN_TAG.py)"
k4run "$RUN_DIR/validate_$RUN_TAG.py"

echo "[pipeline] $RUN_TAG: done"
