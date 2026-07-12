#!/usr/bin/env bash
# Environment setup for the SiW-ECAL TB2026CERN analysis suite.
#
# Source it (do not execute):   source setup.sh
#
# It loads the key4hep stack (numpy/scipy/pandas/sklearn/matplotlib/uproot/ROOT)
# from cvmfs and, if present, activates the local virtualenv that adds the
# event viewer's web stack (dash/plotly). See README.md for first-time setup.

# --- repo root (also needed to read the shared key4hep release file) --------
REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"

# --- key4hep release (numpy, scipy, uproot, ROOT, ...) ----------------------
# Single source of truth: .key4hep-release (shared with install.sh); override
# per-shell with the KEY4HEP_RELEASE env var.
KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-$(cat "${REPO_ROOT}/.key4hep-release" 2>/dev/null || echo 2026-04-08)}"
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${KEY4HEP_RELEASE}"

# --- repo on PYTHONPATH so the top-level packages import anywhere ------------
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH}"

# --- optional k4SiWEcalReco build (so k4run finds the compiled plugin) -------
# Built once by install.sh (or by hand, see k4SiWEcalReco/README.md). Without
# these run_pid_batch.py's `k4run` cannot load libk4SiWEcalRecoPlugins.so.
if [ -d "${REPO_ROOT}/k4SiWEcalReco/build" ]; then
    export LD_LIBRARY_PATH="${REPO_ROOT}/k4SiWEcalReco/build:${LD_LIBRARY_PATH}"
    export PYTHONPATH="${REPO_ROOT}/k4SiWEcalReco/build/genConfDir:${PYTHONPATH}"
fi

# --- optional viewer virtualenv (dash/plotly on top of key4hep) -------------
if [ -f "${REPO_ROOT}/.venv-viewer/bin/activate" ]; then
    source "${REPO_ROOT}/.venv-viewer/bin/activate"
    echo "[setup] key4hep ${KEY4HEP_RELEASE} + .venv-viewer active"
else
    echo "[setup] key4hep ${KEY4HEP_RELEASE} active (no .venv-viewer;"
    echo "        the event viewer needs dash/plotly -- see README.md)"
fi
