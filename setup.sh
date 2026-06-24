#!/usr/bin/env bash
# Environment setup for the SiW-ECAL TB2026CERN analysis suite.
#
# Source it (do not execute):   source setup.sh
#
# It loads the key4hep stack (numpy/scipy/pandas/sklearn/matplotlib/uproot/ROOT)
# from cvmfs and, if present, activates the local virtualenv that adds the
# event viewer's web stack (dash/plotly). See README.md for first-time setup.

# --- key4hep release (numpy, scipy, uproot, ROOT, ...) ----------------------
KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-2026-04-08}"
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${KEY4HEP_RELEASE}"

# --- repo root on PYTHONPATH so the top-level packages import anywhere -------
REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH}"

# --- optional viewer virtualenv (dash/plotly on top of key4hep) -------------
if [ -f "${REPO_ROOT}/.venv-viewer/bin/activate" ]; then
    source "${REPO_ROOT}/.venv-viewer/bin/activate"
    echo "[setup] key4hep ${KEY4HEP_RELEASE} + .venv-viewer active"
else
    echo "[setup] key4hep ${KEY4HEP_RELEASE} active (no .venv-viewer;"
    echo "        the event viewer needs dash/plotly -- see README.md)"
fi
