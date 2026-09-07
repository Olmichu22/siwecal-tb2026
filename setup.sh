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
if [ -z "${KEY4HEP_RELEASE:-}" ]; then
    source "${REPO_ROOT}/key4hep_release.sh"
    KEY4HEP_RELEASE="$(k4_release "${REPO_ROOT}")" || return 1
fi
# key4hep refuses to be sourced twice in one shell ("already set up, please
# start a new shell"), and returns non-zero when it does. That made re-sourcing
# setup.sh -- and launching a Jupyter kernel (jupyter_kernel.sh) from a shell
# that already had the stack -- fail for no good reason. So: if the release we
# want is already loaded, skip the cvmfs step and just (re)apply the repo
# wiring below; if a *different* release is loaded, say so instead of mixing
# two releases in one process.
#
# KEY4HEP_RELEASE may be the moving name "latest"; k4_setup_release turns it
# into what `-r` actually understands and k4_stack_is resolves it for the
# comparison (see key4hep_release.sh).
if ! command -v k4_stack_is >/dev/null 2>&1; then
    source "${REPO_ROOT}/key4hep_release.sh"
fi
if k4_stack_is "${KEY4HEP_RELEASE}"; then
    :   # same release already loaded, nothing to do
elif [ -n "${KEY4HEP_STACK:-}" ]; then
    echo "ERROR: this shell already has a different key4hep loaded:" >&2
    echo "         ${KEY4HEP_STACK}" >&2
    echo "       Requested: ${KEY4HEP_RELEASE}. Start a new shell." >&2
    return 1
else
    source /cvmfs/sw.hsf.org/key4hep/setup.sh \
           -r "$(k4_setup_release "${KEY4HEP_RELEASE}")" || return 1
fi

# --- repo on PYTHONPATH so the top-level packages import anywhere ------------
export PYTHONPATH="${REPO_ROOT}:${PYTHONPATH}"

# --- optional k4SiWEcalReco build (so k4run finds the compiled plugin) -------
# Built once by install.sh (or by hand, see gaudi_source/README.md). Without
# these `k4run` cannot load libk4SiWEcalRecoPlugins.so. The source tree is
# gaudi_source/ (the batch jobs export exactly this path).
#
# There is deliberately NO fallback to the pre-rename k4SiWEcalReco/build. That
# directory survives the rename as an untracked leftover, so a fallback fires
# precisely when gaudi_source/build is missing -- silently putting a months-old
# .so on LD_LIBRARY_PATH at the one moment you need to notice the build is gone.
if [ -d "${REPO_ROOT}/gaudi_source/build" ]; then
    export LD_LIBRARY_PATH="${REPO_ROOT}/gaudi_source/build:${LD_LIBRARY_PATH}"
    export PYTHONPATH="${REPO_ROOT}/gaudi_source/build/genConfDir:${PYTHONPATH}"
fi

# --- what to call the release in the messages below -------------------------
# With KEY4HEP_RELEASE=latest the interesting number is the one latest resolved
# to, which key4hep exports as key4hep_stack_version -- printing "latest" alone
# would hide the fact that it moved under you.
if [ "${KEY4HEP_RELEASE}" = "latest" ] && [ -n "${key4hep_stack_version:-}" ]; then
    _K4_SHOWN="latest (${key4hep_stack_version})"
else
    _K4_SHOWN="${KEY4HEP_RELEASE}"
fi

# --- optional viewer virtualenv (dash/plotly on top of key4hep) -------------
if [ -f "${REPO_ROOT}/.venv-viewer/bin/activate" ]; then
    # Skip if it is already the active venv, so re-sourcing setup.sh does not
    # stack another copy of it on PATH.
    if [ "${VIRTUAL_ENV:-}" != "${REPO_ROOT}/.venv-viewer" ]; then
        source "${REPO_ROOT}/.venv-viewer/bin/activate"
    fi
    echo "[setup] key4hep ${_K4_SHOWN} + .venv-viewer active"
else
    echo "[setup] key4hep ${_K4_SHOWN} active (no .venv-viewer;"
    echo "        the event viewer needs dash/plotly -- see README.md)"
fi
unset _K4_SHOWN
