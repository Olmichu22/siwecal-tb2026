#!/usr/bin/env bash
# One-time installer for the SiW-ECAL TB2026CERN analysis suite.
#
# Run it (do NOT source):   ./install.sh [options]
# Afterwards, in each shell:   source setup.sh
#
# It sets up the three pieces of the stack, each skippable:
#   1) key4hep        sourced from cvmfs (provides numpy/scipy/uproot/ROOT/...)
#   2) .venv-viewer   the event viewer's web stack (dash/plotly) on top of it
#   3) k4SiWEcalReco  the compiled Gaudi plugin that `k4run` loads
#
# This is the heavy, occasional step (venv creation + C++ build). The per-shell
# environment is set up separately by `source setup.sh` (which install.sh does
# NOT replace).

set -eo pipefail

REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
# Single source of truth: .key4hep-release (shared with setup.sh); override
# with the KEY4HEP_RELEASE env var.
KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-$(cat "${REPO_ROOT}/.key4hep-release" 2>/dev/null || echo 2026-04-08)}"

DO_VIEWER=1
DO_K4=1
FORCE_VENV=0
JOBS="$(nproc 2>/dev/null || echo 4)"

usage() {
    cat <<EOF
Usage: ./install.sh [options]

  --no-viewer      skip the .venv-viewer (dash/plotly) step
  --no-k4          skip the k4SiWEcalReco C++ build
  --force-venv     delete and recreate .venv-viewer from scratch
  -j, --jobs N     parallel build jobs (default: ${JOBS})
  -h, --help       show this help

Environment:
  KEY4HEP_RELEASE  key4hep release to source (default: ${KEY4HEP_RELEASE})
EOF
}

while [ $# -gt 0 ]; do
    case "$1" in
        --no-viewer)  DO_VIEWER=0 ;;
        --no-k4)      DO_K4=0 ;;
        --force-venv) FORCE_VENV=1 ;;
        -j|--jobs)    JOBS="$2"; shift ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 2 ;;
    esac
    shift
done

# --- 1/3: key4hep -----------------------------------------------------------
echo "==> [1/3] key4hep ${KEY4HEP_RELEASE} (from cvmfs)"
if [ ! -r /cvmfs/sw.hsf.org/key4hep/setup.sh ]; then
    echo "ERROR: /cvmfs/sw.hsf.org/key4hep/setup.sh not readable -- is cvmfs mounted?" >&2
    exit 1
fi
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${KEY4HEP_RELEASE}"
echo "    python: $(command -v python)"

# --- 2/3: viewer virtualenv -------------------------------------------------
VENV="${REPO_ROOT}/.venv-viewer"
if [ "${DO_VIEWER}" = 1 ]; then
    echo "==> [2/3] event-viewer virtualenv (.venv-viewer)"
    if [ "${FORCE_VENV}" = 1 ] && [ -d "${VENV}" ]; then
        echo "    --force-venv: removing existing ${VENV}"
        rm -rf "${VENV}"
    fi
    if [ ! -d "${VENV}" ]; then
        python -m venv --system-site-packages "${VENV}"
    else
        echo "    reusing existing ${VENV} (pass --force-venv to rebuild)"
    fi
    # shellcheck disable=SC1091
    source "${VENV}/bin/activate"
    pip install --no-input -r "${REPO_ROOT}/requirements.txt"
    deactivate
else
    echo "==> [2/3] skipped (--no-viewer)"
fi

# --- 3/3: k4SiWEcalReco build -----------------------------------------------
if [ "${DO_K4}" = 1 ]; then
    echo "==> [3/3] building k4SiWEcalReco (-j${JOBS})"
    cmake -S "${REPO_ROOT}/k4SiWEcalReco" -B "${REPO_ROOT}/k4SiWEcalReco/build" \
          -DCMAKE_BUILD_TYPE=Release
    cmake --build "${REPO_ROOT}/k4SiWEcalReco/build" -j"${JOBS}"
else
    echo "==> [3/3] skipped (--no-k4)"
fi

cat <<EOF

Done. The stack is installed. In every new shell, load the environment with:

    source setup.sh

That sources key4hep, puts the repo + k4SiWEcalReco build on PYTHONPATH, and
activates .venv-viewer if present.
EOF
