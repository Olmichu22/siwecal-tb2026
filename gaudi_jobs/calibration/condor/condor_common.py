"""
Shared helpers for the Condor pipeline generators (generate_convert_jobs.py,
generate_dag.py): repo/environment paths, EOS scratch layout, and small
HTCondor submit-file writing utilities.

Plain sibling module (not a package) -- both consumers live directly in
gaudi_jobs/calibration/condor/, so ``import condor_common`` resolves via
Python's automatic inclusion of a running script's own directory in
sys.path, same convention already used by calib_run_utils.py one directory
up.
"""
import os

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, "..", "..", ".."))
OPTIONS_DIR = os.path.join(REPO_ROOT, "gaudi_source", "options")

# Must match the release used elsewhere in this repo (gaudi_source/README.md
# "Build & run" instructions) -- Condor worker nodes source the same CVMFS
# release so the AFS-resident gaudi_source/build is ABI-compatible.
KEY4HEP_RELEASE = "2026-04-08"

# EOS scratch area for the Condor pipeline's intermediate files. Lives
# INSIDE rundata_converted_test/ (already proven writable all session)
# rather than as its own top-level sibling under .../Data/. NEVER
# auto-deleted by any script here (see the repo-wide "never delete
# anything under /Data/" rule).
DEFAULT_FILL_SCRATCH_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_test/calib_fill_scratch"

# IMPORTANT: shell `mkdir -p` reproducibly FAILS on this EOS FUSE mount --
# even creating a single new level under an already-existing, writable
# parent -- with a misleading "Permission denied" pointing at the mount
# root, not the real target. Plain single-level `mkdir` and Python's
# os.makedirs() (which internally does one os.mkdir() per path component,
# not a `mkdir -p` syscall) both work fine. Every generated wrapper script
# below MUST create directories via `mkdirs_line()`, never via `mkdir -p`.


def mkdirs_line(dir_expr):
    """Shell line that creates `dir_expr` (a shell expression, e.g. a
    variable or `$(dirname "$X")`) via Python's os.makedirs instead of
    `mkdir -p`, which is broken on the EOS FUSE mount used here (see above)."""
    return f'python3 -c "import os, sys; os.makedirs(sys.argv[1], exist_ok=True)" "{dir_expr}"\n'


def decoded_dir(fill_scratch_dir):
    """Where CONVERT (D.1) writes per-chunk decoded siwecaldecoded ROOT files."""
    return os.path.join(fill_scratch_dir, "decoded")


def hist_dir(fill_scratch_dir):
    """Where CALIBRACIÓN/Fill (D.2) writes per-chunk/per-run/per-threshold histogram ROOT files."""
    return os.path.join(fill_scratch_dir, "hist")


def env_wrapper_preamble(repo_root=REPO_ROOT, release=KEY4HEP_RELEASE):
    """Shell snippet every executable wrapper (.sh) sources to set up key4hep
    + LD_LIBRARY_PATH/PYTHONPATH pointing at the already-built gaudi_source/build.

    Every generated wrapper must use `set -eo pipefail` (NOT `-u`/nounset):
    key4hep's /cvmfs/sw.hsf.org/key4hep/setup.sh is sourced (not executed),
    so it runs in the caller's own shell, and it internally references an
    unset variable partway through its own argument handling -- harmless
    normally, but fatal under `set -u` ("unbound variable"). Verified this
    is specifically about `-u`, not about the wrapper's own positional
    parameters leaking into the sourced script.
    """
    return (
        f"source /cvmfs/sw.hsf.org/key4hep/setup.sh -r {release}\n"
        f"REPO={repo_root}\n"
        'export LD_LIBRARY_PATH="$REPO/gaudi_source/build:${LD_LIBRARY_PATH:-}"\n'
        'export PYTHONPATH="$REPO/gaudi_source/build/genConfDir:$REPO:${PYTHONPATH:-}"\n'
    )


def write_executable(path, content):
    with open(path, "w") as fh:
        fh.write(content)
    os.chmod(path, 0o755)


def write_text(path, content):
    with open(path, "w") as fh:
        fh.write(content)
