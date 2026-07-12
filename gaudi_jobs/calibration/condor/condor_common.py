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

# EOS scratch area for the Condor pipeline's intermediate files. NEVER
# auto-deleted by any script here (see the repo-wide "never delete anything
# under /Data/" rule).
#
# Deliberately still under rundata_converted_test/ even though converted
# events now go to rundata_converted_gaudi/: this holds ~1.1 TB of already
# validated per-chunk/per-run/per-threshold histograms (incl. the merged
# merged_th220/th230.root the current calibration tables were fitted from).
# It is calibration scratch, not converted events, and repointing it would
# orphan all of it and force a full re-fill for no gain.
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


# Where a fully converted run is PUBLISHED, one directory per run holding one
# ROOT file named after it: <converted_dir>/<RUN>/<RUN>.root. This is the
# collaboration's own layout (cf. rundata_converted_old/, and the runs written
# by gaudi_jobs/run000012, run000013 and run_calibration_batch.py) -- calibration
# runs follow it too, they are not a special case.
DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"


def converted_run_file(converted_dir, run):
    """The published converted file for `run`: <converted_dir>/<RUN>/<RUN>.root."""
    return os.path.join(converted_dir, run, f"{run}.root")


def decoded_dir(fill_scratch_dir):
    """Where CONVERT (D.1) writes per-chunk decoded siwecaldecoded ROOT files.

    These per-chunk files are an artefact of parallelising the decode across
    Condor jobs, NOT the converted run itself: the canonical converted run is
    the single <RUN>.root that generate_publish_converted_jobs.py hadds out of
    them into DEFAULT_CONVERTED_DIR.
    """
    return os.path.join(fill_scratch_dir, "decoded")


def hist_dir(fill_scratch_dir):
    """Where CALIBRACIÓN/Fill (D.2) writes per-chunk/per-run/per-threshold histogram ROOT files."""
    return os.path.join(fill_scratch_dir, "hist")


def condor_log_dir(dag_dir):
    """Where the Condor per-job log/output/error files go: a ``logs/`` subdir of
    the (AFS) submit directory.

    They MUST stay on AFS, not EOS: CERN's standard batch schedds reject /eos
    paths for log/output/error/executable outright ("Standard batch schedds
    cannot use /eos paths directly within the submit file"). The defence against
    thousands of logs filling the small AFS home is therefore keeping each log
    tiny -- Fill/Fit run Gaudi at WARNING and Fill skips already-done chunks, so
    their .out is a line or two; merges run ``hadd -v 0`` instead of listing
    every source file (which for a 2000-chunk run was thousands of lines)."""
    return os.path.join(dag_dir, "logs")


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
