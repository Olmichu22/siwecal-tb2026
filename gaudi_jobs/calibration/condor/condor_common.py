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
import re
import sys

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(_THIS_DIR, "..", "..", ".."))
OPTIONS_DIR = os.path.join(REPO_ROOT, "gaudi_source", "options")

_RELEASE_RE = re.compile(r"^[A-Za-z0-9._-]+$")


def _key4hep_release():
    """The key4hep release the generated Condor wrappers source on the worker.

    Read from ``.key4hep-release`` at the repo root -- the SAME file setup.sh,
    install.sh and key4hep_release.sh read, so the farm and an interactive shell
    cannot end up on different releases. This used to be a literal here, which
    made the "single source of truth" the README advertises quietly untrue:
    bumping the file and rebuilding gaudi_source/build against the new release
    left every Condor job still sourcing the old one, loading an AFS-resident
    .so with a mismatched ABI. That does not fail at submit time -- it fails on
    the worker, obscurely.

    Raises rather than falling back to a baked-in default: generating a whole
    DAG against a release nobody chose is exactly the silent drift this is meant
    to prevent, and a generator is a good place to fail -- nothing has been
    submitted yet. Mirrors k4_release() in key4hep_release.sh, including the
    conservative charset (the value is interpolated into the `source ... -r`
    line of every generated wrapper).
    """
    path = os.path.join(REPO_ROOT, ".key4hep-release")
    try:
        with open(path) as fh:
            lines = [ln.strip() for ln in fh.read().splitlines()]
    except FileNotFoundError:
        raise SystemExit(
            f"ERROR: key4hep release file not found: {path}\n"
            "       It pins the CVMFS release the generated Condor jobs source, and it is\n"
            "       version-controlled -- if it is gone, the checkout is incomplete.\n"
            "       Restore it with:  git checkout -- .key4hep-release")
    except OSError as e:
        raise SystemExit(f"ERROR: cannot read key4hep release file {path}: {e}")

    values = [ln for ln in lines if ln]
    if not values:
        raise SystemExit(
            f"ERROR: key4hep release file is empty: {path}\n"
            "       Expected a single release name, e.g. 2026-04-08")

    release = values[0]
    if not _RELEASE_RE.match(release):
        raise SystemExit(
            f"ERROR: invalid key4hep release in {path}: {release!r}\n"
            "       Only letters, digits, dot, underscore and hyphen are allowed (the value\n"
            "       is interpolated into the `source ... -r <release>` line of every\n"
            "       generated Condor wrapper).\n"
            "       Expected something like: 2026-04-08")

    cvmfs = f"/cvmfs/sw.hsf.org/key4hep/releases/{release}"
    if not os.path.isdir(cvmfs):
        # Warn, do not fail: CVMFS may not be mounted on the submit host, and
        # nightlies do not live under releases/. The name is well-formed; let
        # the worker's setup script have the final word.
        print(f"WARNING: key4hep release {release!r} (from {path}) is not under "
              f"/cvmfs/sw.hsf.org/key4hep/releases/ -- continuing.", file=sys.stderr)
    return release


KEY4HEP_RELEASE = _key4hep_release()

# EOS scratch area for the Fill stage's histograms. NEVER auto-deleted by any
# script here (see the repo-wide "never delete anything under /Data/" rule).
#
# Sits at the top level of Data/, NOT under a rundata_converted_* area: these
# are calibration histograms, not converted events, and the previous home
# (rundata_converted_test/) was deleted wholesale in July 2026 precisely
# because it looked like a scratch area for converted events. That took the
# merged_th220/th230.root the current MIP/pedestal tables were fitted from with
# it -- the tables themselves survived only because they live in git, under
# calibration/MuonCalib_gaudi/. Re-fitting now means re-running Fill from the
# chunks in rundata_converted_gaudi/<RUN>/chunks/.
DEFAULT_FILL_SCRATCH_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/calib_fill_scratch"

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


# Where a run's decoded data lives: <converted_dir>/<RUN>/chunks/.
#
# There is deliberately no merged per-run <RUN>.root any more. Every run used to
# get one, to match the collaboration's one-file-per-run layout, but nothing reads
# it: the calibration Fill reads the chunks and so does the event builder
# (EcalEventBuilder chains them via InputFiles). It was ~223 GB of duplicate data,
# and for run_000004 (~176 GB of chunks) it cannot even exist -- that crosses
# ROOT's 100 GB TTree::fgMaxTreeSize, which is what killed the hadd.
DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"


def chunks_dir(converted_dir, run):
    """Where CONVERT writes a run's decoded siwecaldecoded chunks:
    <converted_dir>/<RUN>/chunks/.

    ONE decoded area, shared by every consumer: the calibration Fill reads these
    chunks, and so does the event builder (EcalEventBuilder's InputFiles chains
    them). They used to live off in the calibration Condor scratch, which made
    the decoded data look like a calibration by-product -- it is not, it is the
    run's converted data, and it belongs next to the run's reconstructed output.

    There is deliberately no merged per-run siwecaldecoded file. Nothing needs
    one: it would duplicate ~390 GB across the calibration runs for something no
    analysis opens, and for run_000004 (~176 GB of chunks) it cannot even exist,
    since that crosses ROOT's 100 GB TTree::fgMaxTreeSize.
    """
    return os.path.join(converted_dir, run, "chunks")


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
