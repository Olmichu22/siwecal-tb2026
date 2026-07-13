#!/usr/bin/env python
"""
PIPELINE 2/2 -- reconstruct + validate TB2026CERN_run_000012 (th230).

A thin wrapper over gaudi_jobs/run_full_pipeline_batch.py (decode -> event
building -> PID), plus the validation pass. The threshold, the calibration
tables and the geometry are all resolved by that driver from the run's own
Run_Settings.txt -- they are not restated here, so this run cannot drift out of
step with every other one.

The decode stage is idempotent: chunks already decoded (by pipeline 1, or by the
Condor DAG) are reused, so re-running this to try different calibration files
costs nothing extra.

Outputs land in Reconstruction/<run>/ -- the events and the plots made from them,
together, outside Data/.

Usage::

    source setup.sh
    export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
    export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
    python gaudi_jobs/run000012/process_run000012.py
"""
import os
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from run_full_pipeline_batch import main as reconstruct  # noqa: E402

from siwecal_common import paths  # noqa: E402

_RUN = "TB2026CERN_run_000012"

if __name__ == "__main__":
    rc = reconstruct(["--run", _RUN] + sys.argv[1:])
    if rc:
        sys.exit(rc)

    ecal = os.path.join(paths.reconstruction_run_dir(_RUN), f"ecal_{_RUN}.root")
    print(f"[validation] {_RUN}: plots from {ecal}")
    result = subprocess.run([sys.executable, "-m", "siwecal_validation",
                             "--file", ecal, "--run", _RUN])
    sys.exit(result.returncode)
