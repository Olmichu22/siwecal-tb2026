#!/usr/bin/env python
"""
Full pipeline for TB2026CERN_run_000013 (52 GeV electrons, beam point P1_52GeV
-- see configs/data/data_reference_base.yml): decode -> event building ->
PID/EDM4hep.

A thin wrapper over gaudi_jobs/run_full_pipeline_batch.py, which does the real
work. The threshold, the calibration tables and the geometry are all resolved
there from the run's own Run_Settings.txt.

This script used to carry its own copy of the pipeline, whose first stage handed
all 304 raw chunks to a single k4run process (``run_raw2root_and_eventbuilder.py``,
``TopAlg=[EcalRawDecoder, EcalEventBuilder]``). That pattern silently drops
acquisitions on some runs. run_000013 happens to survive it intact -- which is
precisely why it went unnoticed for so long: spot-checking this run "confirms"
the code is fine while run_000012 was losing 75% of its data. See
gaudi_jobs/decode_chunks.py.

Usage::

    source setup.sh
    export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
    export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
    python gaudi_jobs/run000013/run_full_pipeline_000013.py
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from run_full_pipeline_batch import main  # noqa: E402

_RUN = "TB2026CERN_run_000013"

if __name__ == "__main__":
    sys.exit(main(["--run", _RUN] + sys.argv[1:]))
