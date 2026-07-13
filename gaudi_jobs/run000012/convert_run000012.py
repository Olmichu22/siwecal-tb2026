#!/usr/bin/env python
"""
PIPELINE 1/2 -- decode TB2026CERN_run_000012 (th230) to siwecaldecoded chunks.

A thin wrapper over gaudi_jobs/run_convert_batch.py, which does the real work:
one k4run per raw chunk into <converted-dir>/<run>/chunks/, then a health check
that the acquisitions add up.

It used to be a k4run steering file that handed EVERY raw chunk of the run to a
single EcalRawDecoder. On this run that silently dropped 75% of the acquisitions
(30,020 decoded out of 119,211) and produced a perfectly valid-looking ROOT file;
the reconstruction built from it had sigma/mu = 1.06, a width larger than the
mean. Nothing failed, nothing warned. See gaudi_jobs/decode_chunks.py.

There is no per-run decoding logic left here on purpose: a copy of the pipeline
per run is how one copy gets fixed and the others quietly do not.

Usage::

    source setup.sh
    export LD_LIBRARY_PATH=$PWD/gaudi_source/build:$LD_LIBRARY_PATH
    export PYTHONPATH=$PWD/gaudi_source/build/genConfDir:$PWD:$PYTHONPATH
    python gaudi_jobs/run000012/convert_run000012.py      # NOT k4run
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from run_convert_batch import main  # noqa: E402

_RUN = "TB2026CERN_run_000012"

if __name__ == "__main__":
    sys.exit(main(["--runs", _RUN] + sys.argv[1:]))
