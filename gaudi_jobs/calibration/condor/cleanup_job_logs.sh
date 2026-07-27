#!/bin/bash
# DAGMan SCRIPT POST helper. Deletes a successful node's captured stdout/
# stderr (.out/.err) so the AFS logs/ dir doesn't hit its ~2 MB directory-
# entry limit -- with thousands of nodes the per-job .out/.err otherwise fill
# the directory object and new long-named log files fail with "errno 27 File
# too large" (see gaudi_jobs/calibration/condor/README.md). Nodes that FAIL
# keep their .out/.err for debugging; the per-node .log (which DAGMan itself
# reads for job events) is never touched.
#
# CRITICAL: with a POST script DAGMan takes the NODE's success from the POST's
# exit code, not the job's -- so this must propagate the node's real return
# code (exit "$NODE_RC"). Exiting 0 unconditionally would mark a failed job as
# succeeded and skip its RETRY; exiting nonzero on success would fail a good
# node. Propagating $NODE_RC keeps RETRY semantics identical to having no POST.
#
# Usage: cleanup_job_logs.sh <node_return_code> <logbase>
#   where logbase is the .../logs/<JOB> stem (no extension); <logbase>.out and
#   <logbase>.err are removed on success -- AND SO ARE <logbase>_<N>.out/.err.
#
# The _<N> forms matter: a node that queues MANY processes (the CONVERT stage
# submits one per chunk group) names its files logs/<JOB>_$(Process).out, so the
# bare "${LOGBASE}.out" matched nothing and those nodes were never cleaned at
# all. That went unnoticed until a 297-node campaign left 1889 .out + 1889 .err
# behind on AFS (2026-07-25) -- exactly the directory-entry pressure this script
# exists to prevent. Single-process nodes (the LCIO path, the calibration Fits)
# keep matching the bare form, so both layouts are covered.
NODE_RC="$1"
LOGBASE="$2"

if [ "$NODE_RC" = "0" ]; then
  rm -f "${LOGBASE}.out" "${LOGBASE}.err"
  rm -f "${LOGBASE}"_[0-9]*.out "${LOGBASE}"_[0-9]*.err
fi

exit "$NODE_RC"
