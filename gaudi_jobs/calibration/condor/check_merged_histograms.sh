#!/bin/bash
# DAGMan SCRIPT POST for the merge_threshold node: refuses to let the Fit
# jobs run against an empty/corrupt merge. DAGMan passes the node's own
# exit status as $1 (POST scripts run regardless of the node's success, so
# this is checked explicitly) followed by the VARS-supplied arguments.
#
# Usage: check_merged_histograms.sh <node_return_code> <merged_root_file>
# NOTE: no `-u` (nounset) -- key4hep's setup.sh isn't nounset-safe and
# blows up with "unbound variable" if sourced under `set -u`.
set -eo pipefail
NODE_RC="$1"
MERGED="$2"

if [ "${NODE_RC}" != "0" ]; then
  echo "check_merged_histograms: merge_threshold node itself failed (rc=${NODE_RC})" >&2
  exit 1
fi

source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08 >/dev/null 2>&1

NKEYS=$(root -l -b -q -e 'TFile f("'"${MERGED}"'"); if (f.IsZombie()) { std::cout << -1 << std::endl; } else { std::cout << f.GetListOfKeys()->GetSize() << std::endl; }' 2>/dev/null | tail -1)

if [ -z "${NKEYS}" ] || [ "${NKEYS}" -le 0 ]; then
  echo "check_merged_histograms: ${MERGED} could not be opened or has zero histograms -- refusing to run Fit" >&2
  exit 1
fi

echo "check_merged_histograms: ${MERGED} has ${NKEYS} histograms, looks sane"
exit 0
