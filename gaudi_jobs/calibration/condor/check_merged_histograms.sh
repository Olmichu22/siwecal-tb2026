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

# Release from .key4hep-release (repo root, three levels up from this script),
# the same file setup.sh and condor_common.py read.
_K4REL="$(cat "$(dirname "${BASH_SOURCE[0]:-$0}")/../../../.key4hep-release" 2>/dev/null || echo 2026-04-08)"
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "${_K4REL}" >/dev/null 2>&1

# The merge job's own `cp -f "$TMP_OUT" "$OUT"` (see generate_dag.py's
# _merge_sh) has already returned by the time this POST script starts, but
# it runs on the worker node while this check runs on the submit/AFS host --
# a different EOS FUSE client. Observed directly on run_000089 AND again on
# th220's own merge (twice, on 07/08 and 07/09): the check reports
# IsZombie/0-keys for a solid ~75-80s after job termination, while the
# exact same `root` open on the exact same file minutes-to-hours later (and
# every time since) succeeds cleanly with the full histogram count. That's
# EOS metadata propagation lag across hosts, not a corrupt merge. A first
# attempt at retrying (5x15s = 75s total) was NOT enough -- it failed
# consistently at exactly that ~80s mark, meaning the real lag exceeds it.
# Retry much more generously (up to 10 minutes) before giving up for real;
# `ls -la` on the parent dir before each open attempts to nudge the FUSE
# client into refreshing its cached directory listing.
NKEYS=-1
PARENT_DIR="$(dirname "${MERGED}")"
for attempt in $(seq 1 20); do
  ls -la "${PARENT_DIR}" >/dev/null 2>&1 || true
  NKEYS=$(root -l -b -q -e 'TFile f("'"${MERGED}"'"); if (f.IsZombie()) { std::cout << -1 << std::endl; } else { std::cout << f.GetListOfKeys()->GetSize() << std::endl; }' 2>/dev/null | tail -1)
  if [ -n "${NKEYS}" ] && [ "${NKEYS}" -gt 0 ] 2>/dev/null; then
    break
  fi
  echo "check_merged_histograms: attempt ${attempt}/20 could not open ${MERGED} (or 0 keys) -- retrying in 30s" >&2
  sleep 30
done

if [ -z "${NKEYS}" ] || [ "${NKEYS}" -le 0 ]; then
  echo "check_merged_histograms: ${MERGED} could not be opened or has zero histograms after 20 attempts (~10 min) -- refusing to run Fit" >&2
  exit 1
fi

echo "check_merged_histograms: ${MERGED} has ${NKEYS} histograms, looks sane"
exit 0
