#!/bin/bash
# trains_out_of_spill: the isolated above-cut acquisitions that sit BETWEEN
# spill trains, extracted and event-built so they can be opened in event_viewer.
#
# Two stages, both cheap (three acquisitions): a pre-scan over the run's chunks
# copying just the wanted entries, then the ordinary event builder over that.
# Runs in a few minutes on one core -- no Condor needed.
set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"

RUN="${RUN:-/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi/TB2026CERN_run_000060}"
ACQS="${ACQS:-2631,16038,26053}"
DATA="$HERE/data"
DECODED="$DATA/decoded_trains_out_of_spill.root"
ECAL="$DATA/ecal_trains_out_of_spill.root"

source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "$(cat "$REPO/.key4hep-release")"
export LD_LIBRARY_PATH="$REPO/gaudi_source/build:${LD_LIBRARY_PATH:-}"
export PYTHONPATH="$REPO/gaudi_source/build/genConfDir:$REPO:${PYTHONPATH:-}"
mkdir -p "$DATA"

echo "=== stage 1: extract acquisitions $ACQS from $(basename "$RUN") ==="
RUN="$RUN" ACQS="$ACQS" OUT="$DECODED" \
    nice -n 19 python3 "$HERE/extract_acquisitions.py"

echo
echo "=== stage 2: event build ==="
TOS_INPUT="$DECODED" TOS_OUTPUT="$ECAL" \
    nice -n 19 k4run "$HERE/steer_trains_out_of_spill.py"

echo
echo "=== done ==="
echo "open it with:"
echo "  python3 -m event_viewer --file $ECAL"
