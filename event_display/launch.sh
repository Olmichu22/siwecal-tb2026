#!/usr/bin/env bash
# SiW-ECAL event display launcher
#
# Usage:
#   ./launch.sh <file.valcache.root> [start_event]
#
# Example:
#   ./launch.sh ../data/TB2026CERN_run_000013/ecal_TB2026CERN_run_000013.valcache.root 42
#
# The ROOT file must have hit_x, hit_y, hit_slab, hit_energy branches.
# start_event defaults to 0; relaunch with a different event index to browse.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_FILE="${1:?Usage: $0 <ecal_xyz.root> [start_event]}"
EVENT="${2:-0}"

python3 "${SCRIPT_DIR}/ecal_event_display.py" \
    --input "${ROOT_FILE}" \
    --event "${EVENT}"
