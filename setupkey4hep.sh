# Bare key4hep environment, without the repo's PYTHONPATH/build wiring.
# For the full environment use `source setup.sh` instead.
#
# Release comes from .key4hep-release (the same file setup.sh, install.sh and
# the Condor generators read), so this cannot drift from the rest of the repo.
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r "$(cat "$(dirname "${BASH_SOURCE[0]:-$0}")/.key4hep-release" 2>/dev/null || echo 2026-04-08)"
