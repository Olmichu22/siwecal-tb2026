# Bare key4hep environment, without the repo's PYTHONPATH/build wiring.
# For the full environment use `source setup.sh` instead.
#
# Release comes from .key4hep-release via k4_release (the same helper setup.sh,
# install.sh and the batch wrappers use), so this cannot drift from the rest of
# the repo, and a missing or malformed file says so instead of guessing.
_K4_REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
source "${_K4_REPO_ROOT}/key4hep_release.sh"
# As in setup.sh: sourcing key4hep twice in one shell is refused by key4hep
# itself, so treat "the release we want is already loaded" as success. And as
# there, the release may be the moving name "latest" -- k4_setup_release /
# k4_stack_is handle it (they know latest is latest-opt, not releases/latest).
_K4_RELEASE="${KEY4HEP_RELEASE:-$(k4_release "${_K4_REPO_ROOT}")}" \
    && { k4_stack_is "${_K4_RELEASE}" \
         || source /cvmfs/sw.hsf.org/key4hep/setup.sh \
                   -r "$(k4_setup_release "${_K4_RELEASE}")" ; }
unset _K4_REPO_ROOT _K4_RELEASE
