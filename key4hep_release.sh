# Resolve the key4hep release this checkout is pinned to.
#
# Sourced by setup.sh, install.sh, setupkey4hep.sh and the batch wrappers, so
# there is ONE implementation of "read .key4hep-release and check it" instead of
# five copies that can disagree. The Python side (gaudi_jobs/calibration/condor/
# condor_common.py) applies the same rules for the Condor generators.
#
# Usage:
#     source "<repo_root>/key4hep_release.sh"
#     RELEASE="$(k4_release "<repo_root>")" || return 1   # or `exit 1`
#
# Prints the release on stdout, or explains the problem on stderr and returns 1.
# It deliberately does NOT fall back to a baked-in default: a silent fallback is
# how a checkout ends up building against one release and running against
# another, which surfaces as an ABI error on a worker node hours later rather
# than as a message here.

k4_release() {
    local repo_root="${1:?k4_release: repo root argument required}"
    local f="${repo_root}/.key4hep-release"

    if [ ! -f "$f" ]; then
        echo "ERROR: key4hep release file not found: $f" >&2
        echo "       It pins the CVMFS release that this repo builds and runs against," >&2
        echo "       and it is version-controlled -- if it is gone, the checkout is" >&2
        echo "       incomplete. Restore it with:  git checkout -- .key4hep-release" >&2
        echo "       (or write one:  echo 2026-04-08 > $f )" >&2
        echo "       To override for one shell instead:  export KEY4HEP_RELEASE=<release>" >&2
        return 1
    fi

    # First non-blank line, stripped of surrounding whitespace and any CR from a
    # file that took a trip through Windows.
    local value
    value="$(sed -e 's/\r$//' -e 's/^[[:space:]]*//' -e 's/[[:space:]]*$//' "$f" \
             | grep -v '^$' | head -n 1)"

    if [ -z "$value" ]; then
        echo "ERROR: key4hep release file is empty: $f" >&2
        echo "       Expected a single release name, e.g. 2026-04-08" >&2
        return 1
    fi

    # Conservative charset: this value is interpolated into a `source ... -r
    # <value>` command line, so anything outside [A-Za-z0-9._-] is rejected
    # rather than passed to a shell.
    case "$value" in
        *[!A-Za-z0-9._-]*)
            echo "ERROR: invalid key4hep release in $f: '$value'" >&2
            echo "       Only letters, digits, dot, underscore and hyphen are allowed" >&2
            echo "       (this value is passed to the CVMFS setup script)." >&2
            echo "       Expected something like: 2026-04-08" >&2
            return 1
            ;;
    esac

    # "latest" is a moving target, not a directory to check -- see
    # k4_setup_release for why it is NOT the releases/latest directory.
    if [ "$value" = "latest" ]; then
        printf '%s\n' "$value"
        return 0
    fi

    if [ ! -d "/cvmfs/sw.hsf.org/key4hep/releases/$value" ]; then
        # A warning, not an error: /cvmfs may simply not be mounted yet on this
        # host, and nightlies do not live under releases/. The release name
        # itself is well-formed, so let the caller proceed and let the CVMFS
        # setup script have the final word.
        echo "WARNING: key4hep release '$value' (from $f) is not under" >&2
        echo "         /cvmfs/sw.hsf.org/key4hep/releases/ -- continuing, but if the" >&2
        echo "         setup below fails, check the release name or that CVMFS is mounted." >&2
    fi

    printf '%s\n' "$value"
}

# Translate a release name into the argument to give to `key4hep/setup.sh -r`.
#
# Only "latest" needs translating, and it is a genuine trap:
#
#     /cvmfs/sw.hsf.org/key4hep/releases/latest      -> 2024-04-12  (dead since)
#     /cvmfs/sw.hsf.org/key4hep/releases/latest-opt  -> the real latest
#
# `-r latest` therefore silently gives you a 2024 stack with ROOT 6.28 instead
# of the current one. key4hep's own default (no -r) is `latest-$build_type`,
# i.e. latest-opt, which is what we pass.
k4_setup_release() {
    case "${1:?k4_setup_release: release argument required}" in
        latest) printf 'latest-opt\n' ;;
        *)      printf '%s\n' "$1" ;;
    esac
}

# Is the key4hep already loaded in THIS shell the release we want?
#
# key4hep refuses to be sourced twice in one process, so setup.sh has to tell
# "already loaded, nothing to do" from "loaded, but a different release". The
# loaded stack is identified by KEY4HEP_STACK:
#
#     /cvmfs/.../releases/<release>/<platform>/key4hep-stack/<...>/setup.sh
#
# For "latest" the comparison needs the concrete release behind the moving
# name, which is exactly what the latest-opt symlink for the *loaded* platform
# resolves to -- so the platform is taken from KEY4HEP_STACK rather than
# re-deriving the OS/compiler that key4hep already decided on.
#
# Returns 0 if they match, 1 otherwise (including "nothing loaded").
k4_stack_is() {
    local wanted="${1:?k4_stack_is: release argument required}"
    local stack="${KEY4HEP_STACK:-}"
    [ -n "$stack" ] || return 1

    local rest="${stack#*/key4hep/releases/}"
    [ "$rest" != "$stack" ] || return 1
    local loaded_rel="${rest%%/*}"
    rest="${rest#*/}"
    local platform="${rest%%/*}"

    if [ "$wanted" = "latest" ]; then
        local target
        target="$(readlink -f "/cvmfs/sw.hsf.org/key4hep/releases/latest-opt/${platform}" 2>/dev/null)" \
            || return 1
        target="${target%/*}"          # strip <platform>
        wanted="${target##*/}"         # concrete release behind latest-opt
        [ -n "$wanted" ] || return 1
    fi

    [ "$loaded_rel" = "$wanted" ]
}
