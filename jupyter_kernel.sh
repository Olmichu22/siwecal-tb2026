#!/usr/bin/env bash
# Jupyter kernel for the SiW-ECAL TB2026CERN suite, running on the key4hep stack.
#
# Jupyter cannot `source` anything: a kernel is a fixed command line, so a
# notebook started from a plain `python3` kernel has no key4hep (no ROOT, no
# uproot) and no repo on PYTHONPATH. The fix is a kernelspec whose command is
# *this same script* in launch mode -- it sources the environment exactly like a
# terminal would and then execs the ipykernel launcher. One file, so the kernel
# can never drift from `setup.sh`.
#
#   ./jupyter_kernel.sh install      register the kernel (once per machine/user)
#   ./jupyter_kernel.sh lab          start a JupyterLab server in this env
#   ./jupyter_kernel.sh list         show installed kernels
#   ./jupyter_kernel.sh remove       unregister the kernel
#
# Run it, do NOT source it.

set -eo pipefail

REPO_ROOT="$( cd "$( dirname "${BASH_SOURCE[0]:-$0}" )" && pwd )"
SELF="${REPO_ROOT}/jupyter_kernel.sh"

# Kernelspecs live in the Jupyter data dir; the --user location is the one every
# jupyter on this account picks up (local JupyterLab, notebook, JupyterHub with
# the default paths), which is why the kernel is usable from a Jupyter that
# knows nothing about this repo.
JUPYTER_KERNELS_DIR="${JUPYTER_DATA_DIR:-${HOME}/.local/share/jupyter}/kernels"

usage() {
    cat <<EOF
Usage: ./jupyter_kernel.sh <command> [options]

Commands:
  install            register the Jupyter kernel for this checkout
  remove             delete the registered kernel
  list               list the kernels Jupyter can see
  lab [args...]      start JupyterLab inside the environment (extra args are
                     passed straight to \`jupyter lab\`)

install options:
  --bare             kernel loads only key4hep (setupkey4hep.sh) instead of the
                     full repo environment (setup.sh: PYTHONPATH + gaudi build
                     + .venv-viewer)
  --release REL      key4hep release the kernel loads (default: latest -- the
                     kernel then follows latest release by release, resolving it
                     at every start)
  --pinned           use the repo's pinned release (.key4hep-release) instead of
                     latest, so the kernel matches what gaudi_source was built
                     against
  --name NAME        kernelspec directory name
                     (default: siwecal-key4hep, or key4hep-<release> with --bare)
  --display-name S   name shown in the Jupyter launcher
  --force            overwrite an existing kernel with the same name
  -h, --help         this help

Environment:
  KEY4HEP_RELEASE    same as --release (--release wins if both are given)
EOF
}

# --- release ----------------------------------------------------------------
# The kernel deliberately does NOT default to .key4hep-release: that file pins
# what gaudi_source is *compiled* against, and for interactive notebooks the
# useful default is the current stack. --pinned opts back into the repo pin
# (needed if a notebook has to load the compiled k4SiWEcalReco plugin, whose
# ABI follows the release it was built with).
resolve_pinned_release() {
    # shellcheck disable=SC1091
    source "${REPO_ROOT}/key4hep_release.sh"
    k4_release "${REPO_ROOT}"
}

# --- environment loading, shared by `launch` and `lab` ----------------------
# Everything the setup scripts print goes to stderr: on the kernel's stdout it
# would sit in front of the ZMQ handshake, and a chatty stdout is exactly the
# kind of thing that turns into an unexplained "kernel died" in the UI.
load_env() {
    local mode="$1"   # full | bare
    if [ "${mode}" = "bare" ]; then
        # shellcheck disable=SC1091
        source "${REPO_ROOT}/setupkey4hep.sh" 1>&2
    else
        # shellcheck disable=SC1091
        source "${REPO_ROOT}/setup.sh" 1>&2
    fi
}

# --- commands ---------------------------------------------------------------
cmd_install() {
    local mode="full" name="" display="" force=0
    local release="${KEY4HEP_RELEASE:-latest}"

    while [ $# -gt 0 ]; do
        case "$1" in
            --bare)         mode="bare" ;;
            --release)      release="$2"; shift ;;
            --pinned)       release="pinned" ;;
            --name)         name="$2"; shift ;;
            --display-name) display="$2"; shift ;;
            --force)        force=1 ;;
            -h|--help)      usage; exit 0 ;;
            *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
        esac
        shift
    done

    if [ "${release}" = "pinned" ]; then
        release="$(resolve_pinned_release)" || exit 1
    fi

    # Release names end up in a directory name and on a `-r` command line.
    case "${release}" in
        *[!A-Za-z0-9._-]*|"")
            echo "ERROR: invalid release '${release}'" >&2
            exit 1 ;;
    esac

    if [ "${mode}" = "bare" ]; then
        name="${name:-key4hep-${release}}"
        display="${display:-key4hep ${release}}"
    else
        name="${name:-siwecal-key4hep}"
        display="${display:-SiW-ECAL (key4hep ${release})}"
    fi

    # Kernelspec names are used as directory names and in URLs; keep them tame.
    case "${name}" in
        *[!A-Za-z0-9._-]*|"")
            echo "ERROR: invalid kernel name '${name}'" >&2
            echo "       Only letters, digits, dot, underscore and hyphen." >&2
            exit 1 ;;
    esac

    local kernel_dir="${JUPYTER_KERNELS_DIR}/${name}"
    if [ -d "${kernel_dir}" ] && [ "${force}" != 1 ]; then
        echo "ERROR: kernel already installed: ${kernel_dir}" >&2
        echo "       Re-run with --force to overwrite, or pick another --name." >&2
        exit 1
    fi

    # Check that the environment can actually host a kernel *before* writing the
    # kernelspec, so a broken setup shows up here and not as a dead kernel in a
    # notebook. Runs in a subshell: this process must not inherit key4hep.
    echo "==> checking key4hep ${release} (${mode} environment)"
    local resolved
    resolved="$(
        set +e
        export KEY4HEP_RELEASE="${release}"
        load_env "${mode}"
        python -c 'import ipykernel' 2>/dev/null || exit 1
        # What "latest" meant today -- reported below, not baked in.
        printf '%s\n' "${key4hep_stack_version:-${release}}"
    )" || {
        echo "ERROR: the ${mode} environment has no importable ipykernel." >&2
        echo "       Is cvmfs mounted and is '${release}' a valid release?" >&2
        exit 1
    }

    mkdir -p "${kernel_dir}"
    # The release is baked in as the *name* the user asked for, not as what it
    # resolves to today: with "latest" the kernel re-resolves it at every start
    # and so follows the stack, while an explicit date stays put forever.
    cat > "${kernel_dir}/kernel.json" <<EOF
{
 "argv": [
  "${SELF}",
  "launch",
  "--mode", "${mode}",
  "--release", "${release}",
  "-f", "{connection_file}"
 ],
 "display_name": "${display}",
 "language": "python",
 "metadata": {
  "debugger": true
 }
}
EOF

    echo "==> installed kernel '${name}' -> ${kernel_dir}"
    echo "    display name : ${display}"
    if [ "${release}" = "latest" ]; then
        echo "    environment  : ${mode}, key4hep latest (today: ${resolved};"
        echo "                   re-resolved every time the kernel starts)"
    else
        echo "    environment  : ${mode}, key4hep ${release} (fixed)"
    fi
    echo
    echo "Start a server with  ./jupyter_kernel.sh lab  (or use your own Jupyter)"
    echo "and pick '${display}' in the launcher / Kernel > Change kernel."
}

cmd_remove() {
    local name="${1:-siwecal-key4hep}"
    local kernel_dir="${JUPYTER_KERNELS_DIR}/${name}"
    if [ ! -d "${kernel_dir}" ]; then
        echo "No such kernel: ${kernel_dir}" >&2
        exit 1
    fi
    rm -rf "${kernel_dir}"
    echo "==> removed ${kernel_dir}"
}

# `lab` and `list` default to latest as well, so a server started this way hosts
# the default kernel without the clean-environment restart in cmd_launch.
cmd_list() {
    export KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-latest}"
    load_env "full"
    exec jupyter kernelspec list
}

cmd_lab() {
    export KEY4HEP_RELEASE="${KEY4HEP_RELEASE:-latest}"
    load_env "full"
    # --no-browser because this normally runs on a headless work node; reach it
    # over an SSH tunnel, same as the event viewer (see README).
    exec jupyter lab --no-browser --notebook-dir="${REPO_ROOT}" "$@"
}

# `launch` is the kernel's own entry point: Jupyter calls it, humans do not.
cmd_launch() {
    local mode="full" release="" args=()
    while [ $# -gt 0 ]; do
        case "$1" in
            --mode)    mode="$2"; shift ;;
            --release) release="$2"; shift ;;
            *)         args+=("$1") ;;
        esac
        shift
    done
    if [ -n "${release}" ]; then export KEY4HEP_RELEASE="${release}"; fi

    # A kernel inherits the environment of the Jupyter server that spawns it, so
    # a server started from a shell with a *different* key4hep would poison it
    # (and key4hep cannot be swapped inside a live process). Rather than die
    # with "already has a different key4hep loaded", start over from a clean
    # environment, keeping only what a kernel legitimately needs: the account,
    # a usable PATH, locale, Kerberos/X11/grid credentials and Jupyter's own
    # settings. K4_KERNEL_CLEAN stops this from ever looping.
    # shellcheck disable=SC1091
    source "${REPO_ROOT}/key4hep_release.sh"
    if [ -n "${KEY4HEP_STACK:-}" ] && [ -z "${K4_KERNEL_CLEAN:-}" ] \
       && [ -n "${release}" ] && ! k4_stack_is "${release}"; then
        local keep=(K4_KERNEL_CLEAN=1 "PATH=/usr/bin:/bin" "HOME=${HOME}")
        local v
        for v in USER LOGNAME SHELL TERM LANG TMPDIR DISPLAY XAUTHORITY \
                 KRB5CCNAME X509_USER_PROXY X509_CERT_DIR \
                 JUPYTER_DATA_DIR JUPYTER_CONFIG_DIR JUPYTER_RUNTIME_DIR; do
            if [ -n "${!v:-}" ]; then keep+=("${v}=${!v}"); fi
        done
        echo "[kernel] inherited key4hep ${KEY4HEP_STACK} != ${release};" \
             "restarting in a clean environment" >&2
        exec env -i "${keep[@]}" "${SELF}" launch \
             --mode "${mode}" --release "${release}" "${args[@]}"
    fi

    load_env "${mode}"
    # -Xfrozen_modules=off: silences the frozen-modules warning and lets the
    # notebook debugger step into stdlib frames.
    exec python -Xfrozen_modules=off -m ipykernel_launcher "${args[@]}"
}

case "${1:-}" in
    install) shift; cmd_install "$@" ;;
    remove)  shift; cmd_remove "$@" ;;
    list)    shift; cmd_list "$@" ;;
    lab)     shift; cmd_lab "$@" ;;
    launch)  shift; cmd_launch "$@" ;;
    -h|--help|"") usage ;;
    *) echo "Unknown command: $1" >&2; usage >&2; exit 2 ;;
esac
