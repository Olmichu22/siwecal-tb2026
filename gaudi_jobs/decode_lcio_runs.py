#!/usr/bin/env python
"""
Decode EUDAQ LCIO runs (``ROC_run_<N>_tp.slcio``) into the SAME
``siwecaldecoded`` chunk layout the raw decoder produces, via the
``EcalLcioDecoder`` component (gaudi_source/options/run_lcio_decode.py).

Why this exists
---------------
13 of the 14 th210 muon runs (EUDAQ 146-178) exist only as the EUDAQ container
the raw decoder cannot parse -- but they also exist, already decoded, as LCIO in
``Data/eudaq/ROC_run_<N>_tp.slcio``. This is the batch launcher that feeds those
files through EcalLcioDecoder so th210 stops being blocked and enters the
existing event-building / calibration pipeline with no special-casing.

A LCIO run is a SINGLE .slcio file, not chunked, so it decodes to one
``chunk_0000.root`` under the same ``<converted-dir>/<run>/chunks/`` directory a
raw run's chunks go to (decode_chunks.chunks_dir). Downstream -- the calibration
Fill and the event builder's InputFiles -- reads it identically; nothing knows or
cares that it came from LCIO. Idempotent: a run whose chunk already OPENS and
holds the tree is left alone unless --force, so an interrupted batch resumes by
re-running.

SAFETY: the EUDAQ LCIO area and the raw ``rundata`` (Run_Settings.txt) are
READ-ONLY. This script only ever writes under --converted-dir, and refuses to
write inside either read-only tree.

Examples::

    # the seven th210 muon runs (the default set)
    python gaudi_jobs/decode_lcio_runs.py

    # explicit runs, by number or by full run name
    python gaudi_jobs/decode_lcio_runs.py --runs 152,153,TB2026CERN_eudaq_run_000158
"""
import argparse
import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
_LCIO_STEERING = os.path.join(_REPO, "gaudi_source", "options", "run_lcio_decode.py")

sys.path.insert(0, _HERE)
from decode_chunks import assert_not_under, chunks_dir, is_valid_chunk  # noqa: E402

# EUDAQ already-decoded LCIO files (read-only), and the raw rundata tree that
# holds each run's Run_Settings.txt (also read-only).
DEFAULT_EUDAQ_BASE = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/eudaq"
DEFAULT_RAW_BASE = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata"
DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"

# The th210 muon runs that exist only as LCIO (orders: 152-155 and 158-160).
DEFAULT_RUNS = [152, 153, 154, 155, 158, 159, 160]


def run_name_for(spec):
    """Normalise a --runs entry to (run_name, run_number).

    Accepts a bare number (``152``) or a full run name
    (``TB2026CERN_eudaq_run_000152``); returns both the canonical run name and
    the integer run number parsed from it.
    """
    spec = str(spec).strip()
    if spec.isdigit():
        num = int(spec)
        return f"TB2026CERN_eudaq_run_{num:06d}", num
    m = re.search(r"(\d+)\s*$", spec)
    if not m:
        raise SystemExit(f"ERROR: cannot parse a run number from --runs entry: {spec!r}")
    return spec, int(m.group(1))


def lcio_path_for(num, eudaq_base):
    """``<eudaq_base>/ROC_run_<num>_tp.slcio`` (zero-padded to 6 digits)."""
    return os.path.join(eudaq_base, f"ROC_run_{num:06d}_tp.slcio")


def decode_lcio_run(run, num, eudaq_base, raw_base, converted_dir,
                    force=False, dry_run=False, tree="siwecaldecoded", verbose=True):
    """Decode one LCIO run into ``<converted_dir>/<run>/chunks/chunk_0000.root``.

    Returns the output chunk path. Raises SystemExit on failure.
    """
    assert_not_under(converted_dir, eudaq_base)
    assert_not_under(converted_dir, raw_base)

    lcio_in = lcio_path_for(num, eudaq_base)
    if not os.path.exists(lcio_in):
        raise SystemExit(f"ERROR: LCIO file not found for run {run}: {lcio_in}")
    run_settings = os.path.join(raw_base, run, "Run_Settings.txt")
    if not os.path.exists(run_settings):
        # Not fatal: the decoder leaves the settings branches at their -1
        # sentinel, exactly as run_lcio_decode.py documents. Warn and go on.
        if verbose:
            print(f"[lcio] {run}: WARNING no Run_Settings.txt at {run_settings} "
                  f"(threshold/holdDelay branches will be -1)")
        run_settings = ""

    cdir = chunks_dir(converted_dir, run)
    os.makedirs(cdir, exist_ok=True)
    out = os.path.join(cdir, "chunk_0000.root")

    if not force and is_valid_chunk(out, tree):
        if verbose:
            print(f"[lcio] {run}: already decoded -> {out}")
        return out
    if verbose:
        print(f"[lcio] {run}: {os.path.basename(lcio_in)} -> {out}")
    if dry_run:
        return out

    env = {
        **os.environ,
        "LCIO_IN": lcio_in,
        "LCIO_OUT": out,
        "LCIO_TREE": tree,
        "LCIO_RUN_SETTINGS_FILE": run_settings,
    }
    result = subprocess.run(["k4run", _LCIO_STEERING], env=env,
                            stdout=subprocess.DEVNULL if verbose else None)
    if result.returncode != 0:
        raise SystemExit(f"ERROR: run_lcio_decode k4run failed on {lcio_in} (exit {result.returncode})")
    if not is_valid_chunk(out, tree):
        raise SystemExit(f"ERROR: {run} decoded but {out} does not open / lacks the '{tree}' tree")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description="Decode EUDAQ LCIO runs into siwecaldecoded chunks.")
    p.add_argument("--runs", default=None,
                   help="Comma-separated run numbers or names. Default: the seven th210 muon "
                        "runs (152-155, 158-160).")
    p.add_argument("--eudaq-base", default=DEFAULT_EUDAQ_BASE,
                   help=f"Read-only dir of ROC_run_<N>_tp.slcio files. Default: {DEFAULT_EUDAQ_BASE}")
    p.add_argument("--raw-base", default=DEFAULT_RAW_BASE,
                   help=f"Read-only rundata dir holding each run's Run_Settings.txt. "
                        f"Default: {DEFAULT_RAW_BASE}")
    p.add_argument("--converted-dir", default=DEFAULT_CONVERTED_DIR,
                   help=f"Where decoded chunks go. Default: {DEFAULT_CONVERTED_DIR}")
    p.add_argument("--tree", default="siwecaldecoded", help="Output TTree name.")
    p.add_argument("--force", action="store_true", help="Re-decode even if a valid chunk exists.")
    p.add_argument("--dry-run", action="store_true", help="Print what would be decoded; do nothing.")
    args = p.parse_args(argv)

    if args.runs:
        specs = [s for s in args.runs.split(",") if s.strip()]
    else:
        specs = [str(n) for n in DEFAULT_RUNS]

    outputs = []
    for spec in specs:
        run, num = run_name_for(spec)
        outputs.append(decode_lcio_run(run, num, args.eudaq_base, args.raw_base,
                                       args.converted_dir, force=args.force,
                                       dry_run=args.dry_run, tree=args.tree))
    print(f"[lcio] {len(outputs)} run(s) "
          f"{'planned' if args.dry_run else 'decoded'}.")
    return outputs


if __name__ == "__main__":
    main()
