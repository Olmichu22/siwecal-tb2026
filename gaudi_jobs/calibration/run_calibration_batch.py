#!/usr/bin/env python
"""
Batch driver for the Gaudi raw2root + pedestal/MIP calibration stages.

Chains three ``k4run`` invocations per request:

1. ``options/run_raw2root.py``, ONE k4run per raw chunk (skipped for chunks
   already decoded, unless ``--force-raw2root``): decodes the raw
   ``<run>_raw.bin_NNNN`` files into ``<converted-dir>/<run>/chunks/``, then
   checks the acquisitions add up. It used to decode the whole run in a single
   process, which silently lost ~75% of the acquisitions on some runs -- and a
   calibration fitted on a quarter of the muons is not obviously wrong, just
   quietly noisier. See gaudi_jobs/decode_chunks.py.
2. ``options/run_pedestal_mip.py`` in ``Mode=Pedestal``.
3. ``options/run_pedestal_mip.py`` in ``Mode=Mip`` (self-contained: it
   recomputes its own on-the-fly pedestal internally, see
   PedestalMipCalibrator.cpp -- no separate pedestal file is read).

Two separate k4run invocations per stage (not one combined TopAlg) so a run
can be re-calibrated without re-decoding its raw binary, and so a failure in
one stage doesn't require re-running the others.

Pedestal and MIP calibration are computed from the SAME data-taking run(s):
the reference tool's single event loop
(SiWECAL-TB-analysis/SLBperformance/DecodedSLBAnalysis.cc::NSlabsAnalysis)
fills the pedestal histograms from samples where hitbit==0 and the MIP
histograms from samples where hitbit==1, in one pass over one input file --
there is no separate "pedestal run" vs. "MIP run" in the deployed
calibration (MuonCalib_it2/{pedestals,mips}/th<N>/ always pairs a Pedestal_*
and a MIP_* file for the same run number). ``--runs`` is therefore a SINGLE
list, used as input for both stages; ALL runs given are pooled into one
combined pedestal calibration and one combined MIP calibration (statistics
across all of them are accumulated together -- see
PedestalMipCalibrator::forEachSelectedChannel, which already loops over
every InputFiles entry into the same histograms -- this maximizes usable
statistics per ASU, matching the deployed cumulative
"..._run_000th<N>_..." MIP files). Use ``--pedestal-only``/``--mip-only`` to
skip generating one of the two tables from the same run list.

``--runs`` accepts a comma-separated list of run FOLDER paths (or bare run
names, resolved against --raw-base/_DEFAULT_RAW_BASE for convenience) -- see
calib_run_utils.parse_run_folder_list. Before decoding anything, every run's
Run_Settings.txt is read to check they all share the same ThresholdDAC; that
shared value fixes the output th<N> subdirectory automatically (override
with --th). Use ``gaudi_jobs/calibration/list_run_thresholds.py`` to see which runs
share a threshold without opening each Run_Settings.txt by hand.

Output layout mirrors the existing ``calibration/MuonCalib_it2*`` directories
(``{pedestals,mips}/th<N>/...``) but under a new ``MuonCalib_gaudi`` directory
by default, so these newly-generated tables never silently overwrite the
existing reference-tool-produced calibration data.

Example::

    python gaudi_jobs/calibration/run_calibration_batch.py \\
        --runs TB2026CERN_run_000060,TB2026CERN_run_000061,TB2026CERN_run_000062

    # explicit folder paths, --th override, pedestal-only
    python gaudi_jobs/calibration/run_calibration_batch.py \\
        --runs /eos/.../rundata/TB2026CERN_run_000060 --th 220 --pedestal-only

SAFETY: the raw data directory (``_DEFAULT_RAW_BASE`` below, .../Data/rundata/)
is READ-ONLY -- this script must never write, overwrite, or delete anything
under it, only read the ``<run>_raw.bin*`` chunk files. Decoded siwecaldecoded
ROOT files always go to a separate ``--converted-dir``
(default ``_DEFAULT_CONVERTED_DIR``, .../Data/rundata_converted_gaudi/).
This script (and any driver like it) must also never delete anything under
.../Data/ -- if a converted file needs regenerating, --force-raw2root
overwrites it in place (TFile RECREATE), it is never removed first.
"""
import argparse
import os
import subprocess
import sys

from siwecal_common import paths

from calib_run_utils import DEFAULT_RAW_BASE as _DEFAULT_RAW_BASE
from calib_run_utils import parse_run_folder_list, run_threshold

# decode_chunks lives one level up, in gaudi_jobs/: the "how to decode a run"
# logic is shared with run_convert_batch.py / run_full_pipeline_batch.py so the
# single-process decode bug cannot come back in one driver but not the others.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from decode_chunks import assert_healthy, chunk_files, chunks_dir, decode_run  # noqa: E402

_OPTIONS_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "gaudi_source", "options")
_CALIB_STEERING = os.path.join(_OPTIONS_DIR, "run_pedestal_mip.py")

# Where decoded siwecaldecoded ROOT files are written. Separate from
# _DEFAULT_RAW_BASE by design (see module docstring SAFETY note).
_DEFAULT_CONVERTED_DIR = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata_converted_gaudi"

# Literal naming convention already used by the deployed cumulative MIP files
# (MuonCalib_it2/mips/th220/MIP_pedestalsubmode1_TB2026CERN_run_000th220_*)
# and expected by siwecal_eventbuilder.cli.resolve_muon_calib_files(th).
# Applied here to both pedestal and MIP outputs whenever >1 run is combined.
_CUMULATIVE_LABEL_PREFIX = "TB2026CERN_run_000th"


def _ensure_chunks(run, raw_dir, converted_dir, force):
    """Decode `run` into <converted_dir>/<run>/chunks/, one k4run per raw chunk.

    This used to decode the whole run in ONE k4run process, into a single
    <run>.root. That silently dropped ~75% of the acquisitions on some runs (see
    gaudi_jobs/decode_chunks.py) -- and a calibration fitted on a quarter of the
    muons is not obviously wrong, just quietly noisier. Now it decodes per chunk
    and verifies the acquisitions add up before anything is fitted.

    Same layout the Condor DAGs write, so an already-decoded run is reused
    instead of decoded a second time.
    """
    chunks = decode_run(run, raw_dir, converted_dir, force=force)
    assert_healthy(chunks, run)
    return chunks


def _run_calibration(mode, input_files, gain, out_path, max_nhit, nslabs_hit, diagnostics):
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    env = {
        **os.environ,
        "CALIB_INPUT_FILES": ",".join(input_files),
        "CALIB_MODE": mode,
        "CALIB_GAIN": gain,
        "CALIB_MAX_NHIT": str(max_nhit),
        "CALIB_NSLABS_HIT": str(nslabs_hit),
    }
    if mode == "Pedestal":
        env["CALIB_OUTPUT_PEDESTAL_FILE"] = out_path
    else:
        env["CALIB_OUTPUT_MIP_FILE"] = out_path
    if diagnostics:
        # Cross-check plots (2D maps + fit-status codes + 1D distributions,
        # see PedestalMipCalibrator.cpp's DiagnosticsFile doc), named after
        # the main output so it's easy to find alongside it.
        diag_path = os.path.splitext(out_path)[0] + ".diagnostics.root"
        env["CALIB_DIAGNOSTICS_FILE"] = diag_path
        print(f"[{mode.lower()}] {gain}-gain -> {out_path} (+ diagnostics: {diag_path})")
    else:
        print(f"[{mode.lower()}] {gain}-gain -> {out_path}")
    result = subprocess.run(["k4run", _CALIB_STEERING], env=env)
    if result.returncode != 0:
        raise SystemExit(f"ERROR: {mode} k4run failed ({out_path})")


def _resolve_siwecaldecoded_files(run_folders, converted_dir, no_raw2root, force_raw2root):
    """Every run's decoded chunks, flattened into one list for CALIB_INPUT_FILES
    (the calibrator chains them, exactly as the Fill stage on Condor does)."""
    files = []
    for run_name, raw_dir in run_folders:
        if no_raw2root:
            existing = chunk_files(converted_dir, run_name)
            if not existing:
                raise SystemExit(f"ERROR: --no-raw2root but no decoded chunks in "
                                 f"{chunks_dir(converted_dir, run_name)}")
            assert_healthy(existing, run_name)
            files.extend(existing)
        else:
            files.extend(_ensure_chunks(run_name, raw_dir, converted_dir, force_raw2root))
    return files


def _check_threshold_consistency(run_folders, th_override):
    """Read ThresholdDAC from each run's Run_Settings.txt (BEFORE any
    decoding) and print it. Returns the threshold label to use for the
    output th<N> subdirectory.

    If th_override is given, the consistency requirement is relaxed to a
    warning (the user explicitly forced a label); otherwise all runs must
    agree or this aborts, since there would be no single sensible th<N> to
    write into.
    """
    print(f"[th-check] Checking ThresholdDAC for {len(run_folders)} run(s) from Run_Settings.txt...")
    per_run = []
    for run_name, raw_dir in run_folders:
        th = run_threshold(raw_dir)
        per_run.append((run_name, th))
        if th == -1:
            print(f"[th-check] {run_name}: ThresholdDAC=UNKNOWN (Run_Settings.txt missing/unreadable at {raw_dir})")
        else:
            print(f"[th-check] {run_name}: ThresholdDAC={th}")

    distinct = sorted({th for _, th in per_run})
    unknown = [name for name, th in per_run if th == -1]

    if th_override is not None:
        if unknown or len(distinct) > 1:
            print(f"[th-check] WARNING: runs do not agree on a single ThresholdDAC ({per_run}); "
                  f"proceeding anyway with the explicit --th {th_override} override.")
        else:
            print(f"[th-check] Detected ThresholdDAC={distinct[0]}; using explicit --th {th_override} override.")
        return th_override

    if unknown:
        raise SystemExit(
            f"ERROR: could not determine ThresholdDAC for run(s) {unknown} (Run_Settings.txt missing/unreadable). "
            "Pass --th to override, or check the run folder with gaudi_jobs/calibration/list_run_thresholds.py."
        )
    if len(distinct) > 1:
        table = "\n".join(f"  {name}: th{th}" for name, th in per_run)
        raise SystemExit(
            "ERROR: --runs mixes different ThresholdDAC values, so there is no single th<N> to write output "
            f"into:\n{table}\n"
            "Split into separate invocations (one per threshold), or pass --th to force a label. "
            "Use gaudi_jobs/calibration/list_run_thresholds.py to group runs by threshold."
        )

    th = distinct[0]
    print(f"[th-check] All {len(run_folders)} run(s) share ThresholdDAC={th} -> output th{th}")
    return str(th)


def _combined_label(run_folders, th):
    """Output-file label: the bare run name for a single run (matches the
    individual Pedestal_<run>_*/MIP_..._<run>_* files already deployed), or
    the literal cumulative-per-threshold label (matches the deployed
    MIP_..._run_000th<N>_* files and resolve_muon_calib_files(th)) when
    pooling statistics from more than one run."""
    if len(run_folders) == 1:
        return run_folders[0][0]
    return f"{_CUMULATIVE_LABEL_PREFIX}{th}"


def main(argv=None):
    p = argparse.ArgumentParser(description="Batch raw2root + pedestal/MIP calibration.")
    p.add_argument("--runs", required=True,
                   help="Comma-separated run folder path(s) (e.g. "
                        f"{_DEFAULT_RAW_BASE}/TB2026CERN_run_000060) or bare run name(s) resolved against "
                        "--raw-base. Used as input for BOTH pedestal and MIP calibration (they come from the "
                        "same data-taking runs); statistics from all given runs are pooled together. "
                        "See gaudi_jobs/calibration/list_run_thresholds.py to find which runs share a threshold.")
    p.add_argument("--th", default=None,
                   help="Threshold label used for the output th<N> subdirectory. Default: auto-derived from "
                        "each run's Run_Settings.txt (ThresholdDAC), which must all agree.")
    p.add_argument("--raw-base", default=_DEFAULT_RAW_BASE,
                   help="Base directory used to resolve bare run names given in --runs "
                        f"(<raw-base>/<run>/). Default: {_DEFAULT_RAW_BASE}")
    p.add_argument("--converted-dir", default=_DEFAULT_CONVERTED_DIR,
                   help="Where the decoded chunks are written: <converted-dir>/<run>/chunks/ "
                        f"(separate from the raw run folders). Default: {_DEFAULT_CONVERTED_DIR}")
    p.add_argument("--outdir", default=None,
                   help="Calibration output base directory. "
                        "Default: settings.yml calib_dir()/MuonCalib_gaudi")
    p.add_argument("--gain", choices=("high", "low", "both"), default="high")
    p.add_argument("--no-raw2root", action="store_true",
                   help="Never decode raw binaries; the decoded chunks must already exist")
    p.add_argument("--force-raw2root", action="store_true",
                   help="Re-decode every raw chunk, even those already decoded")
    stage = p.add_mutually_exclusive_group()
    stage.add_argument("--pedestal-only", action="store_true")
    stage.add_argument("--mip-only", action="store_true")
    p.add_argument("--max-nhit", type=int, default=1)
    p.add_argument("--nslabs-hit", type=int, default=8)
    p.add_argument("--diagnostics", action="store_true",
                    help="Also write a <output>.diagnostics.root cross-check file per calibration table "
                         "(2D layer/chip vs. sca/channel maps, fit-status codes, 1D value distributions)")
    args = p.parse_args(argv)

    do_pedestal = not args.mip_only
    do_mip = not args.pedestal_only

    run_folders = parse_run_folder_list(args.runs, default_base=args.raw_base)
    th = _check_threshold_consistency(run_folders, args.th)
    label = _combined_label(run_folders, th)

    out_base = args.outdir or os.path.join(paths.calib_dir(), "MuonCalib_gaudi")
    gains = ["high", "low"] if args.gain == "both" else [args.gain]

    siwecaldecoded_files = _resolve_siwecaldecoded_files(
        run_folders, args.converted_dir, args.no_raw2root, args.force_raw2root)

    if do_pedestal:
        for gain in gains:
            ped_out = os.path.join(out_base, "pedestals", f"th{th}", f"Pedestal_{label}_{gain}gain.txt")
            _run_calibration("Pedestal", siwecaldecoded_files, gain, ped_out, args.max_nhit, args.nslabs_hit,
                              args.diagnostics)

    if do_mip:
        for gain in gains:
            mip_out = os.path.join(out_base, "mips", f"th{th}",
                                    f"MIP_pedestalsubmode1_{label}_{gain}gain.txt")
            _run_calibration("Mip", siwecaldecoded_files, gain, mip_out, args.max_nhit, args.nslabs_hit,
                              args.diagnostics)

    print("\n[Done]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
