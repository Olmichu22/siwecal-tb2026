"""
Shared helpers for the calibration batch scripts (``run_calibration_batch.py``
and ``list_run_thresholds.py``): parsing a "list of run folders" CLI argument
and looking up each run's ThresholdDAC from its ``Run_Settings.txt``.

Kept as a plain sibling module (not a package) since both consumers live
directly in ``gaudi_jobs/`` and Python already adds a running script's own
directory to ``sys.path``, so ``import calib_run_utils`` resolves without
any extra packaging.
"""
import glob
import os

from siwecal_eventbuilder.run_settings import read_threshold_dac

# Read-only source of raw DAQ chunk files (per-run subdirectories). NEVER
# write or delete anything here (see run_calibration_batch.py SAFETY note).
DEFAULT_RAW_BASE = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/Data/rundata"


def parse_run_folder_list(spec, default_base=DEFAULT_RAW_BASE):
    """Parse a comma-separated ``--runs`` value into ``[(run_name, raw_dir), ...]``.

    Each entry can be either a full path to a run's raw data folder (contains
    a path separator) or a bare run name (e.g. ``TB2026CERN_run_000060``),
    which is resolved as ``<default_base>/<name>``. Mixing both styles in one
    list is fine. Order is preserved (duplicates are not removed: repeating a
    run just accumulates its statistics twice, which is harmless but the
    caller may want to dedupe upstream if that matters).
    """
    entries = [e.strip() for e in spec.split(",") if e.strip()]
    if not entries:
        raise SystemExit("ERROR: --runs is empty")
    resolved = []
    for entry in entries:
        if os.sep in entry:
            raw_dir = entry.rstrip(os.sep)
            run_name = os.path.basename(raw_dir)
        else:
            run_name = entry
            raw_dir = os.path.join(default_base, entry)
        resolved.append((run_name, raw_dir))
    return resolved


def run_threshold(raw_dir):
    """ThresholdDAC (mode across the file) from ``<raw_dir>/Run_Settings.txt``.

    Returns -1 if the file is missing/unreadable or has no ThresholdDAC
    token (see read_threshold_dac).
    """
    return read_threshold_dac(os.path.join(raw_dir, "Run_Settings.txt"))


def raw_chunk_files(raw_dir, run):
    """Sorted list of a run's raw DAQ chunk files.

    Two naming conventions coexist in the raw data area depending on the
    run's trigger mode: plain runs use ``<run>_raw.bin``/``<run>_raw.bin_NNNN``,
    but ``..._eudaq_run_...`` runs use ``<run>.bin``/``<run>.bin_NNNN`` (no
    ``_raw`` infix) -- e.g. TB2026CERN_eudaq_run_000147..000151. Try the
    plain pattern first and fall back to the eudaq one; without this, every
    eudaq run raises "no raw chunk files found" even though its data is
    right there under a differently-shaped filename.

    The first chunk has no numeric suffix, later ones do (_0001, _0002,
    ...); lexicographic sort already orders them correctly since the bare
    filename is a strict prefix of the suffixed ones and the suffixes are
    zero-padded.
    """
    pattern = os.path.join(raw_dir, f"{run}_raw.bin*")
    files = sorted(glob.glob(pattern))
    if files:
        return files
    eudaq_pattern = os.path.join(raw_dir, f"{run}.bin*")
    files = sorted(glob.glob(eudaq_pattern))
    if not files:
        raise SystemExit(f"ERROR: no raw chunk files found for {run}: tried {pattern} and {eudaq_pattern}")
    return files


def discover_run_folders(base_dir=DEFAULT_RAW_BASE):
    """``[(run_name, raw_dir), ...]`` for every immediate subdirectory of
    base_dir, sorted by name."""
    if not os.path.isdir(base_dir):
        raise SystemExit(f"ERROR: raw base directory not found: {base_dir}")
    names = sorted(
        name for name in os.listdir(base_dir)
        if os.path.isdir(os.path.join(base_dir, name))
    )
    return [(name, os.path.join(base_dir, name)) for name in names]
