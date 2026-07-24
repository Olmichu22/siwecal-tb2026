#!/usr/bin/env python
"""
Decode a run's raw chunks to ``siwecaldecoded`` ROOT -- ONE k4run process per
chunk -- and check that nothing was lost.

Why this module exists
----------------------
Decoding a run's raw chunks one-per-process is what makes the Condor fan-out
possible: N chunks decode on N workers instead of one process grinding through
them. That is the ONLY reason. Multi-file decoding in a single process is
correct -- run_000012's 292 chunks decoded in one k4run give 119,211 entries,
exactly what 292 separate processes give.

(For a while this module claimed the single-process path silently dropped ~75%
of the acquisitions. That was wrong: the number came from reading a ROOT file
that was still being written. It does not reproduce. Retracted 2026-07-13.)

assert_healthy() below is kept anyway, because a decode CAN go wrong in ways
nothing else would notice, and the check costs one scan of the chunk headers.
Note what it does and does not see: it compares tree entries against the largest
acqNumber, so it catches a decode that skipped acquisitions -- but not one that
simply stopped early, which stays internally consistent (ratio ~1.0) while
holding a fraction of the run. Against truncation, compare with the raw.

Chunks are the decoded data. They are NOT an intermediate to be merged away:
the calibration Fill reads them, and the event builder chains them
(EcalEventBuilder's ``InputFiles``). See gaudi_jobs/condor/README.md.
"""
import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, ".."))
_RAW2ROOT_STEERING = os.path.join(_REPO, "gaudi_source", "options", "run_raw2root.py")

sys.path.insert(0, os.path.join(_HERE, "calibration"))
from calib_run_utils import raw_chunk_files  # noqa: E402


def chunks_dir(converted_dir, run):
    """``<converted_dir>/<run>/chunks/`` -- where a run's decoded chunks live.

    Matches gaudi_jobs/calibration/condor/condor_common.chunks_dir(), which the
    Condor DAGs use; both must agree or the DAG and the interactive drivers would
    decode into different places.
    """
    return os.path.join(converted_dir, run, "chunks")


def chunk_files(converted_dir, run):
    """A run's already-decoded chunks, in order. Empty list if none."""
    cdir = chunks_dir(converted_dir, run)
    if not os.path.isdir(cdir):
        return []
    return sorted(os.path.join(cdir, n) for n in os.listdir(cdir)
                  if n.startswith("chunk_") and n.endswith(".root"))


def is_valid_chunk(path, tree="siwecaldecoded"):
    """Does this decoded chunk actually OPEN and hold the tree?

    Not `os.path.getsize(path) > 0`. A k4run killed mid-write -- OOM, eviction, an
    XRootD Close that timed out -- leaves a multi-megabyte ROOT file with NO KEYS
    in it. A size check waves that through, and it then survives every stage until
    the calibration Fill finally chokes on it ("Tree 'siwecaldecoded' not found"),
    thousands of jobs and hours later. Ask the file, not the filesystem.
    """
    if not (os.path.exists(path) and os.path.getsize(path) > 0):
        return False
    import ROOT  # local: callers that never decode need no ROOT

    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    handle = ROOT.TFile.Open(path)
    if not handle or handle.IsZombie():
        return False
    ok = bool(handle.Get(tree))
    handle.Close()
    return ok


def assert_not_under(path, readonly_root):
    """Refuse to write anywhere inside the read-only raw data tree."""
    real_path = os.path.realpath(path)
    real_root = os.path.realpath(readonly_root)
    if real_path == real_root or real_path.startswith(real_root + os.sep):
        raise SystemExit(f"REFUSING to write inside the read-only raw data directory: "
                         f"{path} is under {readonly_root}")


def decode_run(run, raw_dir, converted_dir, force=False, dry_run=False,
               tree="siwecaldecoded", max_cycle_jump=10, verbose=True):
    """Decode every raw chunk of `run` into ``<converted_dir>/<run>/chunks/``.

    One k4run per chunk. Idempotent: a chunk whose output already exists is left
    alone unless `force`, so an interrupted decode is resumed by re-running.

    Returns the list of decoded chunk paths. Raises SystemExit on failure.
    """
    assert_not_under(converted_dir, raw_dir)
    raw = raw_chunk_files(raw_dir, run)
    cdir = chunks_dir(converted_dir, run)
    os.makedirs(cdir, exist_ok=True)

    run_settings = os.path.join(raw_dir, "Run_Settings.txt")
    # Whether these chunks use the EUDAQ 0xABCD frame sync (bare "<run>.bin")
    # rather than the plain 0xEEEEEEEE one ("<run>_raw.bin") -- this is a
    # property of which FILE PATTERN raw_chunk_files() actually matched, not
    # of the run's name: some "..._eudaq_run_..." runs (e.g. run 000316) still
    # use plain _raw.bin data. Verified against real data (2026-07-20): even
    # with this flag correctly set, EcalRawDecoder decodes 0 acquisitions for
    # every bare-".bin" eudaq run tried -- those need gaudi_jobs/decode_lcio_runs.py
    # (EcalLcioDecoder reading Data/eudaq/ROC_run_<N>_tp.slcio) instead, not a flag.
    # This still sets the flag correctly (for whatever future run may need it),
    # it just does not make bare-".bin" runs decode via this path today.
    eudaq_format = not os.path.basename(raw[0]).startswith(f"{run}_raw.bin")
    outputs, todo = [], []
    for i, raw_chunk in enumerate(raw):
        out = os.path.join(cdir, f"chunk_{i:04d}.root")
        outputs.append(out)
        if force or not is_valid_chunk(out, tree):
            todo.append((raw_chunk, out))

    if verbose:
        done = len(raw) - len(todo)
        print(f"[decode] {run}: {len(raw)} chunk(s) -> {cdir}/"
              + (f"  ({done} already decoded, {len(todo)} to do)" if done else ""))
    if dry_run:
        return outputs

    for n, (raw_chunk, out) in enumerate(todo, 1):
        if verbose:
            print(f"[decode] {run}: [{n}/{len(todo)}] {os.path.basename(raw_chunk)} "
                  f"-> {os.path.basename(out)}")
        env = {
            **os.environ,
            "RAW_FILES": raw_chunk,
            "RAW2ROOT_OUT": out,
            "RAW2ROOT_TREE": tree,
            "RAW2ROOT_MAX_CYCLE_JUMP": str(max_cycle_jump),
            "RAW2ROOT_RUN_SETTINGS_FILE": run_settings,
        }
        if eudaq_format:
            env["RAW2ROOT_EUDAQ"] = "1"
        # Quiet: one process per chunk means a run's worth of Gaudi banners
        # would otherwise drown the progress lines above.
        result = subprocess.run(["k4run", _RAW2ROOT_STEERING], env=env,
                                stdout=subprocess.DEVNULL if verbose else None)
        if result.returncode != 0:
            raise SystemExit(f"ERROR: raw2root k4run failed on {raw_chunk} (exit {result.returncode})")

    return outputs


def decode_health(paths, tree="siwecaldecoded"):
    """Cheap check that no acquisitions were dropped.

    In a healthy decode there is one tree entry per acquisition, so the largest
    ``acqNumber`` is about the entry count. A ratio well above 1 means acquisitions
    were skipped. (Boundary-split acquisitions push it slightly BELOW 1 -- two
    entries for one acqNumber -- so 0.88 is healthy, not suspicious.)

    Returns ``(entries, max_acq, ratio)``; ratio is NaN when the chain is empty.
    """
    import ROOT  # local: the batch drivers that only decode need no ROOT

    chain = ROOT.TChain(tree)
    for path in paths:
        chain.AddFile(path)
    entries = int(chain.GetEntries())
    if entries == 0:
        return 0, 0, float("nan")
    max_acq = int(chain.GetMaximum("acqNumber"))
    return entries, max_acq, max_acq / float(entries)


def assert_healthy(paths, run, tree="siwecaldecoded", max_ratio=1.5):
    """Fail loudly if a decode looks like it dropped acquisitions.

    Catches a decode that SKIPPED acquisitions. It does NOT catch one that stopped
    early: a truncated file stays internally consistent (entries and acqNumber both
    small, ratio ~1.0) and sails through. Only a comparison against the raw would
    see that. Cheap enough -- one pass over the chunk headers -- to run every time.
    """
    entries, max_acq, ratio = decode_health(paths, tree=tree)
    if entries == 0:
        hint = (" This run's raw data may be EUDAQ-container-only (the raw decoder cannot "
                "parse it at all, regardless of RAW2ROOT_EUDAQ) -- try "
                "gaudi_jobs/decode_lcio_runs.py instead if Data/eudaq/ROC_run_<N>_tp.slcio exists."
                if "eudaq" in run else "")
        raise SystemExit(f"ERROR: {run}: decoded chunks hold 0 entries{hint}")
    if ratio > max_ratio:
        raise SystemExit(
            f"ERROR: {run}: decode looks lossy -- {entries:,} entries but the largest "
            f"acqNumber is {max_acq:,} (ratio {ratio:.2f}, expected ~1.0).\n"
            f"       About {100 * (1 - 1 / ratio):.0f}% of the acquisitions are missing. "
            f"Do NOT reconstruct from these chunks.")
    print(f"[decode] {run}: {entries:,} entries, max acqNumber {max_acq:,} "
          f"(ratio {ratio:.2f}) -- healthy")
    return entries
