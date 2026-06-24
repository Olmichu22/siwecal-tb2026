"""
Parallel driver: build one run across several worker processes.

A run is split into contiguous slices of acquisitions, each processed by its own
:class:`~siwecal_eventbuilder.event_builder.EventBuilder` in a separate process,
writing a temporary ROOT file. The temporary files are then merged with ``hadd``.

Multiprocessing note
--------------------
On Linux the default ``fork`` start method lets child processes inherit the
parent's memory copy-on-write. We exploit that to share the (potentially large)
calibration tables without pickling them: the pipeline stashes its components in
module-level globals just before creating the pool, and each worker reads them
back. This mirrors the original script's behaviour.
"""

import multiprocessing
import os
import re
import subprocess

import ROOT

from .calibration import Calibration
from .config import BuilderConfig
from .event_builder import EventBuilder
from .geometry import DetectorGeometry
from .root_io import AcquisitionReader, EcalWriter


# Set by EventBuildingPipeline.build_run before the pool is forked; read by the
# worker processes. Not part of the public API.
_WORKER_CONFIG = None
_WORKER_GEOMETRY = None
_WORKER_CALIBRATION = None
_WORKER_PAD_MAP = None


def _process_chunk(input_path: str, output_path: str,
                   entry_start: int, entry_end: int, worker_id: int,
                   run_id: int = -1) -> int:
    """Build events for acquisitions ``[entry_start, entry_end)`` into one file.

    Runs inside a worker process and relies on the fork-inherited globals for its
    configuration, geometry and calibration. ``run_id`` is stamped on every event
    for provenance. Returns the number of events written.
    """
    ROOT.gROOT.SetBatch(True)
    config, geometry, calibration = _WORKER_CONFIG, _WORKER_GEOMETRY, _WORKER_CALIBRATION
    pad_map = _WORKER_PAD_MAP

    reader = AcquisitionReader(input_path, geometry, config.tree_name)
    builder = EventBuilder(config, geometry, calibration, pad_map=pad_map)
    writer = EcalWriter(output_path, config.max_hits_per_event, run_id=run_id)

    n_written = 0
    n_chunk = entry_end - entry_start
    for local_index, acquisition_index in enumerate(range(entry_start, entry_end)):
        if local_index % 1000 == 0:
            print(f"  [W{worker_id}] {local_index}/{n_chunk}  events: {n_written}",
                  end="\r", flush=True)
        acquisition = reader.load(acquisition_index)
        for event_index, event in enumerate(builder.build(acquisition)):
            if writer.write(event, spill_index=acquisition_index, event_index=event_index):
                n_written += 1

    writer.close()
    reader.close()
    print(f"  [W{worker_id}] done -- {n_written} events in {output_path}")
    return n_written


def _chunk_worker(args) -> tuple:
    """Picklable Pool entry point: unpacks ``args`` and runs :func:`_process_chunk`."""
    input_path, output_path, entry_start, entry_end, worker_id, run_id = args
    return output_path, _process_chunk(input_path, output_path,
                                       entry_start, entry_end, worker_id, run_id)


class EventBuildingPipeline:
    """Coordinates parallel event building and output merging for whole runs."""

    def __init__(self, config: BuilderConfig, geometry: DetectorGeometry,
                 calibration: Calibration, pad_map=None):
        self._config = config
        self._geometry = geometry
        self._calibration = calibration
        self._pad_map = pad_map

    def build_run(self, input_path: str, output_path: str,
                  max_entries=None, n_workers=None, run_id=None) -> int:
        """Build every event of one run and write the merged ``ecal`` file.

        ``run_id`` is the numeric run identifier stamped on every event for
        provenance; if ``None`` it is inferred from the input file name.
        Returns the total number of events written.
        """
        n_workers = n_workers or self._config.default_workers
        if run_id is None:
            run_id = self._run_id_from_path(input_path)

        with AcquisitionReader(input_path, self._geometry, self._config.tree_name) as reader:
            n_total = reader.n_acquisitions
        if max_entries is not None:
            n_total = min(n_total, max_entries)

        n_workers = max(1, min(n_workers, n_total))
        chunk_args = self._partition(input_path, output_path, n_total, n_workers, run_id)

        print(f"[Build] {n_total} acquisitions -> {n_workers} workers "
              f"(~{n_total // n_workers} acq/worker)")
        for _input, tmp_path, start, end, worker_id, _run in chunk_args:
            print(f"  W{worker_id}: [{start}, {end})  ->  {os.path.basename(tmp_path)}")

        # Publish components to the workers via fork-inherited globals.
        global _WORKER_CONFIG, _WORKER_GEOMETRY, _WORKER_CALIBRATION, _WORKER_PAD_MAP
        _WORKER_CONFIG = self._config
        _WORKER_GEOMETRY = self._geometry
        _WORKER_CALIBRATION = self._calibration
        _WORKER_PAD_MAP = self._pad_map

        if n_workers == 1:
            results = [_chunk_worker(chunk_args[0])]
        else:
            with multiprocessing.Pool(processes=n_workers) as pool:
                results = pool.map(_chunk_worker, chunk_args)

        tmp_files = [output for output, _count in results]
        total_events = sum(count for _output, count in results)

        print(f"\n[Merge] {total_events} events total -- merging into {output_path}")
        self._merge(tmp_files, output_path)
        self._cleanup(tmp_files)

        print(f"[Build] Done. {total_events} events -> {output_path}")
        return total_events

    def build_runs(self, input_paths: list, output_path: str,
                   max_entries=None, n_workers=None) -> int:
        """Build several runs and merge **all** their events into one file.

        Each run is built independently into a temporary per-run file (reusing
        :meth:`build_run`), and the per-run files are then concatenated with the
        same ``hadd``-based merge used for worker chunks. Returns the grand total
        of events written across all runs.

        Note
        ----
        The events keep their per-run ``event``/``spill`` numbering, so those
        numbers are no longer unique across the combined file. The physics
        content (one row per reconstructed event) is unaffected.
        """
        if len(input_paths) == 1:
            # Nothing to concatenate: build straight into the final file.
            return self.build_run(input_paths[0], output_path,
                                  max_entries=max_entries, n_workers=n_workers)

        per_run_files = []
        total_events = 0
        for run_position, input_path in enumerate(input_paths):
            part_path = output_path.replace(".root", f"_part{run_position:02d}.root")
            print(f"\n--- Run {run_position + 1}/{len(input_paths)}: "
                  f"{os.path.basename(input_path)} ---")
            total_events += self.build_run(input_path, part_path,
                                           max_entries=max_entries, n_workers=n_workers)
            per_run_files.append(part_path)

        print(f"\n[Concat] {total_events} events across {len(per_run_files)} run(s) "
              f"-- merging into {output_path}")
        self._merge(per_run_files, output_path)
        self._cleanup(per_run_files)

        print(f"[Concat] Done. {total_events} events -> {output_path}")
        return total_events

    # ------------------------------------------------------ internal steps ---

    @staticmethod
    def _partition(input_path, output_path, n_total, n_workers, run_id=-1) -> list:
        """Split ``n_total`` acquisitions into contiguous per-worker slices."""
        chunk_size = n_total // n_workers
        chunk_args = []
        for worker_id in range(n_workers):
            start = worker_id * chunk_size
            end = start + chunk_size if worker_id < n_workers - 1 else n_total
            tmp_path = output_path.replace(".root", f"_chunk{worker_id:02d}.root")
            chunk_args.append((input_path, tmp_path, start, end, worker_id, run_id))
        return chunk_args

    @staticmethod
    def _run_id_from_path(input_path: str) -> int:
        """Extract the numeric run id from a path like ``..._run_000007.root``.

        Returns ``-1`` when no ``run_<number>`` token can be parsed, so the
        provenance branch always holds a defined value.
        """
        match = re.search(r"run_(\d+)", os.path.basename(input_path))
        return int(match.group(1)) if match else -1

    @staticmethod
    def _merge(tmp_files: list, output_path: str) -> None:
        """Merge the per-worker files into ``output_path`` (hadd, with fallback)."""
        try:
            result = subprocess.run(["hadd", "-f", output_path] + tmp_files,
                                    capture_output=True, text=True, check=True)
            print(result.stdout.strip())
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("[Merge] hadd unavailable -- falling back to TFileMerger")
            merger = ROOT.TFileMerger(False)
            merger.OutputFile(output_path)
            for tmp_file in tmp_files:
                merger.AddFile(tmp_file)
            merger.Merge()

    @staticmethod
    def _cleanup(tmp_files: list) -> None:
        """Delete the temporary per-worker files."""
        for tmp_file in tmp_files:
            try:
                os.unlink(tmp_file)
            except OSError:
                pass
