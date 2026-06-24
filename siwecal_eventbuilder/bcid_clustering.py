"""
BCID clustering: turn one acquisition into a list of candidate time windows.

Physics background
------------------
Every triggered SCA carries a *bunch-crossing ID* (BCID): a counter that ticks
with the machine clock and tells *when* the cell was filled. Channels belonging
to the same physical particle shower are read out within a few BCID counts of
each other, but not necessarily at the exact same value. The job of this stage
is therefore to **group nearby BCIDs across all chips into time windows**, each
window being one candidate event.

Two subtleties handled here:

* **Overflow.** The BCID counter is 12-bit and wraps every 4096 counts. We
  cluster on the *raw* BCID (modulo 4096) and reconstruct the true, unwrapped
  timestamp afterwards -- exactly as tkamiyam's reference does. We deliberately
  do NOT use the converter's ``corrected_bcid`` branch: its overflow correction
  is computed only over the (frequently wrong) ``nColumns`` range and injects
  spurious +4096 offsets that fragment real events.

* **min_slabs_hit.** A genuine high-energy shower lights up many slabs. Windows
  spanning fewer than ``min_slabs_hit`` distinct slabs are discarded as noise.
"""

import numpy as np

from .config import BuilderConfig
from .geometry import DetectorGeometry
from .models import BcidWindow
from .root_io import Acquisition


class BcidClusterer:
    """Groups the BCIDs of one acquisition into :class:`BcidWindow` objects."""

    def __init__(self, config: BuilderConfig, geometry: DetectorGeometry):
        self._config = config
        self._geometry = geometry

    # ------------------------------------------------------------------ API ---

    def find_windows(self, acquisition: Acquisition) -> list:
        """Return the accepted BCID windows of ``acquisition``.

        Steps:
        1. Build the overflow-adjustment map from the full raw-BCID matrix.
        2. Collect, per chip, the valid raw BCIDs (after quality cuts).
        3. Merge nearby BCIDs into windows.
        4. Keep windows that span at least ``min_slabs_hit`` slabs and attach the
           overflow-corrected timestamp.
        """
        raw_matrix = acquisition.raw_bcid_matrix()
        overflow_map = self._build_overflow_map(raw_matrix)
        chip_to_scas_bcid = self._collect_valid_bcids(acquisition, raw_matrix)
        if not chip_to_scas_bcid:
            return []

        windows = []
        for start_raw, stop_raw in self._merge_into_windows(chip_to_scas_bcid):
            active = self._chips_in_window(chip_to_scas_bcid, start_raw, stop_raw)
            n_slabs = len({slab for (slab, _chip) in active})
            if n_slabs < self._config.min_slabs_hit:
                continue
            windows.append(BcidWindow(
                bcid_label=self._overflow_corrected_label(start_raw, stop_raw, overflow_map),
                start_raw=start_raw,
                stop_raw=stop_raw,
                chip_to_scas=active,
            ))
        return windows

    # ------------------------------------------------------ internal steps ---

    def _build_overflow_map(self, raw_matrix: np.ndarray) -> dict:
        """Map ``raw_bcid -> overflow-adjusted bcid`` for the whole acquisition.

        Within a chip's memory the BCID must increase from one SCA to the next;
        a decrease signals that the counter wrapped. We count the wraps per chip
        (``n_cycles``), find the SCAs where a *new* cycle begins, and record the
        unwrapped value ``raw + 4096 * n_cycles`` for each such BCID.

        Faithful re-implementation of
        ``BCIDHandler._calculate_overflow_correction``. The map stays empty
        unless the acquisition contains more than one overflow transition.
        """
        bad_value = self._config.bad_value
        overflow = self._config.bcid_overflow

        step_decreases = raw_matrix[:, 1:] - raw_matrix[:, :-1] < 0
        n_cycles = np.cumsum(step_decreases, axis=1)

        # A cell starts a *new* cycle when its cumulative wrap count differs from
        # the previous cell's. The very first step (index 0) is never counted, and
        # cells holding the bad-value sentinel are excluded.
        starts_new_cycle = np.zeros_like(n_cycles, dtype=bool)
        starts_new_cycle[:, 1:] = n_cycles[:, 1:] != n_cycles[:, :-1]
        starts_new_cycle[raw_matrix[:, 1:] == bad_value] = False

        overflow_map = {}
        if starts_new_cycle.sum() > 1:
            adjusted_values = np.unique(
                raw_matrix[:, 1:][starts_new_cycle]
                + overflow * n_cycles[starts_new_cycle]
            )
            for adjusted in adjusted_values:
                overflow_map[int(adjusted) % overflow] = int(adjusted)
        return overflow_map

    def _collect_valid_bcids(self, acquisition: Acquisition, raw_matrix: np.ndarray) -> dict:
        """Gather the surviving raw BCID of every ``(slab, chip)`` SCA.

        A cell is kept only if it has a configured chip, is non-empty
        (``nhits != 0``), and its raw BCID passes the start/drop cuts. We iterate
        over *all* SCAs (the reference ignores ``nColumns``).

        Returns ``{(slab, chip): {sca: raw_bcid}}``.
        """
        config = self._config
        n_chips = self._geometry.n_chips_per_slab
        n_scas = self._geometry.n_scas_per_chip

        chip_to_scas_bcid = {}
        for slab in range(acquisition.n_slab_positions):
            for chip in range(n_chips):
                if acquisition.chip_id(slab, chip) < 0:
                    continue
                for sca in range(n_scas):
                    n_hits = acquisition.n_hits(slab, chip, sca)
                    if n_hits == 0 or n_hits > config.max_hits_per_sca:
                        continue
                    raw_bcid = int(raw_matrix[slab * n_chips + chip, sca])
                    if (raw_bcid < 0
                            or raw_bcid < config.skip_bcid_start
                            or raw_bcid in config.drop_bcids):
                        continue
                    chip_to_scas_bcid.setdefault((slab, chip), {})[sca] = raw_bcid
        return chip_to_scas_bcid

    def _merge_into_windows(self, chip_to_scas_bcid: dict) -> list:
        """Greedy single-pass merge of all raw BCIDs into ``(start, stop)`` windows.

        Sorted unique BCIDs are walked left to right; a new window opens whenever
        the gap to the current window's end reaches ``merge_delta``.
        """
        unique_bcids = sorted({
            raw_bcid
            for scas in chip_to_scas_bcid.values()
            for raw_bcid in scas.values()
        })
        if not unique_bcids:
            return []

        merge_delta = self._config.merge_delta
        windows = []
        start = stop = unique_bcids[0]
        for bcid in unique_bcids[1:]:
            if bcid - stop < merge_delta:
                stop = bcid
            else:
                windows.append((start, stop))
                start = stop = bcid
        windows.append((start, stop))
        return windows

    @staticmethod
    def _chips_in_window(chip_to_scas_bcid: dict, start_raw: int, stop_raw: int) -> dict:
        """Restrict ``chip_to_scas_bcid`` to the SCAs inside ``[start, stop]``."""
        active = {}
        for chip_key, scas in chip_to_scas_bcid.items():
            matched = [sca for sca, raw_bcid in scas.items()
                       if start_raw <= raw_bcid <= stop_raw]
            if matched:
                active[chip_key] = matched
        return active

    def _overflow_corrected_label(self, start_raw: int, stop_raw: int, overflow_map: dict) -> int:
        """Unwrap the window's start BCID into an absolute event timestamp.

        Looks up every raw BCID in the window, takes the largest unwrapped value,
        derives how many overflow cycles it represents, and shifts the window
        start by that many cycles. Matches the reference's ``spill_bcids_overflown``.
        """
        overflow = self._config.bcid_overflow
        unwrapped_max = max(overflow_map.get(bcid, bcid)
                            for bcid in range(start_raw, stop_raw + 1))
        n_cycles = unwrapped_max // overflow
        return start_raw + n_cycles * overflow
