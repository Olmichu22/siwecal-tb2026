"""
Channel calibration: pedestals and MIP scale.

Two numbers are needed to turn a raw ADC count into a physical energy:

* **Pedestal** -- the ADC value a channel reads with *no* signal (electronic
  baseline). It is subtracted from every measurement. Keyed per
  ``(slab_id, chip_id, channel, sca)`` because the baseline depends on the SCA.

* **MIP** -- the most-probable high-gain signal (after pedestal subtraction)
  left by a *minimum-ionising particle*. Dividing the pedestal-subtracted ADC by
  the MIP expresses the deposit in calibrated "MIP" units. Keyed per
  ``(slab_id, chip_id, channel)``.

A :class:`Calibration` can be built three ways:

* :meth:`from_files`   -- read pre-computed text tables (the normal path);
* :meth:`from_data`    -- derive both from raw runs (slow, for commissioning);
* :meth:`disabled`     -- "raw ADC" mode: no subtraction, MIP = 1.
"""

from collections import defaultdict

import numpy as np

from .config import BuilderConfig
from .geometry import DetectorGeometry
from .root_io import AcquisitionReader


class Calibration:
    """Holds pedestal and MIP look-up tables and applies them to channels."""

    def __init__(self, config: BuilderConfig, pedestal_map=None,
                 mip_map=None, default_mip: float = 1.0):
        self._config = config
        self._pedestal_map = pedestal_map      # None => raw mode
        self._mip_map = mip_map                # None => raw mode
        self._default_mip = default_mip

    # ----------------------------------------------------------- accessors ---

    @property
    def enabled(self) -> bool:
        """``True`` unless this is a raw-ADC (disabled) calibration."""
        return self._pedestal_map is not None

    def pedestal(self, slab_id: int, chip_id: int, channel: int, sca: int) -> float:
        """Pedestal for a channel/SCA (0 in raw mode, fallback if unknown)."""
        if not self.enabled:
            return 0.0
        return self._pedestal_map.get(
            (slab_id, chip_id, channel, sca), self._config.pedestal_fallback)

    def mip(self, slab_id: int, chip_id: int, channel: int) -> float:
        """MIP scale for a channel (1 in raw mode, default if uncalibrated)."""
        if not self.enabled:
            return 1.0
        return self._mip_map.get((slab_id, chip_id, channel), self._default_mip)

    # ----------------------------------------------------- constructors -----

    @classmethod
    def disabled(cls, config: BuilderConfig) -> "Calibration":
        """Raw-ADC mode: pedestals = 0, MIP = 1 (no physical calibration)."""
        return cls(config, pedestal_map=None, mip_map=None, default_mip=1.0)

    @classmethod
    def from_files(cls, config: BuilderConfig,
                   pedestal_path: str, mip_path: str) -> "Calibration":
        """Load pedestals and MIPs from the standard text tables."""
        pedestal_map = cls._read_pedestal_file(config, pedestal_path)
        mip_map, default_mip = cls._read_mip_file(config, mip_path)
        return cls(config, pedestal_map, mip_map, default_mip)

    @classmethod
    def from_data(cls, config: BuilderConfig, geometry: DetectorGeometry,
                  pedestal_run_path: str, mip_run_path: str,
                  pedestal_max_entries=None, mip_max_entries=None) -> "Calibration":
        """Compute pedestals and MIPs directly from raw runs (slow path)."""
        pedestal_map = cls._compute_pedestals(
            config, geometry, pedestal_run_path, pedestal_max_entries)
        mip_map, default_mip = cls._compute_mips(
            config, geometry, mip_run_path, pedestal_map, mip_max_entries)
        return cls(config, pedestal_map, mip_map, default_mip)

    # -------------------------------------------------- file loaders --------

    @staticmethod
    def _read_pedestal_file(config: BuilderConfig, path: str) -> dict:
        """Parse ``layer chip channel  mean0 err0 wid0 mean1 ... (15 SCAs)``."""
        pedestal_map = {}
        with open(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                items = line.split()
                layer, chip, channel = int(items[0]), int(items[1]), int(items[2])
                means = items[3::3]
                for sca in range(min(len(means), 15)):
                    pedestal_map[(layer, chip, channel, sca)] = float(means[sca])
        print(f"[Pedestals] Loaded {len(pedestal_map)} entries from {path}")
        return pedestal_map

    @staticmethod
    def _read_mip_file(config: BuilderConfig, path: str):
        """Parse ``layer chip channel mpv ...``; the MPV is the MIP peak."""
        mip_map = {}
        with open(path) as handle:
            for line in handle:
                if line.startswith("#") or not line.strip():
                    continue
                items = line.split()
                layer, chip, channel = int(items[0]), int(items[1]), int(items[2])
                mpv = float(items[3])
                if mpv > 0:
                    mip_map[(layer, chip, channel)] = mpv
        default_mip = float(np.median(list(mip_map.values()))) if mip_map \
            else config.default_mip_fallback
        print(f"[MIP] Loaded {len(mip_map)} entries from {path}")
        print(f"[MIP] Default MIP (median MPV): {default_mip:.1f} ADC counts")
        return mip_map, default_mip

    # --------------------------------------------- compute-from-data --------

    @staticmethod
    def _compute_pedestals(config: BuilderConfig, geometry: DetectorGeometry,
                           path: str, max_entries) -> dict:
        """Estimate per-(slab,chip,channel,sca) pedestals from non-physics SCAs.

        Pedestals are sampled from SCAs flagged ``badbcid != 0`` (retrigger/empty
        cells with no real deposit), then cleaned with a 3-sigma trim.
        """
        print(f"[Pedestals] {path}")
        reader = AcquisitionReader(path, geometry, config.tree_name)
        samples = defaultdict(list)
        n_entries = reader.n_acquisitions if max_entries is None \
            else min(reader.n_acquisitions, max_entries)

        for index in range(n_entries):
            if index % 2000 == 0:
                print(f"  [Pedestals] {index}/{n_entries}", end="\r")
            acquisition = reader.load(index)
            for slab in range(acquisition.n_slab_positions):
                slab_id = acquisition.slboard_id(slab)
                for chip in range(geometry.n_chips_per_slab):
                    chip_id = acquisition.chip_id(slab, chip)
                    if chip_id < 0:
                        continue
                    for sca in range(geometry.n_scas_per_chip):
                        if acquisition.badbcid(slab, chip, sca) == 0:
                            continue  # keep only non-physics cells for the baseline
                        for channel in range(geometry.n_channels_per_chip):
                            adc = acquisition.adc_high(slab, chip, sca, channel)
                            if adc > config.adc_underflow_threshold:
                                samples[(slab_id, chip_id, channel, sca)].append(adc)
        reader.close()

        print(f"\n[Pedestals] Computing stats for {len(samples)} keys ...")
        pedestal_map = {}
        for key, values in samples.items():
            array = np.array(values, dtype=np.float32)
            mean, sigma = array.mean(), array.std()
            if sigma > 0:
                array = array[np.abs(array - mean) < 3 * sigma]
            pedestal_map[key] = float(array.mean())
        print(f"[Pedestals] Done. {len(pedestal_map)} channels.")
        return pedestal_map

    @staticmethod
    def _compute_mips(config: BuilderConfig, geometry: DetectorGeometry,
                      path: str, pedestal_map: dict, max_entries):
        """Estimate the MIP peak per channel as the mode of its signal spectrum."""
        print(f"[MIP] {path}")
        reader = AcquisitionReader(path, geometry, config.tree_name)
        signal_samples = defaultdict(list)
        n_entries = reader.n_acquisitions if max_entries is None \
            else min(reader.n_acquisitions, max_entries)

        for index in range(n_entries):
            if index % 50000 == 0:
                print(f"  [MIP] {index}/{n_entries}", end="\r")
            acquisition = reader.load(index)
            for slab in range(acquisition.n_slab_positions):
                slab_id = acquisition.slboard_id(slab)
                for chip in range(geometry.n_chips_per_slab):
                    chip_id = acquisition.chip_id(slab, chip)
                    if chip_id < 0:
                        continue
                    for sca in range(geometry.n_scas_per_chip):
                        if acquisition.badbcid(slab, chip, sca) > config.badbcid_max_good:
                            continue
                        n_hits = acquisition.n_hits(slab, chip, sca)
                        if n_hits == 0 or n_hits > config.max_hits_per_sca:
                            continue
                        for channel in range(geometry.n_channels_per_chip):
                            if not acquisition.hitbit_high(slab, chip, sca, channel):
                                continue
                            adc = acquisition.adc_high(slab, chip, sca, channel)
                            if adc <= config.adc_underflow_threshold:
                                continue
                            pedestal = pedestal_map.get(
                                (slab_id, chip_id, channel, sca), config.pedestal_fallback)
                            signal = adc - pedestal
                            if signal > 5:
                                signal_samples[(slab_id, chip_id, channel)].append(signal)
        reader.close()

        print(f"\n[MIP] Computing peaks for {len(signal_samples)} channels ...")
        mip_map = {}
        for key, values in signal_samples.items():
            if len(values) < 5:
                continue
            array = np.array(values, dtype=np.float32)
            histogram, edges = np.histogram(array, bins=100, range=(5, 2000))
            peak_bin = int(np.argmax(histogram))
            mip_map[key] = float(0.5 * (edges[peak_bin] + edges[peak_bin + 1]))
        default_mip = float(np.median(list(mip_map.values()))) if mip_map \
            else config.default_mip_fallback
        print(f"[MIP] {len(mip_map)} channels calibrated. Default: {default_mip:.1f}")
        return mip_map, default_mip
