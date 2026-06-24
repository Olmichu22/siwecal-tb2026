"""
Hit collection: turn one BCID window into a list of calibrated channel hits.

Given the chips/SCAs that fall inside a :class:`BcidWindow`, this stage looks at
every channel, keeps the ones that actually triggered, calibrates them, and
returns :class:`Hit` objects.

Multi-SCA handling ("merge within chip")
----------------------------------------
A chip may have several SCAs inside the same window. For a given channel we keep
the single SCA whose high-gain discriminator response (``hitbit_high``) is
largest -- i.e. the earliest/strongest trigger -- so each channel contributes at
most one hit per event.
"""

from .calibration import Calibration
from .config import BuilderConfig
from .geometry import DetectorGeometry
from .models import BcidWindow, Hit
from .root_io import Acquisition


class HitCollector:
    """Builds calibrated :class:`Hit` lists from BCID windows."""

    def __init__(self, config: BuilderConfig, geometry: DetectorGeometry,
                 calibration: Calibration, pad_map=None):
        self._config = config
        self._geometry = geometry
        self._calibration = calibration
        self._pad_map = pad_map

    def collect(self, acquisition: Acquisition, window: BcidWindow) -> list:
        """Return all calibrated hits contained in ``window``."""
        hits = []
        for (slab, chip), scas in window.chip_to_scas.items():
            chip_id = acquisition.chip_id(slab, chip)
            if chip_id < 0:
                continue
            slab_id = acquisition.slboard_id(slab)
            for channel, sca in self._best_sca_per_channel(acquisition, slab, chip, scas).items():
                hit = self._build_hit(acquisition, slab, chip, sca, channel, slab_id, chip_id)
                if hit is not None:
                    hits.append(hit)
        return hits

    # ------------------------------------------------------ internal steps ---

    def _best_sca_per_channel(self, acquisition: Acquisition,
                              slab: int, chip: int, scas: list) -> dict:
        """For each triggered channel, pick the SCA with the largest hitbit_high.

        Returns ``{channel: sca}`` for channels that fired in at least one SCA.
        """
        best_sca = {}
        best_response = {}
        for sca in scas:
            for channel in range(self._geometry.n_channels_per_chip):
                response = acquisition.hitbit_high(slab, chip, sca, channel)
                if response <= 0:
                    continue
                if channel not in best_response or response > best_response[channel]:
                    best_response[channel] = response
                    best_sca[channel] = sca
        return best_sca

    def _build_hit(self, acquisition: Acquisition, slab: int, chip: int, sca: int,
                   channel: int, slab_id: int, chip_id: int):
        """Apply the ADC underflow cut and calibration; return a :class:`Hit`.

        Returns ``None`` if the channel fails the high-gain underflow cut.
        """
        adc_high = acquisition.adc_high(slab, chip, sca, channel)
        if adc_high <= self._config.adc_underflow_threshold:
            return None
        adc_low = acquisition.adc_low(slab, chip, sca, channel)

        pedestal = self._calibration.pedestal(slab_id, chip_id, channel, sca)
        mip = self._calibration.mip(slab_id, chip_id, channel)

        adc_high_pedsub = adc_high - pedestal
        if self._calibration.is_masked(slab_id, chip_id, channel):
            energy_mip = 0.0
        else:
            energy_mip = adc_high_pedsub / mip if mip > 0 else 0.0
        # Transverse position from the pad map (keyed by the *geometric* chip
        # index, 0..15); longitudinal position is fixed by the slab.
        if self._pad_map is not None:
            x, y = self._pad_map.position(slab, chip, channel)
        else:
            x = y = float("nan")
        return Hit(
            slab_position=slab,
            chip_id=chip_id,
            channel=channel,
            sca=sca,
            adc_high_pedsub=adc_high_pedsub,
            adc_low_pedsub=float(adc_low - pedestal),
            energy_mip=energy_mip,
            x=x,
            y=y,
            z=self._geometry.slab_z(slab),
        )
