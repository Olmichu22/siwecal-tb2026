"""
Event builder: combine BCID clustering and hit collection for one acquisition.

This is the small orchestration layer that ties the two physics stages together:
for each accepted BCID window it collects the calibrated hits and emits a
:class:`ReconstructedEvent`. It owns no ROOT state, so it is trivial to unit-test
with a fake :class:`Acquisition`.
"""

from .bcid_clustering import BcidClusterer
from .calibration import Calibration
from .config import BuilderConfig
from .geometry import DetectorGeometry
from .hit_collector import HitCollector
from .models import ReconstructedEvent
from .root_io import Acquisition


class EventBuilder:
    """Turns a single :class:`Acquisition` into a list of reconstructed events."""

    def __init__(self, config: BuilderConfig, geometry: DetectorGeometry,
                 calibration: Calibration, pad_map=None):
        self._clusterer = BcidClusterer(config, geometry)
        self._hit_collector = HitCollector(config, geometry, calibration,
                                           pad_map=pad_map)

    def build(self, acquisition: Acquisition) -> list:
        """Return every reconstructed event found in ``acquisition``.

        Windows whose hit collection yields no surviving hits are dropped (this
        is the hit-level safety net: a window can have valid BCIDs yet no channel
        passing the underflow cut).
        """
        events = []
        for window in self._clusterer.find_windows(acquisition):
            hits = self._hit_collector.collect(acquisition, window)
            if not hits:
                continue
            events.append(ReconstructedEvent(bcid=window.bcid_label, hits=hits))
        return events
