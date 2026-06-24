"""Model layer: detector geometry, single events and per-file datasets."""

from .detector import DetectorModel
from .event import Event
from .dataset import EventDataset

__all__ = ["DetectorModel", "Event", "EventDataset"]
