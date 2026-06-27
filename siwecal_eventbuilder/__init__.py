"""
siwecal_eventbuilder
====================

Object-oriented SiW-ECAL test-beam event builder for the TB2026CERN data.

It reconstructs physics *events* from the per-acquisition, per-chip data stored
in the converted ROOT files: it groups channel hits that share a bunch-crossing
ID (BCID) time window, applies pedestal/MIP calibration, and writes a flat
``ecal`` tree of reconstructed events.

Pipeline overview (one class per stage)
---------------------------------------
* :class:`~siwecal_eventbuilder.geometry.DetectorGeometry` -- array index math.
* :class:`~siwecal_eventbuilder.config.BuilderConfig`      -- all cut values.
* :class:`~siwecal_eventbuilder.calibration.Calibration`   -- pedestals + MIPs.
* :class:`~siwecal_eventbuilder.root_io.AcquisitionReader` -- read input tree.
* :class:`~siwecal_eventbuilder.bcid_clustering.BcidClusterer` -- BCID windows.
* :class:`~siwecal_eventbuilder.hit_collector.HitCollector` -- calibrated hits.
* :class:`~siwecal_eventbuilder.pad_map.PadMap`             -- pad (x,y) positions.
* :class:`~siwecal_eventbuilder.event_builder.EventBuilder` -- per-acquisition glue.
* :class:`~siwecal_eventbuilder.root_io.EcalWriter`        -- write output tree.
* :class:`~siwecal_eventbuilder.pipeline.EventBuildingPipeline` -- parallel driver.
"""

from .bcid_clustering import BcidClusterer
from .calibration import Calibration
from .config import BuilderConfig
from .event_builder import EventBuilder
from .geometry import DetectorGeometry, W_X0_MM, load_slab_w_thickness_mm
from .hit_collector import HitCollector
from .models import BcidWindow, Hit, ReconstructedEvent
from .pad_map import PadMap, load_mapping_file
from .pipeline import EventBuildingPipeline
from .root_io import Acquisition, AcquisitionReader, EcalWriter
from .run_settings import read_threshold_dac, run_settings_path
from .settings import AppSettings, load_config_file

__all__ = [
    "BcidClusterer",
    "Calibration",
    "BuilderConfig",
    "EventBuilder",
    "DetectorGeometry",
    "W_X0_MM",
    "load_slab_w_thickness_mm",
    "HitCollector",
    "PadMap",
    "load_mapping_file",
    "BcidWindow",
    "Hit",
    "ReconstructedEvent",
    "EventBuildingPipeline",
    "Acquisition",
    "AcquisitionReader",
    "EcalWriter",
    "read_threshold_dac",
    "run_settings_path",
    "AppSettings",
    "load_config_file",
]
