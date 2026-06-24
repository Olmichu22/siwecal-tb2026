"""
siwecal_validation
==================

Object-oriented validation plots for the SiW-ECAL reconstructed ``ecal`` tree.

Pipeline (one class per concern)
--------------------------------
* :class:`~siwecal_validation.config.PlotConfig`   -- tunables / constants.
* :class:`~siwecal_validation.event_data.EventData` -- per-event arrays + I/O.
* :class:`~siwecal_validation.selection.CutSet`     -- selection cuts (general
  or per-energy).
* :class:`~siwecal_validation.plots.Plotter`        -- one plot type each
  (subclass + register to add a new plot).
* :class:`~siwecal_validation.output.OutputLayout`  -- structured output tree.
* :class:`~siwecal_validation.results.ResultsWriter` -- signal-rate results table.
* :class:`~siwecal_validation.runner.ValidationRunner` -- orchestrates it all.
"""

from .config import PlotConfig, energy_to_float
from .event_data import EventData
from .fits import fit_gaussian, gauss
from .output import OutputLayout
from .plots import (DEFAULT_PLOTTERS, EnergyHistPlotter, MipLikenessPlotter,
                    NhitBarycenterPlotter, PlotContext, Plotter, plot_grid,
                    plot_type_grids)
from .results import ResultsWriter
from .runner import ValidationRunner
from .selection import CutSet

__all__ = [
    "PlotConfig", "energy_to_float", "EventData", "fit_gaussian", "gauss",
    "OutputLayout", "DEFAULT_PLOTTERS", "EnergyHistPlotter",
    "MipLikenessPlotter", "NhitBarycenterPlotter", "PlotContext", "Plotter",
    "plot_grid", "plot_type_grids", "ResultsWriter", "ValidationRunner",
    "CutSet",
]
