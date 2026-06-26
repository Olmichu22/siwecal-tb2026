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

# Lazy public API (PEP 562): importing this package must NOT pull in the heavy
# I/O / plotting stack (uproot, matplotlib, ...). That keeps the pure-numpy core
# -- notably ``siwecal_validation.metrics`` -- importable (and unit-testable)
# with just numpy installed, so the metrics tests run without the full key4hep
# stack. Each public name is imported from its submodule only on first access.
import importlib

# name -> submodule it lives in.
_EXPORTS = {
    "PlotConfig": ".config",
    "energy_to_float": ".config",
    "EventData": ".event_data",
    "fit_gaussian": ".fits",
    "gauss": ".fits",
    "OutputLayout": ".output",
    "DEFAULT_PLOTTERS": ".plots",
    "EnergyHistPlotter": ".plots",
    "MipLikenessPlotter": ".plots",
    "NhitBarycenterPlotter": ".plots",
    "PlotContext": ".plots",
    "Plotter": ".plots",
    "plot_grid": ".plots",
    "plot_type_grids": ".plots",
    "ResultsWriter": ".results",
    "ValidationRunner": ".runner",
    "CutSet": ".selection",
}

__all__ = list(_EXPORTS)


def __getattr__(name):
    """Resolve a public name on first access (see PEP 562)."""
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(module_name, __name__), name)
    globals()[name] = value          # cache so subsequent access is a plain attr
    return value


def __dir__():
    return sorted(set(globals()) | set(_EXPORTS))
