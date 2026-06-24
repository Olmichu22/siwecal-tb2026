"""
Validation runner: ties together loading, selection, plotting and results.

* :meth:`ValidationRunner.run_sample` processes one ROOT file (a run or an
  energy point): load -> select -> run every plotter -> record a results row.
* :meth:`ValidationRunner.run_all` iterates the ``event_data`` map of a YAML,
  applying per-energy cuts on top of the general ones, then writes the aggregate
  energy-calibration / resolution plots and the results table.
"""

import os

import numpy as np

from siwecal_common import paths

from .config import PlotConfig, energy_to_float
from .event_data import EventData
from .output import OutputLayout
from .plots import (DEFAULT_PLOTTERS, PlotContext, plot_energy_calibration,
                    plot_grid, plot_resolution, plot_type_grids)
from .results import ResultsWriter
from . import vars_cache
from .selection import CutSet


def _summary_metrics(selected) -> dict:
    """Per-sample aggregate discrimination metrics for the results table.

    Averages over the *selected* events: shower fraction, mean tungsten-weighted
    energy, mean Molière radius (over showering events only) and mean shower-start
    layer. Returns ``None`` for a quantity with no valid entries.
    """
    if len(selected) == 0:
        return {}
    moliere = selected.moliere[selected.is_shower]
    moliere = moliere[np.isfinite(moliere)]
    starts = selected.shower_start[np.isfinite(selected.shower_start)]
    return {
        "shower_frac": float(selected.is_shower.mean()),
        "mean_weighte": float(np.mean(selected.weighte)),
        "mean_moliere": float(np.mean(moliere)) if moliere.size else None,
        "mean_shower_start": float(np.mean(starts)) if starts.size else None,
    }


class ValidationRunner:
    """Coordinates the whole validation for one or many samples."""

    def __init__(self, layout: OutputLayout, config: PlotConfig = None,
                 plotters=DEFAULT_PLOTTERS, make_individual=True,
                 make_grid=True, create_tree=False, save_tree=False,
                 cache_dir=None):
        self._layout = layout
        self._config = config or PlotConfig()
        self._plotters = plotters
        self._make_individual = make_individual
        self._make_grid = make_grid
        self._create_tree = create_tree   # force (re)generate the metrics cache
        self._save_tree = save_tree        # also write a cut-applied augmented tree
        self._cache_dir = cache_dir        # redirect *.valcache.root here (or None)
        self._results = ResultsWriter()
        # One run id per invocation, shared by the summary plots and the
        # results table so they stay associated (and never overwrite a
        # previous run with a different cut).
        self._run_id = layout.allocate_run_id()
        print(f"[Run id] {layout.id_token(self._run_id)}")

    def _write_cut_tree(self, events_path, label, data, cutset) -> None:
        """``--save-tree``: augmented tree of the cut-passing events.

        Written to ``<out>/<label>/trees/<label><suffix>.root`` (the cut suffix
        in the name). Best-effort: a write failure only logs a warning.
        """
        mask = cutset.mask(data)
        out_path = os.path.join(self._layout.label_dir(label, "trees"),
                                f"{label}{cutset.suffix}.root")
        try:
            vars_cache.write(events_path, out_path, self._config, data,
                             row_mask=mask)
            print(f"[SaveTree] wrote {int(mask.sum())} cut-passing event(s) "
                  f"-> {out_path}")
        except (OSError, RuntimeError, ValueError) as error:
            print(f"WARNING: could not write --save-tree {out_path}: {error}")

    # ----------------------------------------------------------- one sample -
    def run_sample(self, events_path, label, cutset: CutSet = None,
                   energy_gev=None, collect=None) -> dict:
        """Process a single events file and return its results row.

        If ``collect`` (a list) is given, append ``(ctx, selected_data)`` for
        this sample so the caller can build cross-sample grids afterwards.
        """
        if not os.path.exists(events_path):
            print(f"ERROR: events file not found: {events_path}")
            return {}

        cutset = cutset or CutSet()
        print(f"\nReading: {events_path}")
        data = EventData.from_root(events_path, label, self._config,
                                   create_tree=self._create_tree,
                                   cache_dir=self._cache_dir)
        n_total = len(data)
        print(f"Valid events: {n_total}")

        selected = data.select(cutset)
        n_selected = len(selected)
        if not cutset.is_empty:
            print(f"Events after cuts ({cutset.label.strip()}): {n_selected}")
        if self._save_tree:
            self._write_cut_tree(events_path, label, data, cutset)
        if n_selected == 0:
            print("WARNING: no events pass the cuts; skipping plots.")
            self._results.add(label=label, energy_gev=energy_gev,
                              n_total=n_total, n_selected=0,
                              signal_rate=0.0,
                              cuts=cutset.label.replace("\n", " ").strip())
            return {}

        ctx = PlotContext(layout=self._layout, label=label, suffix=cutset.suffix,
                          cut_label=cutset.label, config=self._config)
        if collect is not None:
            collect.append((ctx, selected))
        row = {}
        if self._make_individual:
            for plotter in self._plotters:
                result = plotter.make(selected, ctx)
                if result:
                    row.update(result)
        if self._make_grid:
            grid_row = plot_grid(self._plotters, selected, ctx)
            if not self._make_individual:   # grid is the only source of results
                row.update(grid_row)

        self._results.add(
            label=label, energy_gev=energy_gev,
            n_total=n_total, n_selected=n_selected,
            signal_rate=n_selected / n_total,
            mu_fit=row.get("mu_fit"), sigma_fit=row.get("sigma_fit"),
            cuts=cutset.label.replace("\n", " ").strip(),
            **_summary_metrics(selected))
        return row

    # ------------------------------------------------------------ many -------
    def run_all(self, event_data_map, base_path, general_cut: CutSet = None):
        """Process every entry of an ``event_data`` map and make summary plots.

        ``event_data_map`` maps an energy label to a dict with at least ``path``
        and optionally ``cuts`` (a per-energy :class:`CutSet` mapping that
        overrides ``general_cut``).
        """
        general_cut = general_cut or CutSet()
        energies, mu_fits, sigma_fits = [], [], []
        mu_errs, sigma_errs = [], []   # fit uncertainties, for the error bars
        samples = []   # (ctx, selected_data) per processed sample, for type grids

        for label, entry in event_data_map.items():
            energy_gev = energy_to_float(label)
            if energy_gev is None:
                print(f"SKIP {label}: no GeV value in label (e.g. muons)")
                continue

            events_path = entry["path"]
            if not os.path.isabs(events_path):
                # Honour the YAML's base_path first, then fall back to the
                # settings.yml data roots (paths.resolve_input).
                candidate = os.path.join(base_path, events_path)
                events_path = (candidate if os.path.exists(candidate)
                               else paths.resolve_input(events_path))
            per_energy_cut = general_cut.merge(CutSet.from_mapping(entry.get("cuts")))

            row = self.run_sample(events_path, label, cutset=per_energy_cut,
                                  energy_gev=energy_gev, collect=samples)
            if row.get("mu_fit") is not None and row.get("sigma_fit") is not None:
                energies.append(energy_gev)
                mu_fits.append(row["mu_fit"])
                sigma_fits.append(row["sigma_fit"])
                mu_errs.append(row.get("mu_err"))
                sigma_errs.append(row.get("sigma_err"))

        if self._make_grid and samples:
            plot_type_grids(self._plotters, samples, self._layout,
                            self._config, self._run_id)

        if energies:
            plot_energy_calibration(
                energies, mu_fits, sigma_fits,
                self._layout.summary_path("measured_e_vs_beam_e", self._run_id),
                dpi=self._config.dpi)
            plot_resolution(
                energies, mu_fits, sigma_fits, mu_errs, sigma_errs,
                self._layout.summary_path("resolution_vs_e", self._run_id),
                dpi=self._config.dpi)

        self.write_results()

    # ----------------------------------------------------------- results -----
    def write_results(self):
        """Write the CSV + text results tables under the output base."""
        self._results.write_csv(
            self._layout.results_path("results.csv", self._run_id))
        self._results.write_txt(
            self._layout.results_path("results.txt", self._run_id))

    @property
    def results(self):
        return self._results
