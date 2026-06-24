"""
Plotters for the SiW-ECAL validation.

Each plot is a small class deriving from :class:`Plotter`. A plotter implements
``draw(data, ctx, ax)``, which draws onto a given Matplotlib ``Axes`` and may
return a dict of scalar results (e.g. fitted mean and sigma). The base class
turns that into a standalone PNG via the shared :meth:`Plotter.make` template,
and :func:`plot_grid` reuses the very same ``draw`` to tile every plotter into
one compact combined figure.

To add a new plot you only need to (1) subclass :class:`Plotter`, setting
``category``/``stem``/``title`` and implementing ``draw``, and (2) add it to
:data:`DEFAULT_PLOTTERS`. It then appears both as an individual PNG and as a
panel of the combined grid with no further changes.
"""

import math
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .config import PlotConfig
from .fits import fit_gaussian, gauss
from .output import OutputLayout

# Resolved relative to this package so it works from any working directory.
PLT_STYLE = os.path.join(os.path.dirname(__file__), "MLPConfig", "newams.mplstyle")
plt.style.use(PLT_STYLE)


@dataclass
class PlotContext:
    """Everything a plotter needs besides the event data itself."""

    layout: OutputLayout
    label: str            # run or energy label (sample name)
    suffix: str           # file-name suffix encoding the active cuts
    cut_label: str        # human-readable cut description for titles
    config: PlotConfig

    def out_path(self, category: str, stem: str) -> str:
        """Full PNG path ``<base>/<label>/<category>/<stem>_<label><suffix>.png``."""
        directory = self.layout.label_dir(self.label, category)
        return f"{directory}/{stem}_{self.label}{self.suffix}.png"


def _sample_title(ctx: "PlotContext", data) -> str:
    """Figure-level title shared by the standalone and combined figures."""
    return f"{ctx.label}  —  {len(data)} events{ctx.cut_label}"


class Plotter(ABC):
    """Base class for all validation plots."""

    #: sub-folder name for this plot type.
    category: str = "plot"
    #: file stem for the standalone PNG.
    stem: str = "plot"
    #: short per-panel title (also used inside the combined grid).
    title: str = ""
    #: figure size of the standalone PNG.
    figsize: tuple = (6, 6)

    @abstractmethod
    def draw(self, data, ctx: PlotContext, ax) -> Optional[dict]:
        """Draw this plot onto ``ax``; optionally return scalar results."""

    def make(self, data, ctx: PlotContext) -> Optional[dict]:
        """Standalone PNG template: build a figure, draw, save, return results."""
        fig, ax = plt.subplots(figsize=self.figsize)
        fig.suptitle(_sample_title(ctx, data), fontsize=13, fontweight="bold")
        result = self.draw(data, ctx, ax)
        fig.tight_layout()
        path = ctx.out_path(self.category, self.stem)
        fig.savefig(path, dpi=ctx.config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {path}")
        return result


# --------------------------------------------------------------------------- #
class MipLikenessPlotter(Plotter):
    """1D histogram of the per-event MIP-likeness score."""

    category = "mip_likeness"
    stem = "mip_likeness"
    title = "MIP likeness"

    def draw(self, data, ctx, ax):
        ax.hist(data.mip_likeness, bins=ctx.config.mip_bins)
        ax.set_title(self.title)
        ax.set_xlabel("MipLike")
        ax.set_ylabel("Events")
        ax.set_xlim(0, 1)
        return None


# --------------------------------------------------------------------------- #
class NhitBarycenterPlotter(Plotter):
    """2D histogram of nHit vs the energy-weighted layer barycenter."""

    category = "nhit_vs_zbary"
    stem = "nhit_vs_zbary"
    title = "Nhit vs energy-weighted <Z>"

    def draw(self, data, ctx, ax):
        cfg = ctx.config
        n_max = int(np.percentile(data.nhit, 99)) + 1
        _h, _xe, _ye, img = ax.hist2d(
            data.zbary, data.nhit,
            bins=[cfg.n_layers, cfg.nhit_zbary_bins],
            range=[[0, cfg.n_layers], [0, n_max]], cmap="plasma")
        ax.figure.colorbar(img, ax=ax, label="Events")
        ax.set_xlabel("<layer> barycenter (energy-weighted)", fontsize=11)
        ax.set_ylabel("N hits (total)", fontsize=11)
        ax.set_title(self.title, fontsize=11)

        ax_top = ax.twiny()
        ax_top.set_xlim(ax.get_xlim())
        ax_top.set_xticks(np.arange(0.5, cfg.n_layers))
        ax_top.set_xticklabels([str(i) for i in range(cfg.n_layers)], fontsize=7)
        ax_top.set_xlabel("Layer", fontsize=9)
        ax.set_yticks(np.arange(0, n_max + 1, max(1, n_max // 10)))
        return None


# --------------------------------------------------------------------------- #
class EnergyHistPlotter(Plotter):
    """Energy spectrum with a Gaussian fit of the peak; returns mu and sigma."""

    category = "energy_dist"
    stem = "energy_dist"
    title = "Energy distribution"

    def draw(self, data, ctx, ax):
        cfg = ctx.config
        # Fit the peak; if it does not converge keep the histogram and carry on
        # (one failed sample must not abort a multi-energy run).
        amplitude = mu_fit = sigma_fit = mu_err = sigma_err = None
        centers = None
        try:
            (amplitude, mu_fit, sigma_fit, mu_err, sigma_err,
             centers, _counts) = fit_gaussian(
                data.energy, percentiles=cfg.energy_fit_percentiles,
                bins=cfg.energy_fit_bins)
        except (RuntimeError, ValueError) as error:
            print(f"WARNING: energy fit failed for {ctx.label}: {error}")

        if centers is not None:
            bin_width = (centers[1] - centers[0]) if len(centers) > 1 else 1.0
            full_bins = np.arange(data.energy.min(),
                                  data.energy.max() + bin_width, bin_width)
        else:
            full_bins = cfg.energy_fit_bins
        ax.hist(data.energy, bins=full_bins, alpha=0.6, label="Data")
        # if mu_fit is not None:
            # x_fit = np.linspace(centers[0], centers[-1], 500)
            # ax.plot(x_fit, gauss(x_fit, amplitude, mu_fit, sigma_fit), "r-", lw=2,
                    # label=f"Gaussian fit (μ={mu_fit:.2f}, σ={sigma_fit:.2f})")
        ax.set_title(self.title)
        ax.set_xlabel("Energy")
        ax.set_ylabel("Events")
        ax.legend()
        return {"mu_fit": None if mu_fit is None else float(mu_fit),
                "sigma_fit": None if sigma_fit is None else float(sigma_fit),
                "mu_err": None if mu_err is None else float(mu_err),
                "sigma_err": None if sigma_err is None else float(sigma_err)}


# --------------------------------------- particle-discrimination plotters ---
def _finite(values: np.ndarray) -> np.ndarray:
    """Drop NaNs (e.g. shower-shape variables of non-shower events)."""
    return values[np.isfinite(values)]


class LongitudinalProfilePlotter(Plotter):
    """Mean energy deposited per layer (longitudinal shower profile)."""

    category = "long_profile"
    stem = "long_profile"
    title = "Longitudinal profile  <E>/layer"

    def draw(self, data, ctx, ax):
        cfg = ctx.config
        if len(data) == 0:
            return None
        mean_profile = data.energy_per_layer.mean(axis=0)
        ax.step(np.arange(cfg.n_layers), mean_profile, where="mid")
        ax.set_xlabel("Layer")
        ax.set_ylabel("<E> per layer [MIP]")
        ax.set_title(self.title)
        return None


class WeighteHistPlotter(Plotter):
    """Tungsten-weighted (sampling-corrected) energy Σ E·W/X0."""

    category = "weighte"
    stem = "weighte"
    title = "Weighted energy (Σ E·W/X0)"

    def draw(self, data, ctx, ax):
        ax.hist(data.weighte, bins=ctx.config.weighte_bins)
        ax.set_xlabel("Weighted energy [MIP·X0]")
        ax.set_ylabel("Events")
        ax.set_title(self.title)
        return None


class MoliereRadiusPlotter(Plotter):
    """Molière radius (90% transverse containment) of showering events."""

    category = "moliere"
    stem = "moliere"
    title = "Molière radius"

    def draw(self, data, ctx, ax):
        vals = _finite(data.moliere[data.is_shower])
        ax.hist(vals, bins=ctx.config.moliere_bins)
        ax.set_xlabel("Molière radius [mm]")
        ax.set_ylabel("Shower events")
        ax.set_title(self.title)
        return None


class TransverseRmsPlotter(Plotter):
    """Energy-weighted transverse RMS radius (compactness proxy)."""

    category = "transverse_rms"
    stem = "transverse_rms"
    title = "Transverse RMS radius"

    def draw(self, data, ctx, ax):
        ax.hist(_finite(data.transverse_rms), bins=ctx.config.moliere_bins)
        ax.set_xlabel("Transverse RMS [mm]")
        ax.set_ylabel("Events")
        ax.set_title(self.title)
        return None


class ShowerStartLayerPlotter(Plotter):
    """Layer where the shower starts (non-shower events excluded)."""

    category = "shower_start"
    stem = "shower_start"
    title = "Shower start layer"

    def draw(self, data, ctx, ax):
        ax.hist(_finite(data.shower_start),
                bins=np.arange(0, ctx.config.n_layers + 1))
        ax.set_xlabel("Shower start layer")
        ax.set_ylabel("Shower events")
        ax.set_title(self.title)
        return None


class NLayersHitPlotter(Plotter):
    """Number of layers with at least one hit (longitudinal extent)."""

    category = "n_layers_hit"
    stem = "n_layers_hit"
    title = "N layers hit"

    def draw(self, data, ctx, ax):
        ax.hist(data.n_layers_hit, bins=np.arange(0, ctx.config.n_layers + 2))
        ax.set_xlabel("N layers hit")
        ax.set_ylabel("Events")
        ax.set_title(self.title)
        return None


class DiscriminatorPlotter(Plotter):
    """2-D particle-type discriminator: Molière radius vs shower-start layer."""

    category = "discriminator"
    stem = "discriminator"
    title = "Molière radius vs shower-start layer"

    def draw(self, data, ctx, ax):
        cfg = ctx.config
        keep = data.is_shower & np.isfinite(data.shower_start) \
            & np.isfinite(data.moliere)
        if not keep.any():
            ax.set_title(self.title + " (no showers)")
            return None
        _h, _xe, _ye, img = ax.hist2d(
            data.shower_start[keep], data.moliere[keep],
            bins=[cfg.n_layers, cfg.moliere_bins],
            range=[[0, cfg.n_layers], [0, np.percentile(
                data.moliere[keep], 99) + 1]],
            cmap="plasma")
        ax.figure.colorbar(img, ax=ax, label="Shower events")
        ax.set_xlabel("Shower start layer")
        ax.set_ylabel("Molière radius [mm]")
        ax.set_title(self.title)
        return None


#: The plots produced for every sample. Append a new Plotter here to enable it.
DEFAULT_PLOTTERS = (
    MipLikenessPlotter(),
    NhitBarycenterPlotter(),
    EnergyHistPlotter(),
    LongitudinalProfilePlotter(),
    WeighteHistPlotter(),
    MoliereRadiusPlotter(),
    TransverseRmsPlotter(),
    ShowerStartLayerPlotter(),
    NLayersHitPlotter(),
    DiscriminatorPlotter(),
)


# ------------------------------------------------------ combined (grid) ------
def plot_grid(plotters, data, ctx: PlotContext) -> dict:
    """Tile every plotter into one compact figure and save it.

    Reuses each plotter's ``draw`` so the grid stays in sync automatically.
    Returns the merged results dict (e.g. ``mu_fit``/``sigma_fit``).
    """
    n = len(plotters)
    ncols = max(1, ctx.config.grid_ncols)
    nrows = math.ceil(n / ncols)
    panel = ctx.config.grid_panel_inches
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * panel, nrows * panel),
                             squeeze=False)
    flat = axes.flatten()

    results = {}
    for ax, plotter in zip(flat, plotters):
        result = plotter.draw(data, ctx, ax)
        if result:
            results.update(result)
    for ax in flat[n:]:                      # hide unused cells
        ax.axis("off")

    fig.suptitle(_sample_title(ctx, data), fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = ctx.out_path("combined", "combined")
    fig.savefig(path, dpi=ctx.config.dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {path}")
    return results


def plot_type_grids(plotters, samples, layout, config, run_id):
    """For each plot type, a grid with one panel per sample (energy).

    ``samples`` is a list of ``(ctx, data)`` pairs (one per processed energy).
    Produces one figure per plotter — e.g. all energies' MIP-likeness together,
    all energies' energy distribution together — reusing each plotter's
    ``draw``. Each panel is titled with its sample label; the plot type goes in
    the figure suptitle. Saved under ``<base>/summary/all_<category>_id<NN>.png``.
    """
    n = len(samples)
    if n == 0:
        return
    ncols = max(1, config.grid_ncols)
    nrows = math.ceil(n / ncols)
    panel = config.grid_panel_inches
    for plotter in plotters:
        fig, axes = plt.subplots(nrows, ncols,
                                 figsize=(ncols * panel, nrows * panel),
                                 squeeze=False)
        flat = axes.flatten()
        for ax, (ctx, data) in zip(flat, samples):
            plotter.draw(data, ctx, ax)
            ax.set_title(ctx.label)          # panel title = energy (not plot type)
        for ax in flat[n:]:                  # hide unused cells
            ax.axis("off")
        fig.suptitle(plotter.title, fontsize=14, fontweight="bold")
        fig.tight_layout()
        path = layout.summary_path(f"all_{plotter.category}", run_id)
        fig.savefig(path, dpi=config.dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {path}")


# ------------------------------------------------------ aggregate (summary) --
def plot_energy_calibration(energies, mu_fits, sigma_fits, out_path, dpi=150):
    """Measured (fitted) energy vs true beam energy, with sigma error bars."""
    fig = plt.figure(figsize=(6, 6))
    plt.errorbar(energies, mu_fits, yerr=sigma_fits, capsize=2, fmt="o",
                 linestyle="--")
    plt.title("Measured energy vs beam energy (e- P1)")
    plt.xlabel("Beam energy [GeV]")
    plt.ylabel("Measured energy (fitted μ)")
    plt.grid()
    plt.savefig(out_path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_resolution(energies, mu_fits, sigma_fits, mu_errs, sigma_errs,
                    out_path, dpi=150):
    """Relative energy resolution R = sigma/mu vs true beam energy.

    Statistical error bars come from propagating the fit uncertainties of
    ``sigma`` and ``mu`` (treated as independent) through ``R = sigma / mu``.
    For a ratio the relative errors add in quadrature::

        sigma_R / R = sqrt( (sigma_err / sigma)^2 + (mu_err / mu)^2 )
    =>  sigma_R     = R * sqrt( (sigma_err / sigma)^2 + (mu_err / mu)^2 )

    (sigma and mu of the same fit are in fact correlated; the exact form adds a
    -2 * (1/mu) * (sigma/mu^2) * cov(sigma, mu) term. mu is very well determined
    here, so the sigma term dominates and the independent approximation holds.)
    """
    # float arrays; any missing value (None) becomes NaN and is dropped by
    # matplotlib, so a single failed fit does not break the whole plot.
    mu = np.asarray(mu_fits, dtype=float)
    sigma = np.asarray(sigma_fits, dtype=float)
    mu_err = np.asarray(mu_errs, dtype=float)
    sigma_err = np.asarray(sigma_errs, dtype=float)

    resolution = np.where(mu != 0, sigma / mu, np.nan)
    res_err = resolution * np.sqrt((sigma_err / sigma) ** 2
                                   + (mu_err / mu) ** 2)

    fig = plt.figure(figsize=(8, 6))
    plt.errorbar(energies, resolution, yerr=res_err, fmt="o", capsize=2)
    plt.title("Energy resolution vs beam energy (e- P1)")
    plt.xlabel("Beam energy [GeV]")
    plt.ylabel("Resolution σ/μ")
    plt.grid()
    plt.savefig(out_path, dpi=dpi)
    plt.close(fig)
    print(f"Saved: {out_path}")
