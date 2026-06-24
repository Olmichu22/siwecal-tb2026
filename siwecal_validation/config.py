"""
Configuration constants and tunables for the SiW-ECAL validation plots.

Everything that used to be a module-level magic number in ``plot_val_plots.py``
lives here, so the plotting code never hard-codes a value.
"""

import os
import re
from dataclasses import dataclass

import numpy as np
import yaml

from siwecal_common import paths

# --------------------------------------------------------------- locations ---
# Filesystem locations come from the shared settings.yml (siwecal_common.paths);
# nothing is hard-coded. BASE_PATH is the primary data root, but inputs are
# resolved against all data roots via paths.resolve_input (see cli/runner).
BASE_PATH = paths.data_dir()
DEFAULT_CONFIG = paths.config_file("data_reference_base_event_data.yml")
# Plots/results are written here by default, OUTSIDE the event directories.
DEFAULT_OUTPUT_DIR = paths.output_dir()
DEFAULT_RUN = "TB2026CERN_run_000013"

# Tungsten radiation length [mm]; the per-hit weight is E * W[slab] / X0.
W_X0_MM = 3.5
# Per-slab tungsten thickness [mm], resolved from Tungsten_thickness.yml
# (slab 0 behind 2.8 mm, slabs 1-7 behind 4.2 mm, slabs 8-14 behind 5.6 mm).
W_THICKNESSES_DEFAULT = (2.8,) + (4.2,) * 7 + (5.6,) * 7
# Per-slab physical z position [mm], mirrors
# event_display/conversion/slab_z_positions.yml (also the source of ``hit_z``).
SLAB_Z_MM_DEFAULT = (0.0, 11.0, 22.0, 33.0, 44.0, 55.0, 66.0, 77.0,
                     88.0, 99.0, 110.0, 132.0, 143.0, 154.0, 165.0)


def load_slab_z_mm(path: str) -> tuple:
    """Read the per-slab z positions [mm] from a ``slab_z_positions`` YAML."""
    with open(path) as handle:
        document = yaml.safe_load(handle) or {}
    return tuple(float(z) for z in document.get("slab_z_mm", ()))


def load_w_thicknesses(path: str) -> tuple:
    """Resolve ``Tungsten_thickness.yml`` to a per-slab thickness tuple [mm].

    The YAML lists, in beam order, absorber plates (``t1``/``t2``/``t3``, whose
    thicknesses are under ``W_thick:``) interleaved with silicon slabs
    (``sb0``, ``sb1``, ...). The thickness assigned to a slab is the sum of the
    absorber plate(s) crossed since the previous slab. Returns the thicknesses
    ordered by slab index.
    """
    with open(path) as handle:
        cfg = yaml.safe_load(handle) or {}
    thick = cfg.get("W_thick", {})
    structure = cfg.get("structure", [])
    per_slab, accumulated = {}, 0.0
    for item in structure:
        if item in thick:
            accumulated += float(thick[item])
        elif isinstance(item, str) and item.startswith("sb"):
            per_slab[int(item[2:])] = accumulated
            accumulated = 0.0
    return tuple(per_slab[i] for i in sorted(per_slab))


@dataclass(frozen=True)
class PlotConfig:
    """Immutable bag of plotting/analysis tunables shared by all plotters."""

    n_layers: int = 15
    """Number of ECAL layers (slabs)."""

    layer_z_mm: float = 10.0
    """Approximate z spacing between layers [mm] (for reference axes)."""

    tree_name: str = "ecal"
    """Name of the input TTree holding the reconstructed events."""

    energy_fit_percentiles: tuple = (5, 95)
    """Percentile window of the energy spectrum used for the Gaussian fit."""

    energy_fit_bins: int = 50
    """Number of bins used when fitting the energy peak."""

    mip_bins: int = 100
    """Number of bins for the MIP-likeness histogram."""

    nhit_zbary_bins: int = 50
    """Number of nHit bins in the Nhit-vs-<Z> 2D histogram."""

    dpi: int = 150
    """Resolution of the saved PNG files."""

    grid_ncols: int = 3
    """Number of columns in the combined (grid) figure."""

    grid_panel_inches: float = 5.0
    """Side length [inch] of each panel in the combined (grid) figure."""

    # ----------------------------------------- particle-discrimination metrics
    w_thicknesses: tuple = W_THICKNESSES_DEFAULT
    """Per-slab tungsten thickness [mm] (length must be ``n_layers``)."""

    w_x0_mm: float = W_X0_MM
    """Tungsten radiation length [mm]; per-hit weight is ``E * W[slab] / X0``."""

    slab_z_mm: tuple = SLAB_Z_MM_DEFAULT
    """Per-slab physical z position [mm] (matches the ``hit_z`` written by the
    converter). Not used by the current metrics — kept here, alongside the
    tungsten map, for the event display and any future z-based variable."""

    shower_profile: str = "nhit"
    """Per-layer profile used for the shower flag / start-end: nhit|sume|weighte."""

    shower_e_threshold: float = 5.0
    """Min per-layer activity to count as shower material (raised above the
    noise+MIP floor: with noise a MIP layer has several hits, so 3 is too low)."""

    shower_max_min: float = 10.0
    """Min height of the profile peak for the event to qualify as a shower."""

    shower_start_frac: float = 0.1
    """Fraction of the peak required by the ``*_10`` start/end variants."""

    moliere_containment: float = 0.90
    """Energy fraction defining the Molière containment radius."""

    moliere_bins: int = 60
    """Number of bins for the Molière-radius / transverse-RMS histograms."""

    weighte_bins: int = 50
    """Number of bins for the tungsten-weighted energy histogram."""

    def w_over_x0(self) -> "np.ndarray":
        """Per-slab absorber depth in radiation lengths (``W[slab] / X0``)."""
        return np.asarray(self.w_thicknesses, dtype=float) / self.w_x0_mm


def energy_to_float(label):
    """Convert an energy label to its value in GeV (float).

    The numeric part may use an underscore as the decimal separator, so the
    7.5 GeV point is handled correctly::

        'P1_74GeV'  -> 74.0
        'P1_7_5GeV' -> 7.5

    Returns ``None`` if the label has no GeV value (e.g. 'muons').
    """
    match = re.search(r"_(\d+(?:_\d+)?)GeV", label)
    if not match:
        return None
    return float(match.group(1).replace("_", "."))
