"""Visualisation layer: Plotly figure builders (3-D scene, 2-D layers, distributions)."""

from .scene3d import DetectorScene3D
from .layers2d import LayerGrid2D
from .distributions import DistributionPlots

__all__ = ["DetectorScene3D", "LayerGrid2D", "DistributionPlots"]
