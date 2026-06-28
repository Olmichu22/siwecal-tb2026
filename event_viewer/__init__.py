"""
Interactive event viewer for the SiW-ECAL TB2026CERN prototype.

A self-contained OOP package (read with ``uproot``, drawn with Plotly Dash) to
inspect reconstructed ``ecal`` events one by one, explore file-level distributions
with dynamic cuts, and run simple clustering on per-event variables.

All physics and geometry code this viewer needs is bundled as private ``_``
modules inside this package, so it can run independently without installing or
importing any other siwecal package.
"""

import os

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(PACKAGE_DIR)

__all__ = ["PROJECT_ROOT", "PACKAGE_DIR"]
