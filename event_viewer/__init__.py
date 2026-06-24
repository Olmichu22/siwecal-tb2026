"""
Interactive event viewer for the SiW-ECAL TB2026CERN prototype.

A small OOP package (read with ``uproot``, drawn with Plotly Dash) to inspect
reconstructed ``ecal`` events one by one, explore file-level distributions with
dynamic cuts, and run simple clustering on per-event variables.

The package reuses the existing geometry / pad map / metrics code that lives at
the project root, so it must be importable. To make ``import event_viewer`` work
regardless of the current working directory we prepend the project root (the
parent of this package) to ``sys.path`` here.
"""

import os
import sys

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(PACKAGE_DIR)

if PROJECT_ROOT not in sys.path:
    sys.path.insert(0, PROJECT_ROOT)

__all__ = ["PROJECT_ROOT", "PACKAGE_DIR"]
