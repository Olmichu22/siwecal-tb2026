"""
Shared infrastructure for the SiW-ECAL TB2026CERN monorepo.

This tiny package holds code that the three analysis packages
(:mod:`siwecal_eventbuilder`, :mod:`siwecal_validation`, :mod:`event_viewer`)
all need but none of them should own. For now that is just :mod:`paths`, the
single source of truth for *where things live* (data, calibration, geometry,
configs, output), driven by the repo-root ``settings.yml``.

Keeping path resolution here means no module hard-codes an absolute ``/eos``
path: change ``settings.yml`` and every package follows.
"""

from . import paths

__all__ = ["paths"]
