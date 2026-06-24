"""
Cheap, derivable per-event quantities for files that lack the validation
metrics cache.

Only quantities computable directly from the per-hit ``(slab, energy)`` arrays
are produced here -- the ones that are inexpensive over the whole file. The heavy
discrimination metrics (Molière radius, shower shape, …) are intentionally *not*
recomputed; for those the user should run the validation pipeline that writes a
``.valcache.root``.

The formulas are reused verbatim from :mod:`siwecal_validation.metrics`, the
single source of truth for event-level physics quantities.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from siwecal_validation import metrics


# Columns produced here (order defines the DataFrame layout).
CHEAP_COLUMNS = (
    "nhit", "sum_energy", "n_layers_hit", "first_layer", "last_layer",
    "zbary", "mip_likeness", "e_over_nhit",
)


def compute_cheap_metrics(reader) -> pd.DataFrame:
    """Compute the cheap per-event quantities for every event in ``reader``.

    Parameters
    ----------
    reader : event_viewer.io.EventFileReader
        Reader whose tree exposes ``hit_slab`` and ``hit_energy`` per-hit branches.

    Returns
    -------
    pandas.DataFrame
        One row per event, columns :data:`CHEAP_COLUMNS`. Rows are aligned with
        the tree entry order (no events are skipped, so the index matches
        ``reader.read_hits``).
    """
    n_layers = reader.n_layers
    arrays = reader.tree.arrays(["hit_slab", "hit_energy"], library="np")
    slab_col, energy_col = arrays["hit_slab"], arrays["hit_energy"]

    rows = []
    for slab_raw, energy_raw in zip(slab_col, energy_col):
        slab = np.asarray(slab_raw, dtype=np.int64)
        energy = np.asarray(energy_raw, dtype=float)
        nhit = slab.size
        if nhit == 0:
            rows.append((0, 0.0, 0, np.nan, np.nan, np.nan, 0.0, np.nan))
            continue
        hits_layer = metrics.hits_per_layer(slab, n_layers)
        first, last, n_layers_hit = metrics.layer_extent(hits_layer)
        sum_energy = metrics.total_energy(energy)
        rows.append((
            nhit,
            sum_energy,
            n_layers_hit,
            first,
            last,
            metrics.zbary(slab, energy),
            metrics.mip_likeness(hits_layer, n_layers),
            sum_energy / nhit,
        ))

    return pd.DataFrame(rows, columns=list(CHEAP_COLUMNS))
