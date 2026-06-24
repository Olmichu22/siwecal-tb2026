"""
Recompute per-event metrics from cached raw per-hit arrays after applying a
hit_energy threshold.  Mirrors siwecal_validation.event_data so that at
threshold = 0.0 the results match the pre-computed valcache columns.

Called by ViewerController when the interactive hit-energy slider is non-zero.
Results are cached server-side per (path, threshold) so each level is computed
only once per session.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from siwecal_validation import metrics

# Mirror siwecal_validation.config defaults so threshold=0 matches valcache.
_SHOWER_THR   = 5.0   # per-layer nhit threshold (profile_kind = "nhit")
_SHOWER_PEAK  = 10.0  # minimum peak nhit for shower classification
_START_FRAC   = 0.1   # fraction of peak for start/end_layer_10
_MOLIERE_CONT = 0.90  # 90% transverse containment
_W_X0_MM      = 3.5   # tungsten radiation length [mm]


def recompute_all_metrics(
    reader,
    w_thickness_mm: np.ndarray,
    threshold: float,
    n_layers: int,
) -> pd.DataFrame:
    """Per-event DataFrame recomputed from hits surviving hit_energy >= threshold.

    Parameters
    ----------
    reader : EventFileReader
        Reader whose ``all_hits()`` is already cached.
    w_thickness_mm : np.ndarray
        Per-slab tungsten thickness [mm] (from ``DetectorModel.w_thickness_mm``).
    threshold : float
        Minimum hit energy [MIP] to keep.
    n_layers : int
        Number of detector layers.
    """
    all_hits  = reader.all_hits()
    slab_col  = all_hits["hit_slab"]
    en_col    = all_hits["hit_energy"]
    x_col     = all_hits.get("hit_x")
    y_col     = all_hits.get("hit_y")
    has_xy    = (x_col is not None) and (y_col is not None)

    # Pad w_over_x0 to at least n_layers so slab indexing never goes out of bounds.
    w_over_x0 = np.zeros(n_layers, dtype=float)
    src = np.asarray(w_thickness_mm, dtype=float) / _W_X0_MM
    w_over_x0[:min(len(src), n_layers)] = src[:n_layers]

    rows = []
    for i in range(len(slab_col)):
        slab   = np.asarray(slab_col[i], dtype=np.int64)
        energy = np.asarray(en_col[i],   dtype=float)

        mask   = energy >= threshold
        slab   = slab[mask]
        energy = energy[mask]
        x = np.asarray(x_col[i], dtype=float)[mask] if has_xy else None
        y = np.asarray(y_col[i], dtype=float)[mask] if has_xy else None

        rows.append(_metrics_one(slab, energy, x, y, w_over_x0, n_layers))

    return pd.DataFrame(rows)


def _metrics_one(slab, energy, x, y, w_over_x0, n_layers):
    NAN = float("nan")
    if slab.size == 0:
        return {
            "nhit": 0, "sum_energy": 0.0, "n_layers_hit": 0,
            "first_layer": NAN, "last_layer": NAN, "zbary": NAN,
            "mip_likeness": 0.0, "e_over_nhit": NAN, "weighte": 0.0,
            "bar_x": NAN, "bar_y": NAN, "bar_r": NAN,
            "transverse_rms": NAN, "moliere": NAN,
            "is_shower": False,
            "shower_start": NAN, "shower_max": NAN, "shower_end": NAN,
            "shower_start_10": NAN, "shower_end_10": NAN, "shower_length": NAN,
        }

    # Clip out-of-range slab indices before any indexing into w_over_x0.
    valid_slab = (slab >= 0) & (slab < n_layers)
    if not valid_slab.all():
        slab   = slab[valid_slab]
        energy = energy[valid_slab]
        if x is not None:
            x = x[valid_slab]
            y = y[valid_slab]

    nhit       = slab.size
    sum_energy = metrics.total_energy(energy)
    hits_layer = metrics.hits_per_layer(slab, n_layers)
    first, last, n_layers_hit = metrics.layer_extent(hits_layer)
    zbary_val  = metrics.zbary(slab, energy)
    mip_score  = metrics.mip_likeness(hits_layer, n_layers)
    weighte    = metrics.weighte_total(slab, energy, w_over_x0)

    # Transverse moments — use only energy > 0 hits (valcache convention).
    bar_x = bar_y = bar_r = transverse_rms_val = moliere_val = NAN
    if x is not None and energy.size:
        pos = energy > 0
        if pos.any():
            pw = metrics.hit_weights(slab[pos], energy[pos], w_over_x0)
            bar_x, bar_y, bar_r = metrics.barycenter_xy(x[pos], y[pos], pw)
            if np.isfinite(bar_x):
                transverse_rms_val = metrics.transverse_rms(
                    x[pos], y[pos], pw, bar_x, bar_y)

    # Longitudinal shower shape (profile_kind = "nhit" matches valcache default).
    sh = metrics.shower_features(hits_layer, _SHOWER_THR, _SHOWER_PEAK, _START_FRAC)

    if sh.is_shower and x is not None and np.isfinite(bar_x):
        pos = energy > 0
        if pos.any():
            pw = metrics.hit_weights(slab[pos], energy[pos], w_over_x0)
            moliere_val = metrics.moliere_radius(
                x[pos], y[pos], pw, bar_x, bar_y, _MOLIERE_CONT)

    return {
        "nhit": nhit,
        "sum_energy": sum_energy,
        "n_layers_hit": n_layers_hit,
        "first_layer": first,
        "last_layer": last,
        "zbary": zbary_val,
        "mip_likeness": mip_score,
        "e_over_nhit": (sum_energy / nhit) if nhit else NAN,
        "weighte": weighte,
        "bar_x": bar_x, "bar_y": bar_y, "bar_r": bar_r,
        "transverse_rms": transverse_rms_val,
        "moliere": moliere_val,
        "is_shower": sh.is_shower,
        "shower_start": sh.start_layer,
        "shower_max": sh.max_layer,
        "shower_end": sh.end_layer,
        "shower_start_10": sh.start_layer_10,
        "shower_end_10": sh.end_layer_10,
        "shower_length": sh.length,
    }
