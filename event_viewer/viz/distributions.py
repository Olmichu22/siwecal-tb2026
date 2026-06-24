"""
File-level distribution plots: histograms (with optional cut shading) and
clustering scatter plots.

These operate on the (already cut-filtered) per-event DataFrame, so the same
cuts that drive the event navigation also shape the distributions.
"""

from __future__ import annotations

from typing import List, Optional

import numpy as np
import plotly.express as px
import plotly.graph_objects as go

from ..analysis.clustering import UNCLUSTERED


class DistributionPlots:
    """Builds histogram and scatter figures for event-level variables."""

    def histogram(self, values: np.ndarray, variable: str,
                  cut_range: Optional[tuple] = None, nbins: int = 60,
                  labels: Optional[np.ndarray] = None) -> go.Figure:
        """1-D histogram of ``values`` with ``nbins`` bins.

        If ``labels`` is given (one cluster label per value), the histogram is
        split into one stacked series per cluster (``barmode="stack"``) so the
        cluster separation is visible; otherwise a single series is drawn. The
        bin edges are shared across series so the stacking aligns. ``cut_range``
        is shaded when provided.
        """
        values = np.asarray(values, dtype=float)
        finite = values[np.isfinite(values)]
        fig = go.Figure()
        if finite.size == 0:
            fig.update_layout(xaxis_title=variable, yaxis_title="events",
                              margin=dict(l=50, r=20, t=20, b=40))
            return fig

        start, end = float(finite.min()), float(finite.max())
        if end <= start:
            end = start + 1.0
        size = (end - start) / max(int(nbins), 1)
        xbins = dict(start=start, end=end + size, size=size)

        if labels is None:
            fig.add_trace(go.Histogram(x=finite, xbins=xbins, name=variable,
                                       marker_color="#4c78a8"))
        else:
            labels = np.asarray(labels)
            for lab in np.unique(labels):
                mask = labels == lab
                series = values[mask]
                series = series[np.isfinite(series)]
                name = "unclustered" if lab == UNCLUSTERED else f"cluster {lab}"
                marker = dict(color="lightgrey") if lab == UNCLUSTERED else None
                fig.add_trace(go.Histogram(
                    x=series, xbins=xbins, name=name, bingroup="dist",
                    opacity=0.85, marker=marker))
            fig.update_layout(barmode="stack")

        if cut_range is not None:
            lo, hi = cut_range
            fig.add_vrect(x0=lo, x1=hi, fillcolor="orange", opacity=0.15,
                          line_width=0)
        fig.update_layout(
            xaxis_title=variable, yaxis_title="events",
            margin=dict(l=50, r=20, t=20, b=40), bargap=0.02,
            uirevision=f"hist-{variable}")
        return fig

    def scatter(self, df, xvar: str, yvar: str,
                labels: Optional[List[int]] = None) -> go.Figure:
        """2-D scatter of ``xvar`` vs ``yvar``, coloured by cluster label."""
        x = df[xvar].to_numpy(dtype=float)
        y = df[yvar].to_numpy(dtype=float)
        if labels is None:
            fig = go.Figure(go.Scattergl(
                x=x, y=y, mode="markers",
                marker=dict(size=5, color="#4c78a8", opacity=0.6)))
        else:
            labels = np.asarray(labels)
            fig = go.Figure()
            for lab in np.unique(labels):
                mask = labels == lab
                name = "unclustered" if lab == UNCLUSTERED else f"cluster {lab}"
                color = "lightgrey" if lab == UNCLUSTERED else None
                fig.add_trace(go.Scattergl(
                    x=x[mask], y=y[mask], mode="markers", name=name,
                    marker=dict(size=5, opacity=0.7,
                                color=color)))
        fig.update_layout(
            xaxis_title=xvar, yaxis_title=yvar,
            margin=dict(l=50, r=20, t=20, b=40),
            uirevision=f"scatter-{xvar}-{yvar}")
        return fig
