"""
2-D per-layer view: a 4x4 grid of the 15 detector layers.

Each subplot shows that layer's pads as a faint grey grid (all channels of the
mapping) with the event's hit pads drawn on top, coloured by energy. Mirrors the
existing ``debug_pad_reflection_4x4.py`` diagnostic, but interactive.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

N_COLS = 4


class LayerGrid2D:
    """Builds the 4x4 grid of per-layer pad maps for one event."""

    def __init__(self, detector, colorscale: str = "Viridis"):
        self.detector = detector
        self.colorscale = colorscale
        self.n_layers = detector.geometry.n_slab_positions
        self.n_rows = int(np.ceil(self.n_layers / N_COLS))

    def build(self, event, color_clip: bool = True) -> go.Figure:
        fig = make_subplots(
            rows=self.n_rows, cols=N_COLS,
            subplot_titles=[f"layer {s}" for s in range(self.n_layers)],
            horizontal_spacing=0.03, vertical_spacing=0.06,
        )

        energy = event.energy if event is not None else np.empty(0)
        if color_clip:
            cmax = float(np.nanmax(energy)) if energy.size else 1.0
            cmin = 0.0
        else:
            bound = float(np.nanmax(np.abs(energy))) if energy.size else 1.0
            cmin, cmax = -bound, bound
        if cmax <= cmin:
            cmax = cmin + 1.0

        # Per-layer hit lookup.
        hit_by_layer = {}
        if event is not None and event.n_hits:
            for s in np.unique(event.slab):
                mask = event.slab == s
                hit_by_layer[int(s)] = (event.x[mask], event.y[mask],
                                        event.energy[mask])

        for slab in range(self.n_layers):
            row, col = divmod(slab, N_COLS)
            row += 1
            col += 1
            pads = self.detector.pads_for_slab(slab)
            if pads.size:
                fig.add_trace(go.Scattergl(
                    x=pads[:, 0], y=pads[:, 1], mode="markers",
                    marker=dict(size=3, color="rgba(180,180,180,0.5)",
                                symbol="square"),
                    hoverinfo="skip", showlegend=False), row=row, col=col)
            if slab in hit_by_layer:
                hx, hy, he = hit_by_layer[slab]
                fig.add_trace(go.Scattergl(
                    x=hx, y=hy, mode="markers",
                    marker=dict(size=9, color=he, colorscale=self.colorscale,
                                cmin=cmin, cmax=cmax, symbol="square",
                                showscale=(slab == 0),
                                colorbar=dict(title="E [MIP]", x=1.02)),
                    text=[f"E = {e:.2f}" for e in he], hoverinfo="text",
                    showlegend=False), row=row, col=col)
            fig.update_xaxes(scaleanchor="y", scaleratio=1, row=row, col=col,
                             showticklabels=False)
            fig.update_yaxes(showticklabels=False, row=row, col=col)

        fig.update_layout(margin=dict(l=0, r=40, t=20, b=0),
                          uirevision="layers2d")
        return fig
