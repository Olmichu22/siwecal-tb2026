"""
2-D per-layer view: a 4x4 grid of the 15 detector layers.

Each subplot shows that layer's pads as a faint grey grid (all channels of the
mapping), outlined sensor by sensor, with the event's hit pads drawn on top,
coloured by energy. The outlines mark the silicon sensors: the cross between
them is dead -- the inactive rim of each sensor, where no pad exists and no hit
can appear. Mirrors the existing ``debug_pad_reflection_4x4.py`` diagnostic, but
interactive.

Pads are drawn as ``Heatmap`` cells, not as scatter markers. A marker is sized in
pixels, so it matches the pad at exactly one zoom level and one window
size: in this 4x4 grid a pad is about 4 px across while the hit markers were 9,
which made neighbouring hits overlap and spill over the guard ring they are
supposed to stop at. Cells are sized in data coordinates and stay honest at any
zoom. See :meth:`DetectorModel.pad_cell_grid` for the spacer cell that keeps the
dead region from being absorbed into its neighbours.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

N_COLS = 4

# Flat colourscale for the inactive pad grid: one grey, whatever the z value.
PAD_GREY = [[0.0, "rgba(200,200,200,0.45)"], [1.0, "rgba(200,200,200,0.45)"]]


def _cell_matrix(grid, xs, ys, values):
    """Place ``values`` on the layer's cell grid; empty cells stay ``NaN``.

    ``NaN`` is what keeps a cell unpainted, so the guard-ring spacer and every
    pad without a hit render as background rather than as a zero-energy cell.
    """
    x_edges, y_edges, x_index, y_index = grid
    z = np.full((y_edges.size - 1, x_edges.size - 1), np.nan)
    for x, y, value in zip(xs, ys, values):
        col = x_index.get(round(float(x), 3))
        row = y_index.get(round(float(y), 3))
        if col is not None and row is not None:
            z[row, col] = value
    return z


def _wafer_outline(rects):
    """Closed rectangle outlines for every wafer, as one NaN-separated trace."""
    xs, ys = [], []
    for x0, x1, y0, y1 in rects:
        xs += [x0, x1, x1, x0, x0, np.nan]
        ys += [y0, y0, y1, y1, y0, np.nan]
    return go.Scattergl(
        x=xs, y=ys, mode="lines",
        line=dict(color="rgba(140,140,140,0.55)", width=1),
        hoverinfo="skip", showlegend=False)


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
            cmin = float(np.nanmin(energy)) if energy.size else 0.0
            cmax = float(np.nanmax(energy)) if energy.size else 1.0
        if cmax <= cmin:
            cmax = cmin + 1.0

        # Per-layer hit lookup.
        hit_by_layer = {}
        if event is not None and event.n_hits:
            for s in np.unique(event.slab):
                mask = event.slab == s
                hit_by_layer[int(s)] = (event.x[mask], event.y[mask],
                                        event.energy[mask])

        wafers = self.detector.wafer_rects
        for slab in range(self.n_layers):
            row, col = divmod(slab, N_COLS)
            row += 1
            col += 1
            pads = self.detector.pads_for_slab(slab)
            grid = self.detector.pad_cell_grid(slab)
            x_edges, y_edges = grid[0], grid[1]

            if pads.size:
                fig.add_trace(go.Heatmap(
                    x=x_edges, y=y_edges,
                    z=_cell_matrix(grid, pads[:, 0], pads[:, 1],
                                   np.zeros(len(pads))),
                    colorscale=PAD_GREY, zmin=0.0, zmax=1.0, showscale=False,
                    hoverinfo="skip"), row=row, col=col)
            if slab in hit_by_layer:
                hx, hy, he = hit_by_layer[slab]
                fig.add_trace(go.Heatmap(
                    x=x_edges, y=y_edges, z=_cell_matrix(grid, hx, hy, he),
                    colorscale=self.colorscale, zmin=cmin, zmax=cmax,
                    showscale=(slab == 0),
                    colorbar=dict(title="E [MIP]", x=1.02),
                    hoverongaps=False,
                    hovertemplate="E = %{z:.2f} MIP<extra></extra>"),
                    row=row, col=col)
            # Drawn last so the outlines stay on top of the filled cells.
            if wafers:
                fig.add_trace(_wafer_outline(wafers), row=row, col=col)
            fig.update_xaxes(scaleanchor="y", scaleratio=1, row=row, col=col,
                             showticklabels=False)
            fig.update_yaxes(showticklabels=False, row=row, col=col)

        fig.update_layout(margin=dict(l=0, r=40, t=20, b=0),
                          uirevision="layers2d")
        return fig
