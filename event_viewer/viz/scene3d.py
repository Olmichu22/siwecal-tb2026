"""
3-D detector scene: silicon planes + tungsten absorbers, with event hits overlaid.

``DetectorScene3D`` builds the static geometry once (a translucent quad per Si
plane, a translucent box per W plate) and caches those traces. ``event_figure``
clones the cached geometry and adds a ``Scatter3d`` of the event's hit pads,
coloured by hit energy.

Note on negative energies: ``hit_energy`` can be slightly negative when the raw
ADC sits below pedestal (a known artefact of the dummy calibration, not physics).
By default the colour range is clipped at 0 so the colour scale tracks real
signal; ``color_clip=False`` switches to a symmetric range to expose the noise.
"""

from __future__ import annotations

from typing import List

import numpy as np
import plotly.graph_objects as go

SILICON_COLOR = "rgba(120, 170, 220, 0.12)"
TUNGSTEN_COLOR = "rgba(150, 150, 150, 0.04)"   # very transparent absorber slabs


def _quad_mesh(x0, x1, y0, y1, z, color, name):
    """A flat rectangle in the z=const plane as a 2-triangle Mesh3d."""
    return go.Mesh3d(
        x=[x0, x1, x1, x0], y=[y0, y0, y1, y1], z=[z, z, z, z],
        i=[0, 0], j=[1, 2], k=[2, 3],
        color=color, opacity=1.0, hoverinfo="skip", name=name,
        showscale=False, flatshading=True,
    )


def _box_mesh(x0, x1, y0, y1, z0, z1, color, name):
    """An axis-aligned box as a 12-triangle Mesh3d."""
    xs = [x0, x1, x1, x0, x0, x1, x1, x0]
    ys = [y0, y0, y1, y1, y0, y0, y1, y1]
    zs = [z0, z0, z0, z0, z1, z1, z1, z1]
    i = [0, 0, 0, 0, 4, 4, 0, 1, 1, 2, 3, 0]
    j = [1, 2, 4, 3, 5, 6, 1, 5, 2, 6, 7, 4]
    k = [2, 3, 5, 7, 6, 7, 5, 6, 6, 7, 4, 7]
    return go.Mesh3d(
        x=xs, y=ys, z=zs, i=i, j=j, k=k,
        color=color, opacity=1.0, hoverinfo="skip", name=name,
        showscale=False, flatshading=True,
    )


class DetectorScene3D:
    """Builds the 3-D detector scene and overlays event hits."""

    def __init__(self, detector, colorscale: str = "Viridis"):
        self.detector = detector
        self.colorscale = colorscale
        self._base_traces = self._build_base()

    def _build_base(self) -> List[go.Mesh3d]:
        (x0, x1), (y0, y1) = self.detector.x_extent, self.detector.y_extent
        traces: List[go.Mesh3d] = []
        for slab, z in self.detector.silicon_quads():
            traces.append(_quad_mesh(x0, x1, y0, y1, z,
                                     SILICON_COLOR, f"Si slab {slab}"))
        for slab, z_down, z_up in self.detector.tungsten_boxes():
            traces.append(_box_mesh(x0, x1, y0, y1, z_up, z_down,
                                    TUNGSTEN_COLOR, f"W slab {slab}"))
        return traces

    # ------------------------------------------------------------- figures --
    def base_figure(self) -> go.Figure:
        """The detector geometry alone (no event)."""
        return self._wrap(list(self._base_traces))

    def event_figure(self, event, color_clip: bool = True,
                     threshold=None) -> go.Figure:
        """Detector geometry + this event's hits coloured by energy.

        ``threshold`` (energy, MIP) fades hits *below* it: they are drawn in a
        separate, almost-transparent mesh so the high-energy core -- the shower
        profile -- stands out. ``None`` or ``0`` shows every hit fully opaque.
        """
        traces = list(self._base_traces)
        if event is not None and event.n_hits:
            energy = np.asarray(event.energy, dtype=float)
            cmin, cmax = self._color_range(energy, color_clip)
            if threshold:
                above = energy >= float(threshold)
            else:
                above = np.ones(energy.shape, dtype=bool)
            below = ~above
            any_above = bool(above.any())
            if below.any():
                faint = self._square_mesh(event, below, cmin, cmax,
                                          opacity=0.06, showscale=not any_above)
                if faint is not None:
                    traces.append(faint)
            if any_above:
                strong = self._square_mesh(event, above, cmin, cmax,
                                           opacity=1.0, showscale=True)
                if strong is not None:
                    traces.append(strong)
        return self._wrap(traces)

    @staticmethod
    def _color_range(energy, color_clip):
        if color_clip:
            cmin, cmax = 0.0, float(np.nanmax(energy)) if energy.size else 1.0
        else:
            bound = float(np.nanmax(np.abs(energy))) if energy.size else 1.0
            cmin, cmax = -bound, bound
        if cmax <= cmin:
            cmax = cmin + 1.0
        return cmin, cmax

    def _square_mesh(self, event, mask, cmin, cmax, opacity, showscale):
        """One ``Mesh3d`` of pad-sized squares for the hits selected by ``mask``.

        4 vertices + 2 triangles per hit; per-vertex ``intensity`` = hit energy so
        every mesh shares the same energy colour scale regardless of ``opacity``.
        """
        half = self.detector.pad_pitch / 2.0
        xs, ys, zs, intensity, hovertext = [], [], [], [], []
        tri_i, tri_j, tri_k = [], [], []
        for n in np.flatnonzero(mask):
            x, y, z, e = event.x[n], event.y[n], event.z[n], event.energy[n]
            if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)):
                continue
            base = len(xs)
            xs += [x - half, x + half, x + half, x - half]
            ys += [y - half, y - half, y + half, y + half]
            zs += [z, z, z, z]
            intensity += [e, e, e, e]
            label = (f"slab {event.slab[n]}, chip {event.chip[n]}, "
                     f"chan {event.chan[n]}<br>E = {e:.2f} MIP")
            hovertext += [label] * 4
            tri_i += [base, base]
            tri_j += [base + 1, base + 2]
            tri_k += [base + 2, base + 3]

        if not xs:
            return None
        return go.Mesh3d(
            x=xs, y=ys, z=zs, i=tri_i, j=tri_j, k=tri_k,
            intensity=intensity, colorscale=self.colorscale,
            cmin=cmin, cmax=cmax, showscale=showscale,
            colorbar=dict(title="E [MIP]", x=1.02) if showscale else None,
            opacity=opacity, flatshading=True, name="hits",
            hovertext=hovertext, hoverinfo="text",
        )

    def _wrap(self, traces) -> go.Figure:
        fig = go.Figure(data=traces)
        # Default orientation: x to the right, y up and z coming *towards* the
        # viewer (out of the screen). This is the natural right-handed view with
        # y as the "up" vector, so the z axis keeps its normal direction.
        fig.update_layout(
            scene=dict(
                xaxis=dict(title="x [mm]"),
                yaxis=dict(title="y [mm]"),
                zaxis=dict(title="z [mm]"),
                aspectmode="data",
                camera=dict(
                    up=dict(x=0, y=1, z=0),
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=1.6, y=1.1, z=1.6),
                ),
            ),
            margin=dict(l=0, r=0, t=0, b=0), showlegend=False,
            uirevision="detector",  # keep camera across event navigation
        )
        return fig
