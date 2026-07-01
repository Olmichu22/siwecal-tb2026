"""
3-D detector scene: silicon planes + tungsten absorbers, with event hits overlaid.

``DetectorScene3D`` builds the static geometry once (a translucent quad per Si
plane, a translucent box per W plate) and caches those traces. ``event_figure``
clones the cached geometry and adds a ``Scatter3d`` of the event's hit pads,
coloured by hit energy.

Note on negative energies: ``hit_energy`` can be slightly negative when the raw
ADC sits below pedestal (a known artefact of the dummy calibration, not physics).
By default the colour range is clipped at 0 so the colour scale tracks real
signal; ``color_clip=False`` switches to the actual ``[min, max]`` of the shown
hits, so the scale reflects the real (incl. negative) energies and its lower end
maps to the faintest displayed pad.
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
                     threshold=None, show_moliere: bool = False,
                     show_axis: bool = False,
                     axis_mode: str = "weighted") -> go.Figure:
        """Detector geometry + this event's hits coloured by energy.

        ``threshold`` (energy, MIP) fades hits *below* it: they are drawn in a
        separate, almost-transparent mesh so the high-energy core -- the shower
        profile -- stands out. ``None`` or ``0`` shows every hit fully opaque.

        ``show_moliere`` overlays the Molière cylinder from the per-event
        ``bar_x``/``bar_y``/``moliere`` metrics; ``show_axis`` overlays the shower
        axis built from the per-layer energy-weighted (``axis_mode="weighted"``)
        or geometric (``"geom"``) barycenters of this event's hits.
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
            if show_moliere:
                traces += self._moliere_traces(event)
            if show_axis:
                axis = self._axis_trace(event, axis_mode)
                if axis is not None:
                    traces.append(axis)
        return self._wrap(traces)

    # ---------------------------------------------------------- shower overlays --
    def _shower_z_range(self, event):
        """(z0, z1) mm spanning the shower, from metrics or the hit z extent."""
        m = event.metrics or {}
        n = self.detector.slab_z_mm.size

        def _z_of_layer(key):
            v = m.get(key)
            if v is None or not np.isfinite(v):
                return None
            i = int(round(float(v)))
            return float(self.detector.slab_z_mm[i]) if 0 <= i < n else None

        z0 = _z_of_layer("shower_start")
        z1 = _z_of_layer("shower_end")
        if z0 is None or z1 is None:
            z0 = z0 if z0 is not None else _z_of_layer("first_layer")
            z1 = z1 if z1 is not None else _z_of_layer("last_layer")
        if z0 is None or z1 is None:
            z = np.asarray(event.z, dtype=float)
            z = z[np.isfinite(z)]
            if z.size == 0:
                return None
            z0, z1 = float(z.min()), float(z.max())
        return (min(z0, z1), max(z0, z1))

    def _moliere_traces(self, event, n_theta: int = 48):
        """Translucent lateral cylinder + two ring outlines at the Molière radius.

        Centred on the transverse barycenter ``(bar_x, bar_y)`` with radius
        ``moliere``, spanning the shower's z range. Returns an empty list when the
        needed metrics are missing (e.g. non-valcache files, or non-showers).
        """
        m = event.metrics or {}
        bx, by, r = m.get("bar_x"), m.get("bar_y"), m.get("moliere")
        if bx is None or by is None or r is None:
            return []
        if not (np.isfinite(bx) and np.isfinite(by) and np.isfinite(r)) or r <= 0:
            return []
        zr = self._shower_z_range(event)
        if zr is None:
            return []
        z0, z1 = zr
        theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=True)
        cx, cy = bx + r * np.cos(theta), by + r * np.sin(theta)
        color = "rgba(224, 87, 46, 0.9)"

        # Lateral surface: two rings (z0, z1) stitched with triangles.
        xs = np.concatenate([cx, cx])
        ys = np.concatenate([cy, cy])
        zs = np.concatenate([np.full(n_theta, z0), np.full(n_theta, z1)])
        tri_i, tri_j, tri_k = [], [], []
        for a in range(n_theta - 1):
            b, lo0, lo1 = a + 1, a, a + 1
            up0, up1 = a + n_theta, a + 1 + n_theta
            tri_i += [lo0, lo1]
            tri_j += [up0, up0]
            tri_k += [up1, lo1]
        surface = go.Mesh3d(
            x=xs, y=ys, z=zs, i=tri_i, j=tri_j, k=tri_k,
            color="rgba(224, 87, 46, 0.12)", opacity=0.12, hoverinfo="skip",
            showscale=False, flatshading=True, name="Molière cylinder")
        rings = go.Scatter3d(
            x=np.concatenate([cx, [np.nan], cx]),
            y=np.concatenate([cy, [np.nan], cy]),
            z=np.concatenate([np.full(n_theta, z0), [np.nan],
                              np.full(n_theta, z1)]),
            mode="lines", line=dict(color=color, width=4),
            name=f"Molière r={r:.1f} mm", hoverinfo="name")
        return [surface, rings]

    def _axis_trace(self, event, axis_mode: str = "weighted"):
        """Polyline through the per-layer barycenters of the event's hits.

        ``axis_mode="weighted"`` uses the energy-weighted transverse barycenter
        of each layer (positive energies only); ``"geom"`` uses the plain mean of
        the hit positions. Returns ``None`` for an empty event.
        """
        x = np.asarray(event.x, dtype=float)
        y = np.asarray(event.y, dtype=float)
        z = np.asarray(event.z, dtype=float)
        e = np.asarray(event.energy, dtype=float)
        good = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
        if not good.any():
            return None
        x, y, z, e = x[good], y[good], z[good], e[good]
        weighted = axis_mode != "geom"
        bx, by, bz = [], [], []
        for zl in np.unique(z):
            sel = z == zl
            w = e[sel]
            if weighted and np.any(w > 0):
                w = np.where(w > 0, w, 0.0)
                tot = w.sum()
                bx.append(float(np.dot(w, x[sel]) / tot))
                by.append(float(np.dot(w, y[sel]) / tot))
            else:
                bx.append(float(x[sel].mean()))
                by.append(float(y[sel].mean()))
            bz.append(float(zl))
        name = "axis (energy-weighted)" if weighted else "axis (geometric)"
        return go.Scatter3d(
            x=bx, y=by, z=bz, mode="lines+markers",
            line=dict(color="#111111", width=5),
            marker=dict(size=3, color="#111111"),
            name=name, hoverinfo="name")

    @staticmethod
    def _color_range(energy, color_clip):
        if color_clip:
            cmin, cmax = 0.0, float(np.nanmax(energy)) if energy.size else 1.0
        else:
            cmin = float(np.nanmin(energy)) if energy.size else 0.0
            cmax = float(np.nanmax(energy)) if energy.size else 1.0
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
