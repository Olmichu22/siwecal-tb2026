"""
Dash callbacks wiring the stores, the controller and the figures together.

Session state lives in light ``dcc.Store``s:

* ``store-file``       : path of the loaded ROOT file.
* ``store-pos``        : absolute entry index of the shown event (stable across
                         threshold changes; position in passing list is derived).
* ``store-event-cuts`` : Event-tab cuts -- limit the one-by-one navigation and
                         drive that tab's bicolor (kept/removed) histogram.
* ``store-cuts``       : Clustering-tab cuts -- limit that tab's histogram and the
                         events entering the clustering fit. Independent from the
                         Event tab, so the two cut sets never interfere.
* ``store-cluster``    : ``{passing: [...], labels: [...]}`` from the last run.

Each cut store has a single writer: its pattern-matching slider callback (resetting
the matching cut-vars to ``[]`` clears it). ``store-pos`` has one primary writer
(navigation) plus a reset on file / event-cut change.
"""

from __future__ import annotations

import os

import numpy as np
import plotly.graph_objects as go
from dash import ALL, MATCH, Input, Output, State, ctx, dcc, html, no_update
from dash.exceptions import PreventUpdate

from ..analysis.cuts import CutModel


def _empty_fig(message: str = "") -> go.Figure:
    fig = go.Figure()
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    if message:
        fig.add_annotation(text=message, showarrow=False,
                           xref="paper", yref="paper", x=0.5, y=0.5)
    return fig


def _pos_of(passing, cur_index) -> int:
    """Position of absolute entry ``cur_index`` in the passing list, else 0."""
    if cur_index is None:
        return 0
    loc = np.where(np.asarray(passing) == cur_index)[0]
    return int(loc[0]) if loc.size else 0


def register_callbacks(app, controller) -> None:
    """Attach every callback to ``app``, all closing over ``controller``."""

    # --------------------------------------------------------- file loading --
    @app.callback(
        Output("store-file", "data"),
        Output("file-status", "children"),
        Output("dist-var", "options"), Output("dist-var", "value"),
        Output("cut-vars", "options"), Output("cut-vars", "value"),
        Output("ev-dist-var", "options"), Output("ev-dist-var", "value"),
        Output("ev-cut-vars", "options"), Output("ev-cut-vars", "value"),
        Output("cluster-features", "options"), Output("cluster-features", "value"),
        Output("scatter-x", "options"), Output("scatter-x", "value"),
        Output("scatter-y", "options"), Output("scatter-y", "value"),
        Output("store-cluster", "data", allow_duplicate=True),
        Output("hit-energy-slider", "disabled"),
        Output("hit-energy-slider", "value"),
        Input("load-btn", "n_clicks"),
        Input("file-dropdown", "value"),
        State("file-path", "value"),
        prevent_initial_call="initial_duplicate",
    )
    def load_file(_n_clicks, dropdown_value, text_value):
        triggered = ctx.triggered_id
        if triggered == "load-btn" and text_value and text_value.strip():
            path = text_value.strip()
        else:
            path = dropdown_value
        empty = []
        if not path:
            return (None, "Select a .root file", empty, None, empty, [],
                    empty, None, empty, [], empty, [],
                    empty, None, empty, None, None, True, 0)
        try:
            ds = controller.dataset(path)
        except Exception as error:  # noqa: BLE001 - surface any I/O failure
            return (no_update, f"Error opening: {error}", empty, no_update,
                    empty, no_update, empty, no_update, empty, no_update,
                    empty, no_update, empty, no_update, empty, no_update,
                    no_update, no_update, no_update)

        cols = ds.feature_columns()
        opts = [{"label": c, "value": c} for c in cols]
        dist_default = "energy" if "energy" in cols else (cols[0] if cols else None)
        x0 = cols[0] if cols else None
        y0 = cols[1] if len(cols) > 1 else x0
        flag = "with metrics" if ds.has_metrics \
            else "WITHOUT metrics (basic quantities only)"
        # Block the interactive MIP cut when it would fall back to the slow
        # in-memory recompute (no pre-computed branches) on a large file.
        max_events = controller.config.max_recompute_events
        block_cut = (not ds.reader.has_mip_thresholds) and ds.n_events > max_events
        if block_cut:
            flag += (f"; MIP cut disabled (>{max_events} events, no precomputed "
                     f"metrics — run validation to build a valcache)")
        status = f"{os.path.basename(path)} — {ds.n_events} events — {flag}"
        slider_value = 0 if block_cut else no_update
        return (path, status, opts, dist_default, opts, [],
                opts, dist_default, opts, [], opts, [],
                opts, x0, opts, y0, None, block_cut, slider_value)

    # -------------------------------------------------- hit-energy threshold --
    @app.callback(
        Output("store-hit-threshold", "data"),
        Input("hit-energy-slider", "value"),
    )
    def update_hit_threshold(value):
        # The "Computing metrics…" note is shown by the dcc.Loading spinners that
        # wrap the heavy graphs, so it appears only while the figures are being
        # (re)built and clears automatically when they finish.
        return float(value or 0.0)

    # ----------------------------------------------------- dynamic cut UI --
    @app.callback(
        Output("cut-sliders", "children"),
        Input("cut-vars", "value"),
        State("store-file", "data"),
        State("store-hit-threshold", "data"),
    )
    def build_sliders(cut_vars, path, hit_threshold):
        if not path or not cut_vars:
            return []
        thr = float(hit_threshold or 0.0)
        children = []
        for var in cut_vars:
            lo, hi = controller.variable_range(path, var, thr)
            step = (hi - lo) / 100 if hi > lo else 1.0
            children.append(html.Div(style={"marginBottom": "14px"}, children=[
                html.Label(var, style={"fontSize": "13px"}),
                dcc.RangeSlider(
                    id={"type": "cut-slider", "index": var},
                    min=lo, max=hi, value=[lo, hi], step=step, allowCross=False,
                    tooltip={"placement": "bottom", "always_visible": False}),
            ]))
        return children

    @app.callback(
        Output("store-cuts", "data"),
        Input({"type": "cut-slider", "index": ALL}, "value"),
        State({"type": "cut-slider", "index": ALL}, "id"),
    )
    def update_cuts(values, ids):
        cuts = []
        for value, ident in zip(values, ids):
            if value is None:
                continue
            cuts.append({"variable": ident["index"],
                         "lo": value[0], "hi": value[1]})
        return cuts

    # ------------------------------------ Event-tab dynamic cut UI (page 1) --
    # A second, independent copy of the cut widget. It writes ``store-event-cuts``
    # (which drives the one-by-one navigation and the bicolor histogram) and is
    # completely decoupled from the clustering tab's ``store-cuts``.
    @app.callback(
        Output("ev-cut-sliders", "children"),
        Input("ev-cut-vars", "value"),
        State("store-file", "data"),
        State("store-hit-threshold", "data"),
    )
    def build_sliders_event(cut_vars, path, hit_threshold):
        if not path or not cut_vars:
            return []
        thr = float(hit_threshold or 0.0)
        children = []
        for var in cut_vars:
            lo, hi = controller.variable_range(path, var, thr)
            step = (hi - lo) / 100 if hi > lo else 1.0
            children.append(html.Div(style={"marginBottom": "14px"}, children=[
                html.Label(var, style={"fontSize": "13px"}),
                dcc.RangeSlider(
                    id={"type": "ev-cut-slider", "index": var},
                    min=lo, max=hi, value=[lo, hi], step=step, allowCross=False,
                    tooltip={"placement": "bottom", "always_visible": False}),
            ]))
        return children

    @app.callback(
        Output("store-event-cuts", "data"),
        Input({"type": "ev-cut-slider", "index": ALL}, "value"),
        State({"type": "ev-cut-slider", "index": ALL}, "id"),
    )
    def update_cuts_event(values, ids):
        cuts = []
        for value, ident in zip(values, ids):
            if value is None:
                continue
            cuts.append({"variable": ident["index"],
                         "lo": value[0], "hi": value[1]})
        return cuts

    # ---------------------------------------------------------- navigation --
    @app.callback(
        Output("store-pos", "data"),
        Input("prev-btn", "n_clicks"),
        Input("next-btn", "n_clicks"),
        Input("event-input", "value"),
        State("store-pos", "data"),
        State("store-file", "data"),
        State("store-event-cuts", "data"),
        State("store-hit-threshold", "data"),
        prevent_initial_call=True,
    )
    def navigate(_prev, _next, event_input, cur_index, path, cuts, hit_threshold):
        if not path:
            raise PreventUpdate
        thr = float(hit_threshold or 0.0)
        passing = controller.passing_indices(path, CutModel.from_store(cuts), thr)
        n_pass = len(passing)
        if n_pass == 0:
            raise PreventUpdate
        pos = _pos_of(passing, cur_index)
        triggered = ctx.triggered_id
        if triggered == "prev-btn":
            pos -= 1
        elif triggered == "next-btn":
            pos += 1
        elif triggered == "event-input":
            pos = (int(event_input) - 1) if event_input else pos
        pos = max(0, min(n_pass - 1, pos))
        new_index = int(passing[pos])
        if new_index == cur_index:
            raise PreventUpdate
        return new_index

    @app.callback(
        Output("store-pos", "data", allow_duplicate=True),
        Input("store-file", "data"),
        Input("store-event-cuts", "data"),
        prevent_initial_call=True,
    )
    def reset_pos(_path, _cuts):
        return None

    # ------------------------------------------------------- event render --
    @app.callback(
        Output("scene3d", "figure"),
        Output("layers2d", "figure"),
        Output("metrics-table", "data"),
        Output("event-label", "children"),
        Output("event-input", "value"),
        Output("event-input", "max"),
        Input("store-file", "data"),
        Input("store-pos", "data"),
        Input("store-event-cuts", "data"),
        Input("color-clip", "value"),
        Input("store-hit-threshold", "data"),
        Input("show-overlays", "value"),
    )
    def render_event(path, cur_index, cuts, clip, hit_threshold, overlays):
        if not path:
            return (_empty_fig("Load a file"), _empty_fig(), [],
                    "no file", None, 1)
        thr = float(hit_threshold or 0.0)
        passing = controller.passing_indices(path, CutModel.from_store(cuts), thr)
        n_pass = len(passing)
        if n_pass == 0:
            return (_empty_fig("No events pass the cuts"), _empty_fig(),
                    [], "0 events passing cuts", None, 1)
        pos = _pos_of(passing, cur_index)
        index = int(passing[pos])
        color_clip = "clip" in (clip or [])
        overlays = overlays or []
        scene, layers, rows = controller.event_figures(
            path, index, color_clip, thr,
            show_moliere="moliere" in overlays, show_axis="axis" in overlays)
        ds = controller.dataset(path)
        label = (f"event {pos + 1} / {n_pass} passing "
                 f"(entry {index}, {ds.n_events} total)")
        return scene, layers, rows, label, pos + 1, n_pass

    # ------------------------------------------------------- distributions --
    @app.callback(
        Output("dist-hist", "figure"),
        Input("dist-var", "value"),
        Input("store-cuts", "data"),
        Input("store-file", "data"),
        Input("dist-nbins", "value"),
        Input("store-cluster", "data"),
        Input("dist-stack", "value"),
        Input("store-hit-threshold", "data"),
    )
    def update_histogram(variable, cuts, path, nbins, cluster, stack, hit_threshold):
        if not path or not variable:
            return _empty_fig("Select a variable")
        thr = float(hit_threshold or 0.0)
        use_cluster = cluster if (cluster and "stack" in (stack or [])) else None
        return controller.histogram(path, variable, CutModel.from_store(cuts),
                                    int(nbins or 60), use_cluster, thr)

    # -------------------------------------------- Event-tab distribution (bicolor)
    @app.callback(
        Output("ev-dist-hist", "figure"),
        Input("ev-dist-var", "value"),
        Input("store-event-cuts", "data"),
        Input("store-file", "data"),
        Input("ev-dist-nbins", "value"),
        Input("store-hit-threshold", "data"),
    )
    def update_event_histogram(variable, cuts, path, nbins, hit_threshold):
        if not path or not variable:
            return _empty_fig("Select a variable")
        thr = float(hit_threshold or 0.0)
        return controller.histogram_split(path, variable,
                                          CutModel.from_store(cuts),
                                          int(nbins or 60), thr)

    # ----------------------------------------------------------- clustering --
    @app.callback(
        Output("cluster-nclusters", "disabled"),
        Output("cluster-eps", "disabled"),
        Output("cluster-minsamples", "disabled"),
        Output("cluster-param-hint", "children"),
        Input("cluster-algo", "value"),
    )
    def toggle_cluster_params(algo):
        """Enable only the parameters each algorithm actually uses + a hint."""
        uses_k = algo in ("kmeans", "gmm", "spectral")
        uses_eps = algo == "dbscan"
        hints = {
            "kmeans": "K-Means uses n_clusters.",
            "gmm": "Gaussian Mixture uses n_clusters (n_components).",
            "spectral": "Spectral uses n_clusters.",
            "dbscan": "DBSCAN uses eps and min_samples (eps is in standardized "
                      "units, since features are z-scored).",
        }
        return (not uses_k, not uses_eps, not uses_eps, hints.get(algo, ""))

    @app.callback(
        Output("store-cluster", "data", allow_duplicate=True),
        Input("store-cuts", "data"),
        Input("store-hit-threshold", "data"),
        prevent_initial_call=True,
    )
    def invalidate_cluster_on_cut(_cuts, _thr):
        """A clustering run is tied to the cut it was computed under; when the cuts
        change its passing-index snapshot is stale, so drop it. This unsticks the
        histogram / scatter / cluster examples (which read the snapshot) and lets
        them follow the current selection again -- re-run clustering to refresh."""
        return None

    @app.callback(
        Output("store-cluster", "data", allow_duplicate=True),
        Input("cluster-reset", "n_clicks"),
        prevent_initial_call=True,
    )
    def reset_cluster(_n_clicks):
        """Clear the last clustering run so the cluster panels, scatter colouring
        and stacked histogram fall back to the plain (uncluster) view."""
        return None

    @app.callback(
        Output("store-cluster", "data"),
        Input("cluster-run", "n_clicks"),
        State("store-file", "data"),
        State("store-cuts", "data"),
        State("cluster-features", "value"),
        State("cluster-algo", "value"),
        State("cluster-nclusters", "value"),
        State("cluster-eps", "value"),
        State("cluster-minsamples", "value"),
        State("store-hit-threshold", "data"),
        prevent_initial_call=True,
    )
    def run_clustering(n_clicks, path, cuts, features, algo, n_clusters,
                       eps, min_samples, hit_threshold):
        if not path or not features:
            raise PreventUpdate
        thr = float(hit_threshold or 0.0)
        passing, labels = controller.run_clustering(
            path, CutModel.from_store(cuts), features, algo,
            int(n_clusters or 3), float(eps or 0.5), int(min_samples or 5), thr)
        # ``token`` identifies this run so the per-cluster accumulation cache and
        # the per-panel threshold callbacks stay consistent.
        return {"token": n_clicks, "passing": passing, "labels": labels}

    @app.callback(
        Output("cluster-examples", "children"),
        Input("store-cluster", "data"),
        State("store-file", "data"),
        prevent_initial_call=True,
    )
    def update_cluster_examples(cluster, path):
        if not cluster or not path:
            return []
        items = []
        for label, n_events, e_max in controller.cluster_panels(path, cluster):
            step = e_max / 100 if e_max > 0 else 0.1
            fig = controller.cluster_scene(path, cluster, label, 0.0)
            name = "unclustered" if label < 0 else f"cluster {label}"
            items.append(html.Div(style={"flex": "0 0 480px"}, children=[
                html.H5(f"{name} — {n_events} events (accumulated)"),
                html.Div(style={"display": "flex", "gap": "8px",
                                "alignItems": "center"}, children=[
                    html.Span("E threshold:"),
                    html.Div(style={"flex": "1"}, children=[
                        dcc.Slider(
                            id={"type": "cluster-thr", "index": label},
                            min=0, max=e_max, value=0, step=step,
                            tooltip={"placement": "bottom",
                                     "always_visible": False}),
                    ]),
                ]),
                dcc.Graph(id={"type": "cluster-graph", "index": label},
                          figure=fig, style={"height": "420px"}),
            ]))
        return items

    @app.callback(
        Output({"type": "cluster-graph", "index": MATCH}, "figure"),
        Input({"type": "cluster-thr", "index": MATCH}, "value"),
        State({"type": "cluster-thr", "index": MATCH}, "id"),
        State("store-cluster", "data"),
        State("store-file", "data"),
        prevent_initial_call=True,
    )
    def update_cluster_threshold(threshold, ident, cluster, path):
        if not cluster or not path:
            raise PreventUpdate
        return controller.cluster_scene(path, cluster, ident["index"],
                                        threshold or 0.0)

    @app.callback(
        Output("cluster-scatter", "figure"),
        Input("scatter-x", "value"),
        Input("scatter-y", "value"),
        Input("store-cluster", "data"),
        Input("store-cuts", "data"),
        Input("store-file", "data"),
        Input("store-hit-threshold", "data"),
    )
    def update_scatter(xvar, yvar, cluster, cuts, path, hit_threshold):
        if not path or not xvar or not yvar:
            return _empty_fig("Select x and y variables")
        thr = float(hit_threshold or 0.0)
        passing = labels = None
        if cluster:
            passing, labels = cluster["passing"], cluster["labels"]
        return controller.cluster_scatter(
            path, xvar, yvar, passing, labels, CutModel.from_store(cuts), thr)

    # ---------------------------------------------- "computing" status badge --
    # Show the badge the instant the MIP slider changes (clientside = no server
    # round-trip), then hide it as soon as the threshold-dependent figures have
    # been rebuilt. This guarantees the note is visible for the whole compute
    # cycle, even when it is fast (pre-computed valcache branches).
    app.clientside_callback(
        """
        function(value) {
            var on = value && value > 0;
            return {display: on ? 'inline-block' : 'none',
                    fontSize: '12px', fontWeight: '600', color: '#b8860b',
                    marginLeft: '8px', whiteSpace: 'nowrap'};
        }
        """,
        Output("compute-status", "style"),
        Input("hit-energy-slider", "value"),
    )

    app.clientside_callback(
        "function(){ return {display: 'none'}; }",
        Output("compute-status", "style", allow_duplicate=True),
        Input("scene3d", "figure"),
        Input("layers2d", "figure"),
        Input("dist-hist", "figure"),
        Input("cluster-scatter", "figure"),
        prevent_initial_call=True,
    )
