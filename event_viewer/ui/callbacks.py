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
the matching cut-vars to ``[]`` clears it), which reads both the range sliders and
the per-cut ``compl.`` checkboxes. Adding or removing a cut variable
rebuilds the whole slider block, so the builder carries the existing selections
*and* their complement flags over instead of resetting every cut to its full
range. The "Complementary" button and the per-variable ``%`` (quantile) toggle
write only the widgets, never the store, so the store keeps its single writer:
the ``%`` toggle rebuilds the slider on the 0-100% scale and the rebuilt slider
is what feeds the store. ``store-pos`` has one primary writer
(navigation) plus a reset on file / event-cut change.
"""

from __future__ import annotations

import os

import numpy as np
import plotly.graph_objects as go
from dash import ALL, MATCH, Input, Output, State, ctx, dcc, html, no_update
from dash.exceptions import PreventUpdate

from .._timing import timed
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


def _fmt_mark(value: float) -> str:
    """Short label for a discrete slider mark (integers without ``.0``)."""
    return str(int(value)) if float(value).is_integer() else f"{value:g}"


def _discrete_marks(values) -> dict:
    """Snap marks for a discrete slider: one per value, labels thinned when
    there are too many to read (every value is still a valid snap point)."""
    n = len(values)
    show_every = 1 if n <= 12 else int(np.ceil(n / 12))
    marks = {}
    for k, value in enumerate(values):
        label = _fmt_mark(value) if k % show_every == 0 or k == n - 1 else ""
        marks[value] = {"label": label,
                        "style": {"fontSize": "10px",
                                  "whiteSpace": "nowrap"}}
    return marks


def _toggled_var(id_type: str):
    """The variable whose checkbox of ``id_type`` fired this callback, if any."""
    trigger = ctx.triggered_id
    if isinstance(trigger, dict) and trigger.get("type") == id_type:
        return trigger.get("index")
    return None


def _cut_slider(controller, path, var, thr, prev, id_type, prev_invert=(),
                prev_quantile=()):
    """One labelled cut slider (discrete if the variable has few unique values,
    otherwise a continuous 100-step range), preserving a previous selection.

    The label carries a ``compl.`` checkbox that moves this cut into the inverted
    group (see :mod:`event_viewer.analysis.cuts`), so a complementary selection
    can be built one variable at a time; the "Complementary" button flips them
    all at once.

    It also carries a ``%`` checkbox: with it on, the slider runs over 0-100%
    of the variable's own distribution instead of its value range, so the same
    selection means the same thing across runs of different beam energy. The
    toggle is disabled for discrete variables, where the slider snaps to the
    unique values and a percentile has no useful meaning."""
    lo, hi, values = controller.variable_domain(path, var, thr)
    discrete = bool(values) and len(values) > 1
    quantile = var in prev_quantile and not discrete
    if quantile:
        # 0-100% of the distribution; the store keeps fractions and the
        # percentiles are resolved against the data at mask() time.
        value = prev.get(var) or [0, 100]
        value = [max(0, min(100, value[0])), max(0, min(100, value[1]))]
        slider = dcc.RangeSlider(
            id={"type": id_type, "index": var},
            min=0, max=100, value=value, step=1, allowCross=False,
            marks={0: "0%", 25: "25%", 50: "50%", 75: "75%", 100: "100%"},
            tooltip={"placement": "bottom", "always_visible": False})
    elif discrete:
        # Discrete: snap to the unique values; carry the previous selection over
        # by snapping each end to the nearest allowed value.
        arr = np.asarray(values, dtype=float)
        prev_val = prev.get(var) or [lo, hi]
        value = [float(arr[np.abs(arr - v).argmin()]) for v in prev_val]
        slider = dcc.RangeSlider(
            id={"type": id_type, "index": var},
            min=lo, max=hi, value=value, step=None, marks=_discrete_marks(values),
            allowCross=False,
            tooltip={"placement": "bottom", "always_visible": False})
    else:
        step = (hi - lo) / 100 if hi > lo else 1.0
        value = prev.get(var) or [lo, hi]
        value = [max(lo, min(hi, value[0])), max(lo, min(hi, value[1]))]
        slider = dcc.RangeSlider(
            id={"type": id_type, "index": var},
            min=lo, max=hi, value=value, step=step, allowCross=False,
            tooltip={"placement": "bottom", "always_visible": False})
    return html.Div(style={"marginBottom": "14px"}, children=[
        html.Div(style={"display": "flex", "gap": "8px",
                        "alignItems": "baseline"}, children=[
            html.Label(var, style={"fontSize": "13px"}),
            dcc.Checklist(
                id={"type": f"{id_type}-invert", "index": var},
                options=[{"label": " compl.", "value": "invert"}],
                value=["invert"] if var in prev_invert else [],
                style={"fontSize": "11px", "color": "#b45309"},
                inputStyle={"marginRight": "3px"}),
            dcc.Checklist(
                id={"type": f"{id_type}-quantile", "index": var},
                options=[{"label": " %", "value": "quantile",
                          "disabled": discrete}],
                value=["quantile"] if quantile else [],
                style={"fontSize": "11px", "color": "#2563eb"},
                inputStyle={"marginRight": "3px"}),
        ]),
        slider])


def _flagged(ids, values):
    """The variables whose checkbox in a pattern-matching group is ticked."""
    return {ident["index"]
            for ident, val in zip(ids or [], values or []) if val}


def _cut_sliders(controller, cut_vars, path, hit_threshold, cur_values, cur_ids,
                 id_type, cur_invert=None, cur_invert_ids=None,
                 cur_quantile=None, cur_quantile_ids=None, toggled=None):
    """Rebuild every cut slider, carrying over existing selections.

    ``toggled`` is the variable whose value/quantile mode just changed, if any:
    its slider is rebuilt on the new scale, so the carried-over range (a value
    range on one scale, percentages on the other) is dropped and it reopens
    fully open instead of keeping a number that now means something else."""
    if not path or not cut_vars:
        return []
    thr = float(hit_threshold or 0.0)
    # Adding/removing a variable rebuilds every slider; carry over the values of
    # variables that were already there so their selections survive.
    prev = {ident["index"]: val
            for ident, val in zip(cur_ids or [], cur_values or [])
            if ident["index"] != toggled}
    prev_invert = _flagged(cur_invert_ids, cur_invert)
    prev_quantile = _flagged(cur_quantile_ids, cur_quantile)
    return [_cut_slider(controller, path, var, thr, prev, id_type, prev_invert,
                        prev_quantile)
            for var in cut_vars]


def _cuts_from_widgets(values, ids, invert_values, invert_ids,
                       quantile_values=None, quantile_ids=None):
    """The cut store, from the slider values and the per-cut checkboxes.

    Matched by variable name rather than by position: the pattern-matching
    groups are built in the same order today, but a store that silently pairs a
    range with another variable's checkbox would be a nasty bug to find.

    A slider in quantile mode runs over 0-100 and is stored as the fractions
    [0, 1] that :class:`~event_viewer.analysis.cuts.Cut` expects; the
    percentiles themselves stay unresolved until the cut is applied."""
    inverted = _flagged(invert_ids, invert_values)
    quantile = _flagged(quantile_ids, quantile_values)
    cuts = []
    for value, ident in zip(values, ids):
        if value is None:
            continue
        var = ident["index"]
        lo, hi = value[0], value[1]
        is_quantile = var in quantile
        if is_quantile:
            lo, hi = lo / 100.0, hi / 100.0
        cuts.append({"variable": var, "lo": lo, "hi": hi,
                     "invert": var in inverted,
                     "mode": "quantile" if is_quantile else "value"})
    return cuts


def _toggle_all(current):
    """All-on if anything is off, else all-off — one click either way."""
    if not current:
        raise PreventUpdate
    turn_on = any(not val for val in current)
    return [["invert"] if turn_on else [] for _ in current]


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
        # The ``%`` toggles rebuild the block so the slider can change scale;
        # this stays the only writer of the sliders, and the store below keeps
        # its single writer.
        Input({"type": "cut-slider-quantile", "index": ALL}, "value"),
        State("store-file", "data"),
        State("store-hit-threshold", "data"),
        State({"type": "cut-slider", "index": ALL}, "value"),
        State({"type": "cut-slider", "index": ALL}, "id"),
        State({"type": "cut-slider-invert", "index": ALL}, "value"),
        State({"type": "cut-slider-invert", "index": ALL}, "id"),
        State({"type": "cut-slider-quantile", "index": ALL}, "id"),
    )
    def build_sliders(cut_vars, cur_quantile, path, hit_threshold, cur_values, cur_ids,
                      cur_invert, cur_invert_ids, cur_quantile_ids):
        return _cut_sliders(controller, cut_vars, path, hit_threshold,
                            cur_values, cur_ids, "cut-slider",
                            cur_invert, cur_invert_ids,
                            cur_quantile, cur_quantile_ids,
                            _toggled_var("cut-slider-quantile"))

    @app.callback(
        Output("store-cuts", "data"),
        Input({"type": "cut-slider", "index": ALL}, "value"),
        Input({"type": "cut-slider-invert", "index": ALL}, "value"),
        State({"type": "cut-slider", "index": ALL}, "id"),
        State({"type": "cut-slider-invert", "index": ALL}, "id"),
        # Read as State: toggling ``%`` rebuilds the slider, and it is the new
        # slider value that re-triggers this callback, mode already set.
        State({"type": "cut-slider-quantile", "index": ALL}, "value"),
        State({"type": "cut-slider-quantile", "index": ALL}, "id"),
    )
    def update_cuts(values, invert_values, ids, invert_ids,
                    quantile_values, quantile_ids):
        return _cuts_from_widgets(values, ids, invert_values, invert_ids,
                                  quantile_values, quantile_ids)

    @app.callback(
        Output({"type": "cut-slider-invert", "index": ALL}, "value"),
        Input("cut-complement-btn", "n_clicks"),
        State({"type": "cut-slider-invert", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def toggle_complement(_n_clicks, current):
        return _toggle_all(current)

    # ------------------------------------ Event-tab dynamic cut UI (page 1) --
    # A second, independent copy of the cut widget. It writes ``store-event-cuts``
    # (which drives the one-by-one navigation and the bicolor histogram) and is
    # completely decoupled from the clustering tab's ``store-cuts``.
    @app.callback(
        Output("ev-cut-sliders", "children"),
        Input("ev-cut-vars", "value"),
        # The ``%`` toggles rebuild the block so the slider can change scale;
        # this stays the only writer of the sliders, and the store below keeps
        # its single writer.
        Input({"type": "ev-cut-slider-quantile", "index": ALL}, "value"),
        State("store-file", "data"),
        State("store-hit-threshold", "data"),
        State({"type": "ev-cut-slider", "index": ALL}, "value"),
        State({"type": "ev-cut-slider", "index": ALL}, "id"),
        State({"type": "ev-cut-slider-invert", "index": ALL}, "value"),
        State({"type": "ev-cut-slider-invert", "index": ALL}, "id"),
        State({"type": "ev-cut-slider-quantile", "index": ALL}, "id"),
    )
    def build_sliders_event(cut_vars, cur_quantile, path, hit_threshold, cur_values, cur_ids,
                            cur_invert, cur_invert_ids, cur_quantile_ids):
        return _cut_sliders(controller, cut_vars, path, hit_threshold,
                            cur_values, cur_ids, "ev-cut-slider",
                            cur_invert, cur_invert_ids,
                            cur_quantile, cur_quantile_ids,
                            _toggled_var("ev-cut-slider-quantile"))

    @app.callback(
        Output("store-event-cuts", "data"),
        Input({"type": "ev-cut-slider", "index": ALL}, "value"),
        Input({"type": "ev-cut-slider-invert", "index": ALL}, "value"),
        State({"type": "ev-cut-slider", "index": ALL}, "id"),
        State({"type": "ev-cut-slider-invert", "index": ALL}, "id"),
        # Read as State: toggling ``%`` rebuilds the slider, and it is the new
        # slider value that re-triggers this callback, mode already set.
        State({"type": "ev-cut-slider-quantile", "index": ALL}, "value"),
        State({"type": "ev-cut-slider-quantile", "index": ALL}, "id"),
    )
    def update_cuts_event(values, invert_values, ids, invert_ids,
                          quantile_values, quantile_ids):
        return _cuts_from_widgets(values, ids, invert_values, invert_ids,
                                  quantile_values, quantile_ids)

    @app.callback(
        Output({"type": "ev-cut-slider-invert", "index": ALL}, "value"),
        Input("ev-cut-complement-btn", "n_clicks"),
        State({"type": "ev-cut-slider-invert", "index": ALL}, "value"),
        prevent_initial_call=True,
    )
    def toggle_complement_event(_n_clicks, current):
        return _toggle_all(current)

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
        Output("tracks-status", "children"),
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
                    "no file", None, 1, "")
        thr = float(hit_threshold or 0.0)
        passing = controller.passing_indices(path, CutModel.from_store(cuts), thr)
        n_pass = len(passing)
        if n_pass == 0:
            return (_empty_fig("No events pass the cuts"), _empty_fig(),
                    [], "0 events passing cuts", None, 1, "")
        pos = _pos_of(passing, cur_index)
        index = int(passing[pos])
        color_clip = "clip" in (clip or [])
        overlays = overlays or []
        scene, layers, rows, tracks_msg = controller.event_figures(
            path, index, color_clip, thr,
            show_moliere="moliere" in overlays, show_axis="axis" in overlays,
            show_tracks="tracks" in overlays)
        ds = controller.dataset(path)
        label = (f"event {pos + 1} / {n_pass} passing "
                 f"(entry {index}, {ds.n_events} total)")
        return scene, layers, rows, label, pos + 1, n_pass, tracks_msg

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
        with timed("update_cluster_examples (all panels)") as info:
            panels = controller.cluster_panels(path, cluster)
            for label, n_events, e_max in panels:
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
            info["panels"] = len(panels)
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
