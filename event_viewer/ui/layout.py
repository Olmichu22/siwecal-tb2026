"""
Dash layout: the file bar, the "Event" tab and the "Distributions" tab.

Only the static structure lives here; all behaviour is wired in
:mod:`event_viewer.ui.callbacks`. Component ids are kept descriptive and stable
because the callbacks reference them by name.
"""

from __future__ import annotations

from typing import Optional

from dash import dash_table, dcc, html

from ..analysis.clustering import ALGORITHMS


def _loading(*children):
    """Wrap heavy graphs so a spinner overlays them *only while* the callbacks
    that build them are running (e.g. live metric recompute), dimming the old
    figure and clearing automatically when they finish. The textual
    ``compute-status`` badge next to the slider complements this."""
    return dcc.Loading(
        type="circle",
        overlay_style={"visibility": "visible", "opacity": 0.5},
        children=list(children),
    )


def _file_bar(controller, initial_path: Optional[str]):
    options = [{"label": _short(p), "value": p} for p in controller.list_files()]
    return html.Div(className="file-bar", children=[
        html.Span("File:", className="label"),
        dcc.Dropdown(id="file-dropdown", options=options, value=initial_path,
                     style={"width": "640px"}, placeholder="Select a .root file"),
        dcc.Input(id="file-path", type="text", placeholder="...or paste a path",
                  style={"width": "320px"}),
        html.Button("Load", id="load-btn", n_clicks=0),
        html.Span(id="file-status", className="status"),
    ], style={"display": "flex", "gap": "10px", "alignItems": "center",
              "padding": "8px", "flexWrap": "wrap"})


def _event_tab():
    return dcc.Tab(label="Event", value="tab-event", children=[
        html.Div(style={"display": "flex", "gap": "12px", "marginTop": "8px"},
                 children=[
            html.Div(style={"flex": "3", "minWidth": "0"}, children=[
                html.Div(className="event-nav", children=[
                    html.Button("◀ Prev", id="prev-btn", n_clicks=0),
                    dcc.Input(id="event-input", type="number", min=1, step=1,
                              value=1, style={"width": "90px"}),
                    html.Span(id="event-label"),
                    html.Button("Next ▶", id="next-btn", n_clicks=0),
                    dcc.Checklist(
                        id="color-clip", options=[
                            {"label": " clip E<0", "value": "clip"}],
                        value=["clip"], style={"display": "inline-block"}),
                ], style={"display": "flex", "gap": "10px",
                          "alignItems": "center", "marginBottom": "6px"}),
                html.Div(style={"display": "flex", "gap": "10px",
                                "alignItems": "center", "marginBottom": "10px"},
                         children=[
                    html.Span("MIP cut  hit_energy ≥",
                              style={"fontSize": "13px", "whiteSpace": "nowrap"}),
                    html.Div(style={"width": "220px"}, children=[
                        dcc.Slider(
                            id="hit-energy-slider",
                            min=0, max=1.0, step=0.5,
                            value=0.0,
                            marks={0: "0", 0.5: "0.5 MIP", 1.0: "1 MIP"},
                            updatemode="mouseup",
                        ),
                    ]),
                    html.Span(id="compute-status", children="⏳ Computing metrics…",
                              style={"display": "none"}),
                ]),
                _loading(
                    dcc.Graph(id="scene3d", style={"height": "460px"}),
                    dcc.Graph(id="layers2d", style={"height": "640px"}),
                ),
            ]),
            html.Div(style={"flex": "1", "minWidth": "260px"}, children=[
                html.H4("Event summary"),
                dash_table.DataTable(
                    id="metrics-table",
                    columns=[{"name": "variable", "id": "variable"},
                             {"name": "value", "id": "value"}],
                    style_cell={"fontSize": "12px", "textAlign": "left"},
                    style_table={"height": "1040px", "overflowY": "auto"},
                    page_size=60),
            ]),
        ]),
    ])


def _distributions_tab():
    algo_options = [{"label": v, "value": k} for k, v in ALGORITHMS.items()]
    return dcc.Tab(label="Distributions", value="tab-dist", children=[
      html.Div(children=[
        html.Div(style={"display": "flex", "gap": "16px", "marginTop": "8px"},
                 children=[
            # Left: histogram + dynamic cuts.
            html.Div(style={"flex": "1"}, children=[
                html.H4("Distribution"),
                dcc.Dropdown(id="dist-var", placeholder="variable"),
                html.Div(style={"display": "flex", "gap": "10px",
                                "alignItems": "center", "margin": "4px 0"},
                         children=[
                    html.Span("bins:"),
                    dcc.Input(id="dist-nbins", type="number", min=5, max=400,
                              step=1, value=60, style={"width": "80px"}),
                    dcc.Checklist(
                        id="dist-stack", options=[
                            {"label": " split by cluster", "value": "stack"}],
                        value=["stack"], style={"display": "inline-block"}),
                ]),
                _loading(dcc.Graph(id="dist-hist", style={"height": "360px"})),
                html.H4("Dynamic cuts"),
                dcc.Dropdown(id="cut-vars", multi=True,
                             placeholder="variables to cut on"),
                html.Div(id="cut-sliders"),
            ]),
            # Right: clustering.
            html.Div(style={"flex": "1"}, children=[
                html.H4("Clustering"),
                html.Div(style={"display": "flex", "gap": "8px",
                                "flexWrap": "wrap", "alignItems": "center"},
                         children=[
                    html.Span("features:"),
                    dcc.Dropdown(id="cluster-features", multi=True,
                                 style={"minWidth": "300px"}),
                ]),
                html.Div(style={"display": "flex", "gap": "8px",
                                "alignItems": "center", "marginTop": "6px",
                                "flexWrap": "wrap"}, children=[
                    dcc.Dropdown(id="cluster-algo", options=algo_options,
                                 value="kmeans", style={"width": "180px"},
                                 clearable=False),
                    html.Span("n_clusters:"),
                    dcc.Input(id="cluster-nclusters", type="number", min=2,
                              step=1, value=3, style={"width": "70px"}),
                    html.Span("eps:"),
                    dcc.Input(id="cluster-eps", type="number", min=0.01,
                              step=0.05, value=0.5, style={"width": "70px"}),
                    html.Span("min_samples:"),
                    dcc.Input(id="cluster-minsamples", type="number", min=1,
                              step=1, value=5, style={"width": "70px"}),
                    html.Button("Run", id="cluster-run", n_clicks=0),
                    html.Button("Reset", id="cluster-reset", n_clicks=0),
                ]),
                html.Div(id="cluster-param-hint",
                         style={"fontSize": "12px", "color": "#666",
                                "marginTop": "4px"}),
                html.Div(style={"display": "flex", "gap": "8px",
                                "marginTop": "8px"}, children=[
                    html.Span("x:"),
                    dcc.Dropdown(id="scatter-x", style={"width": "200px"}),
                    html.Span("y:"),
                    dcc.Dropdown(id="scatter-y", style={"width": "200px"}),
                ]),
                _loading(dcc.Graph(id="cluster-scatter",
                                   style={"height": "460px"})),
            ]),
        ]),
        # Full-width row of accumulated 3-D scenes, one per cluster.
        html.H4("Cluster examples (accumulated hits)"),
        _loading(html.Div(id="cluster-examples",
                          style={"display": "flex", "flexWrap": "wrap",
                                 "gap": "12px"})),
      ]),
    ])


def build_layout(controller, initial_path: Optional[str] = None) -> html.Div:
    """Assemble the full page (file bar + tabs + state stores)."""
    return html.Div(style={"fontFamily": "sans-serif", "padding": "8px"},
                    children=[
        html.H2("SiW-ECAL Event Viewer"),
        _file_bar(controller, initial_path),
        # Lightweight session state.
        dcc.Store(id="store-file", data=initial_path),
        dcc.Store(id="store-pos", data=0),          # position within passing list
        dcc.Store(id="store-cuts", data=[]),        # list of cut dicts
        dcc.Store(id="store-cluster", data=None),   # {passing:[...], labels:[...]}
        dcc.Store(id="store-hit-threshold", data=0.0),  # hit_energy >= threshold
        dcc.Tabs(id="tabs", value="tab-event", children=[
            _event_tab(),
            _distributions_tab(),
        ]),
    ])


def _short(path: str) -> str:
    """Display label for a file path: last two components (parent/filename)."""
    import os
    parts = path.replace("\\", "/").split("/")
    return "/".join(parts[-2:]) if len(parts) >= 2 else path
