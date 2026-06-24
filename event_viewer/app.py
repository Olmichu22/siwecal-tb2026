"""
Dash application factory.

``build_app`` wires the configuration, controller, layout and callbacks into a
ready-to-run :class:`dash.Dash` instance. Keeping it a factory (rather than a
module-level global) makes it easy to launch with a different config or to import
the app object for a WSGI server.
"""

from __future__ import annotations

from typing import Optional

import dash

from .config import ViewerConfig
from .controller import ViewerController
from .ui import build_layout, register_callbacks


def build_app(config: Optional[ViewerConfig] = None,
              initial_path: Optional[str] = None) -> dash.Dash:
    """Create the Dash app for one viewer session."""
    config = config or ViewerConfig()
    controller = ViewerController(config)

    app = dash.Dash(__name__, title="SiW-ECAL Event Viewer",
                    suppress_callback_exceptions=True)
    app.layout = build_layout(controller, initial_path)
    register_callbacks(app, controller)
    return app
