"""
Command-line entry point: ``python -m event_viewer [--file PATH] [--port N]``.

Launches the Dash development server. Access it from your laptop with an SSH
tunnel, e.g.::

    ssh -L 8050:localhost:8050 you@lxplus.cern.ch
    # then open http://localhost:8050
"""

from __future__ import annotations

import argparse

from .app import build_app
from .config import ViewerConfig


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="event_viewer",
        description="Interactive SiW-ECAL event viewer (Plotly Dash).")
    parser.add_argument("--file", default=None,
                        help="ROOT file to open on startup (.valcache.root or ecal_*.root)")
    parser.add_argument("--data-dir", default=None,
                        help="Directory scanned for the file dropdown (default: <project>/data)")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8050)
    parser.add_argument("--debug", action="store_true",
                        help="Run Dash in debug mode (hot reload).")
    args = parser.parse_args()

    config = ViewerConfig()
    if args.data_dir:
        config.data_dir = args.data_dir
    config.host, config.port = args.host, args.port

    app = build_app(config, initial_path=args.file)
    app.run(host=args.host, port=args.port, debug=args.debug)


if __name__ == "__main__":
    main()
