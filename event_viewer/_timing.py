"""Lightweight timing / progress logging for the viewer's hot paths.

The per-cluster accumulated 3-D scenes ("clustering histograms") can be slow to
build the first time, mostly because reading every hit of an EDM4hep run goes
through per-hit podio accessors. To tell a slow-but-progressing build from a hung
one -- and to measure where the time actually goes -- the hot stages are wrapped
in :func:`timed` and the long per-hit read logs periodic :func:`progress` lines.

Logging is on at INFO by default (the CLI calls :func:`configure`). Silence it
with ``EVENT_VIEWER_TIMING=0`` in the environment, or raise detail elsewhere by
setting the ``event_viewer`` logger's level.
"""

from __future__ import annotations

import logging
import os
import time
from contextlib import contextmanager

logger = logging.getLogger("event_viewer")

# A single opt-out switch, read once. Any value other than "0"/"false" keeps the
# timings on; the CLI still has to call configure() to attach the handler.
_ENABLED = os.environ.get("EVENT_VIEWER_TIMING", "1").lower() not in ("0", "false", "no")


def configure(level: int = logging.INFO) -> None:
    """Attach a stderr handler (once) so timing lines show on the console."""
    if not _ENABLED:
        logger.setLevel(logging.WARNING)
        return
    if not any(getattr(h, "_event_viewer_timing", False) for h in logger.handlers):
        handler = logging.StreamHandler()
        handler._event_viewer_timing = True  # type: ignore[attr-defined]
        handler.setFormatter(logging.Formatter(
            "%(asctime)s.%(msecs)03d  event_viewer  %(message)s",
            datefmt="%H:%M:%S"))
        logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False


@contextmanager
def timed(label: str):
    """Time a block, logging its start and elapsed ms.

    Yields a mutable ``dict`` the caller can fill with context (e.g. counts);
    those key=value pairs are appended to the completion line::

        with timed("accumulate") as info:
            ...
            info["events"] = n; info["pads"] = m
    """
    if not _ENABLED:
        yield {}
        return
    info: dict = {}
    logger.info("→ %s ...", label)
    t0 = time.perf_counter()
    try:
        yield info
    finally:
        dt = (time.perf_counter() - t0) * 1000.0
        extra = "  ".join(f"{k}={v}" for k, v in info.items())
        logger.info("✓ %s  %.0f ms%s", label, dt, f"  [{extra}]" if extra else "")


def progress(label: str, done: int, total: int, t0: float, **fields) -> None:
    """Log a one-line progress heartbeat for a long loop (rate + ETA)."""
    if not _ENABLED:
        return
    elapsed = time.perf_counter() - t0
    rate = done / elapsed if elapsed > 0 else 0.0
    eta = (total - done) / rate if rate > 0 else float("inf")
    extra = "  ".join(f"{k}={v}" for k, v in fields.items())
    logger.info("  %s: %d/%d  %.1fs  %.0f/s  eta %.0fs%s",
                label, done, total, elapsed, rate, eta,
                f"  {extra}" if extra else "")
