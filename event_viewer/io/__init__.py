"""I/O layer: reading reconstructed events (legacy ``ecal`` tree or EDM4hep)."""

from .edm4hep_reader import Edm4hepEventReader
from .valcache_reader import EventFileReader


def is_edm4hep(path: str) -> bool:
    """True if ``path`` is a podio/EDM4hep file (has an ``events`` category)
    rather than a flat ``ecal`` TTree."""
    import uproot
    with uproot.open(path) as handle:
        keys = {k.split(";")[0] for k in handle.keys()}
    # podio files carry the per-event 'events' tree; the legacy files an 'ecal'.
    return "events" in keys and "ecal" not in keys


def make_reader(path: str, tree_name: str = "ecal", n_layers: int = 15):
    """Return the right reader for ``path`` (EDM4hep PID file or legacy ecal)."""
    if is_edm4hep(path):
        return Edm4hepEventReader(path, n_layers=n_layers)
    return EventFileReader(path, tree_name, n_layers)


__all__ = ["EventFileReader", "Edm4hepEventReader", "make_reader", "is_edm4hep"]
