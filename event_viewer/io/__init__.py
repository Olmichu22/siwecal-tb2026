"""I/O layer: reading reconstructed ``ecal`` ROOT files with uproot."""

from .valcache_reader import EventFileReader

__all__ = ["EventFileReader"]
