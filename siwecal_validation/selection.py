"""
Event selection cuts for the SiW-ECAL validation plots.

A :class:`CutSet` is an immutable set of optional bounds on the per-event
variables. The same cut is applied to every plot of a sample, and is reflected
in the output file names and titles.

Besides the original ``nhit``/``energy``/``rate`` bounds, **every**
particle-discrimination variable is selectable: tungsten-weighted energy,
Molière radius, transverse RMS, transverse barycenter radius, the longitudinal
shower layers (start/max/end/length), the number of layers hit, the first/last
hit layer and the hit-energy density — plus a boolean ``is_shower`` filter.

Two-level cut philosophy
------------------------
Cuts can be defined **generally** (e.g. from the CLI, applied to all energies)
and **per energy** (from the YAML), so a different selection can be optimised
for each energy point. The effective cut for an energy is::

    general_cut.merge(per_energy_cut)

where the per-energy values override the general ones wherever they are set.
"""

from dataclasses import dataclass, fields
from typing import Mapping, Optional

import numpy as np

# Order matters: it fixes the order of tokens in the file-name suffix.
# (field_min, field_max, short_name, EventData attribute)
_CUT_SPEC = (
    ("nhit_min",           "nhit_max",           "nhit",   "nhit"),
    ("energy_min",         "energy_max",         "E",      "energy"),
    ("rate_min",           "rate_max",           "rate",   "mip_likeness"),
    ("weighte_min",        "weighte_max",        "wE",     "weighte"),
    ("moliere_min",        "moliere_max",        "mol",    "moliere"),
    ("transverse_rms_min", "transverse_rms_max", "trms",   "transverse_rms"),
    ("bar_r_min",          "bar_r_max",          "barR",   "bar_r"),
    ("shower_start_min",   "shower_start_max",   "shStart", "shower_start"),
    ("shower_max_min",     "shower_max_max",     "shMax",  "shower_max"),
    ("shower_end_min",     "shower_end_max",     "shEnd",  "shower_end"),
    ("shower_length_min",  "shower_length_max",  "shLen",  "shower_length"),
    ("n_layers_hit_min",   "n_layers_hit_max",   "nLay",   "n_layers_hit"),
    ("first_layer_min",    "first_layer_max",    "first",  "first_layer"),
    ("last_layer_min",     "last_layer_max",     "last",   "last_layer"),
    ("e_over_nhit_min",    "e_over_nhit_max",    "eN",     "e_over_nhit"),
)


@dataclass(frozen=True)
class CutSet:
    """Immutable set of optional lower/upper bounds on per-event variables."""

    nhit_min: Optional[float] = None
    nhit_max: Optional[float] = None
    energy_min: Optional[float] = None
    energy_max: Optional[float] = None
    rate_min: Optional[float] = None
    rate_max: Optional[float] = None
    weighte_min: Optional[float] = None
    weighte_max: Optional[float] = None
    moliere_min: Optional[float] = None
    moliere_max: Optional[float] = None
    transverse_rms_min: Optional[float] = None
    transverse_rms_max: Optional[float] = None
    bar_r_min: Optional[float] = None
    bar_r_max: Optional[float] = None
    shower_start_min: Optional[float] = None
    shower_start_max: Optional[float] = None
    shower_max_min: Optional[float] = None
    shower_max_max: Optional[float] = None
    shower_end_min: Optional[float] = None
    shower_end_max: Optional[float] = None
    shower_length_min: Optional[float] = None
    shower_length_max: Optional[float] = None
    n_layers_hit_min: Optional[float] = None
    n_layers_hit_max: Optional[float] = None
    first_layer_min: Optional[float] = None
    first_layer_max: Optional[float] = None
    last_layer_min: Optional[float] = None
    last_layer_max: Optional[float] = None
    e_over_nhit_min: Optional[float] = None
    e_over_nhit_max: Optional[float] = None
    # Boolean filter: keep only showers (True) / non-showers (False) / all (None).
    is_shower: Optional[bool] = None

    # ------------------------------------------------------- constructors ---
    @classmethod
    def from_mapping(cls, mapping: Optional[Mapping]) -> "CutSet":
        """Build a cut set from a dict (e.g. a YAML ``cuts:`` block).

        Unknown keys raise ``ValueError`` so typos are caught early. ``is_shower``
        is coerced to bool; every other field to float.
        """
        if not mapping:
            return cls()
        valid = {f.name for f in fields(cls)}
        unknown = [key for key in mapping if key not in valid]
        if unknown:
            raise ValueError(
                f"Unknown cut option(s): {', '.join(sorted(unknown))}. "
                f"Valid options: {', '.join(sorted(valid))}.")

        def _coerce(key, value):
            if value is None:
                return None
            return bool(value) if key == "is_shower" else float(value)

        return cls(**{key: _coerce(key, value) for key, value in mapping.items()})

    def merge(self, override: "CutSet") -> "CutSet":
        """Return a new cut set where ``override``'s set values win."""
        merged = {}
        for f in fields(self):
            other = getattr(override, f.name)
            merged[f.name] = other if other is not None else getattr(self, f.name)
        return CutSet(**merged)

    # ----------------------------------------------------------- queries ---
    @property
    def is_empty(self) -> bool:
        return all(getattr(self, f.name) is None for f in fields(self))

    def mask(self, data) -> np.ndarray:
        """Boolean mask of events passing every bound, given an EventData.

        NaN-valued variables (e.g. shower layers of non-shower events) fail any
        bound comparison, so a shower-shape cut implicitly drops non-showers.
        """
        keep = np.ones(len(data), dtype=bool)
        for field_min, field_max, _name, attr in _CUT_SPEC:
            values = getattr(data, attr)
            low, high = getattr(self, field_min), getattr(self, field_max)
            if low is not None:
                keep &= values >= low
            if high is not None:
                keep &= values <= high
        if self.is_shower is not None:
            keep &= data.is_shower == self.is_shower
        return keep

    @property
    def suffix(self) -> str:
        """Compact file-name suffix encoding the active cuts (``""`` if none)."""
        tokens = []
        for field_min, field_max, name, _attr in _CUT_SPEC:
            low, high = getattr(self, field_min), getattr(self, field_max)
            if low is not None:
                tokens.append(f"{name}Min{low:g}")
            if high is not None:
                tokens.append(f"{name}Max{high:g}")
        if self.is_shower is not None:
            tokens.append("shower" if self.is_shower else "noShower")
        return ("_" + "_".join(tokens)) if tokens else ""

    @property
    def label(self) -> str:
        """Human-readable cut description for plot titles (``""`` if none)."""
        parts = []
        for field_min, field_max, name, _attr in _CUT_SPEC:
            low, high = getattr(self, field_min), getattr(self, field_max)
            if low is not None:
                parts.append(f"{name}≥{low:g}")
            if high is not None:
                parts.append(f"{name}≤{high:g}")
        if self.is_shower is not None:
            parts.append("shower" if self.is_shower else "non-shower")
        return ("\ncuts: " + ", ".join(parts)) if parts else ""
