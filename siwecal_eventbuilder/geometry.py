"""
Detector geometry and flat-array index arithmetic for the SiW-ECAL prototype.

The converted ROOT tree stores its per-cell quantities in flat C-style arrays
whose multidimensional shape is implied by the detector geometry. This module
centralises both the geometry numbers and the index arithmetic needed to map a
``(slab, chip, sca, channel)`` coordinate onto a position in those flat arrays.

Coordinate vocabulary
----------------------
slab     : physical layer slot in the prototype, indexed by ``ib`` in 0..14.
           (The *position* in the stack, not the hardware ID of the board.)
chip     : SKIROC2 readout ASIC on a slab, indexed by ``ic`` in 0..15.
sca      : "switched-capacitor array" memory cell inside a chip, indexed by
           ``isca`` (a.k.a. ``icol``) in 0..14. Each chip can buffer up to 15
           triggered time slices before being read out.
channel  : pixel / pad read out by a chip, indexed by ``ipix`` in 0..63.
"""

from dataclasses import dataclass, fields, replace
from typing import Mapping, Optional

import os
import yaml

# Radiation length of tungsten [mm].
W_X0_MM = 3.5

# Per-slab physical z position [mm] along the beam (hit_slab 0..14). Mirrors
# event_display/conversion/slab_z_positions.yml; the live file can override it
# (see :func:`load_slab_z_mm`).
DEFAULT_SLAB_Z_MM = (0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0, 105.0,
                     120.0, 135.0, 150.0, 180.0, 195.0, 210.0, 225.0)

# Thickness [mm] of the W absorber plate in front of each slab (hit_slab 0..14).
# Mirrors the w_thickness_mm key in slab_z_positions.yml.
DEFAULT_SLAB_W_THICKNESS_MM = (2.8, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2,
                                4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6)


def load_slab_z_mm(path: str) -> tuple:
    """Read the per-slab z positions [mm] from a ``slab_z_positions`` YAML."""
    with open(path) as handle:
        document = yaml.safe_load(handle) or {}
    return tuple(float(z) for z in document.get("slab_z_mm", ()))


def load_slab_w_thickness_mm(path: str) -> tuple:
    """Read the per-slab W absorber thicknesses [mm] from a ``slab_z_positions`` YAML."""
    with open(path) as handle:
        document = yaml.safe_load(handle) or {}
    return tuple(float(t) for t in document.get("w_thickness_mm", ()))


# --------------------------------------------------------------------------- #
# Per-slab technology (the optional block of slab_z_positions.yml)
# --------------------------------------------------------------------------- #

DEFAULT_TECHNOLOGY = "FEV10"
DEFAULT_SENSOR_THICKNESS_UM = 500
DEFAULT_PAD_MAP = "fev10_rotate_chip_channel_x_y_mapping.txt"


@dataclass(frozen=True)
class SlabTechnology:
    """What each slab physically is: board technology, sensor thickness, pad
    map, and whether it ran at a threshold DAC of its own (-1 = the run's).

    Slab 12 of the 2026 stack is the FEV11 chip-on-board with its own pad map
    and its own DAC (243 against 215-230); this is where that fact lives, so
    nothing downstream carries a hard-coded 12.  Without the block in the YAML
    every slab is a FEV10 with a 500 um sensor and the default map.
    """

    technology: tuple
    sensor_thickness_um: tuple
    pad_map: tuple            # per slab, absolute path
    threshold_dac: tuple
    has_technology_block: bool = False

    def default_pad_map(self) -> str:
        counts: dict = {}
        for p in self.pad_map:
            counts[p] = counts.get(p, 0) + 1
        return max(counts, key=lambda k: (counts[k], k == self.pad_map[0]))

    def pad_map_overrides(self) -> dict:
        """``{slab: path}`` for the slabs whose map is not the default one --
        exactly the ``PadMapSlabOverrides`` list of the event builder."""
        default = self.default_pad_map()
        return {s: p for s, p in enumerate(self.pad_map) if p != default}

    def pad_map_overrides_env(self) -> str:
        """The same, as the ``EVBLD_PADMAP_SLAB_OVERRIDES`` value (``12:path,...``)."""
        return ",".join(f"{s}:{p}" for s, p in sorted(self.pad_map_overrides().items()))

    def pad_map_files(self) -> dict:
        """``{"default": path, slab: path, ...}`` -- the ``pad_map_files`` shape
        of the event viewer and the CLI config."""
        out = {"default": self.default_pad_map()}
        out.update(self.pad_map_overrides())
        return out

    def slabs_with_threshold_override(self) -> dict:
        return {s: d for s, d in enumerate(self.threshold_dac) if d >= 0}


def load_slab_technology(path: str, mappings_dir: Optional[str] = None,
                         n_slabs: int = 15) -> SlabTechnology:
    """Read the technology block of a ``slab_z_positions`` YAML; pad maps are
    resolved against ``mappings_dir`` (default: the YAML's own directory)."""
    with open(path) as handle:
        document = yaml.safe_load(handle) or {}
    mappings_dir = mappings_dir or os.path.dirname(os.path.abspath(path))
    technologies = document.get("technologies") or {}

    def per_slab(key, default):
        values = document.get(key)
        if values is None:
            return [default] * n_slabs
        if len(values) != n_slabs:
            raise ValueError(f"{path}: {key} has {len(values)} entries, expected {n_slabs}")
        return list(values)

    tech = per_slab("slab_technology", DEFAULT_TECHNOLOGY)
    for s, t in enumerate(tech):
        if t not in technologies and t != DEFAULT_TECHNOLOGY:
            raise ValueError(f"{path}: slab {s} technology '{t}' is not declared under 'technologies'")

    def tech_default(t, key, fallback):
        return (technologies.get(t) or {}).get(key, fallback)

    thickness = document.get("sensor_thickness_um")
    if thickness is None:
        thickness = [tech_default(t, "sensor_thickness_um", DEFAULT_SENSOR_THICKNESS_UM) for t in tech]
    thickness = [int(v) for v in per_slab("sensor_thickness_um", None)] if "sensor_thickness_um" in document \
        else [int(v) for v in thickness]
    pad_map = [tech_default(t, "pad_map", DEFAULT_PAD_MAP) for t in tech]
    pad_map = [p if os.path.isabs(p) else os.path.join(mappings_dir, p) for p in pad_map]
    dac = [int(v) for v in per_slab("threshold_dac", -1)]
    return SlabTechnology(technology=tuple(tech), sensor_thickness_um=tuple(thickness),
                         pad_map=tuple(pad_map), threshold_dac=tuple(dac),
                         has_technology_block="slab_technology" in document)


@dataclass(frozen=True)
class DetectorGeometry:
    """Immutable description of the detector's logical array dimensions.

    The default values match the TB2026CERN 15-slab configuration. They are kept
    here (rather than scattered as module constants) so that every component
    derives its index arithmetic from a single, explicit source of truth.
    """

    n_slab_positions: int = 15      # physical slab slots (ib)
    n_chips_per_slab: int = 16      # SKIROC2 ASICs per slab (ic)
    n_scas_per_chip: int = 15       # SCA memory cells per chip (isca / icol)
    n_channels_per_chip: int = 64   # pixels per chip (ipix)
    slab_z_mm: tuple = DEFAULT_SLAB_Z_MM              # per-slab z [mm] -> hit_z
    slab_w_thickness_mm: tuple = DEFAULT_SLAB_W_THICKNESS_MM  # W absorber thickness [mm] per slab

    @classmethod
    def from_mapping(cls, overrides: Optional[Mapping] = None) -> "DetectorGeometry":
        """Build a geometry from defaults, overriding only the keys provided.

        Mirrors :meth:`BuilderConfig.from_mapping` for the optional ``geometry:``
        section of ``config.yml``. Absent keys keep their default; ``None`` or an
        empty mapping returns the plain defaults.

        Raises
        ------
        ValueError
            If a key does not name a geometry field (lists the valid names).
        """
        if not overrides:
            return cls()

        valid_names = {f.name for f in fields(cls)}
        unknown = [key for key in overrides if key not in valid_names]
        if unknown:
            raise ValueError(
                "Unknown DetectorGeometry option(s) in config file: "
                f"{', '.join(sorted(unknown))}. "
                f"Valid options are: {', '.join(sorted(valid_names))}."
            )
        converted = {}
        for key, value in overrides.items():
            if key in ("slab_z_mm", "slab_w_thickness_mm"):
                converted[key] = tuple(float(item) for item in value)
            else:
                converted[key] = int(value)
        return replace(cls(), **converted)

    def slab_z(self, slab: int) -> float:
        """Physical z [mm] of a slab along the beam, or NaN if out of range."""
        if 0 <= slab < len(self.slab_z_mm):
            return float(self.slab_z_mm[slab])
        return float("nan")

    def slab_x0(self, slab: int) -> float:
        """Cumulative radiation lengths of W traversed up to and including slab ``slab``.

        Sums the W absorber thicknesses for slabs 0..slab and divides by X0(W)=3.5 mm.
        Returns NaN if ``slab`` is out of range.
        """
        if 0 <= slab < len(self.slab_w_thickness_mm):
            return sum(self.slab_w_thickness_mm[:slab + 1]) / W_X0_MM
        return float("nan")

    def sca_index(self, slab: int, chip: int, sca: int) -> int:
        """Flat index into the ``(slab, chip, sca)`` arrays.

        Applies to the per-SCA branches: ``nhits``, ``bcid``, ``corrected_bcid``
        and ``badbcid``. Equivalent to the original ``flat3`` helper.
        """
        return (slab * self.n_chips_per_slab + chip) * self.n_scas_per_chip + sca

    def channel_index(self, slab: int, chip: int, sca: int, channel: int) -> int:
        """Flat index into the ``(slab, chip, sca, channel)`` arrays.

        Applies to the per-channel branches: ``adc_high``, ``adc_low``,
        ``hitbit_high``, etc. Equivalent to the original ``flat4`` helper.
        """
        return self.sca_index(slab, chip, sca) * self.n_channels_per_chip + channel

    @property
    def n_chip_rows(self) -> int:
        """Total number of ``(slab, chip)`` rows, used to shape the BCID matrix."""
        return self.n_slab_positions * self.n_chips_per_slab
