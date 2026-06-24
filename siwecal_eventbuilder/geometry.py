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

import yaml

# Per-slab physical z position [mm] along the beam (hit_slab 0..14). Mirrors
# event_display/conversion/slab_z_positions.yml; the live file can override it
# (see :func:`load_slab_z_mm`).
DEFAULT_SLAB_Z_MM = (0.0, 11.0, 22.0, 33.0, 44.0, 55.0, 66.0, 77.0,
                     88.0, 99.0, 110.0, 132.0, 143.0, 154.0, 165.0)


def load_slab_z_mm(path: str) -> tuple:
    """Read the per-slab z positions [mm] from a ``slab_z_positions`` YAML."""
    with open(path) as handle:
        document = yaml.safe_load(handle) or {}
    return tuple(float(z) for z in document.get("slab_z_mm", ()))


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
    slab_z_mm: tuple = DEFAULT_SLAB_Z_MM   # per-slab z position [mm] -> hit_z

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
            if key == "slab_z_mm":
                converted[key] = tuple(float(item) for item in value)
            else:
                converted[key] = int(value)
        return replace(cls(), **converted)

    def slab_z(self, slab: int) -> float:
        """Physical z [mm] of a slab along the beam, or NaN if out of range."""
        if 0 <= slab < len(self.slab_z_mm):
            return float(self.slab_z_mm[slab])
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
