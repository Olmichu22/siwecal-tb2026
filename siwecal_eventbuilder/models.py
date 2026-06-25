"""
Plain data containers passed between the event-building stages.

Keeping these as small dataclasses (instead of bare tuples or dicts) makes the
data flow self-documenting: every field has a name and a meaning, and the
derived quantities written to the output tree are computed by clearly named
properties rather than ad-hoc expressions.
"""

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple


# A ``(slab_position, chip_index)`` pair used as a dictionary key.
ChipKey = Tuple[int, int]


@dataclass
class Hit:
    """One reconstructed channel hit inside an event.

    Attributes mirror the per-hit branches written to the output ``ecal`` tree.
    """

    slab_position: int          # physical slab slot (ib, 0..14)  -> hit_slab
    chip_id: int                # hardware chip ID (chipid)        -> hit_chip
    channel: int                # pixel index (ipix, 0..63)        -> hit_chan
    sca: int                    # SCA memory cell (isca, 0..14)    -> hit_sca
    adc_high_pedsub: float      # high-gain ADC minus pedestal     -> hit_hg
    adc_low_pedsub: float       # low-gain ADC minus pedestal      -> hit_lg
    energy_mip: float           # high-gain signal in MIP units    -> hit_energy
    x: float = math.nan         # transverse position in mm        -> hit_x
    y: float = math.nan         # transverse position in mm        -> hit_y
    z: float = math.nan         # longitudinal position in mm      -> hit_z
    is_masked: bool = False     # masked in ped or MIP file         -> hit_ismasked


@dataclass
class BcidWindow:
    """A merged BCID time window and the chip/SCA cells that fall inside it.

    This is the intermediate product of the BCID clustering stage: it knows
    *when* an event happened (its BCID label) and *which* readout cells were
    active, but not yet the calibrated channel hits.
    """

    bcid_label: int                      # overflow-corrected BCID (event timestamp)
    start_raw: int                       # first raw BCID in the merged window
    stop_raw: int                        # last raw BCID in the merged window
    chip_to_scas: Dict[ChipKey, List[int]]  # active SCAs per (slab, chip)

    @property
    def n_slabs(self) -> int:
        """Number of distinct slabs contributing to this window."""
        return len({slab for (slab, _chip) in self.chip_to_scas})


@dataclass
class ReconstructedEvent:
    """A fully reconstructed event: a BCID timestamp plus its calibrated hits."""

    bcid: int
    hits: List[Hit]

    @property
    def n_channels(self) -> int:
        """Total number of channel hits (``nhit_chan``)."""
        return len(self.hits)

    @property
    def n_slabs(self) -> int:
        """Number of distinct slabs with at least one hit (``nhit_slab``)."""
        return len({hit.slab_position for hit in self.hits})

    @property
    def n_chips(self) -> int:
        """Number of distinct ``(slab, chip)`` pairs with hits (``nhit_chip``)."""
        return len({(hit.slab_position, hit.chip_id) for hit in self.hits})

    @property
    def sum_adc_high(self) -> float:
        """Sum of pedestal-subtracted high-gain ADC over all hits (``sum_hg``)."""
        return sum(hit.adc_high_pedsub for hit in self.hits)

    @property
    def sum_energy(self) -> float:
        """Sum of calibrated energy over all hits (``sum_energy``)."""
        return sum(hit.energy_mip for hit in self.hits)
