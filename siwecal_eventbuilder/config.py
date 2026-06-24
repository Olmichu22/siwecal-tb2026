"""
Tunable parameters for the SiW-ECAL event builder.

Every threshold, cut value and processing knob lives in :class:`BuilderConfig`
so that the algorithmic code never hard-codes a magic number. The defaults
reproduce the behaviour validated against tkamiyam's reference event builder
(3191 events on run 7 / 74 GeV).
"""

import math
from dataclasses import dataclass, field, fields, replace
from typing import FrozenSet, Mapping, Optional


@dataclass(frozen=True)
class BuilderConfig:
    """Immutable bag of configuration values shared by all components."""

    # ------------------------------------------------------------------ I/O ---
    tree_name: str = "siwecaldecoded"
    """Name of the input TTree produced by the RAW2ROOT converter."""

    # ------------------------------------------------------- BCID selection ---
    skip_bcid_start: int = 50
    """Reject BCIDs below this value: the start of an acquisition is noisy."""

    drop_bcids: FrozenSet[int] = field(default_factory=lambda: frozenset({0, 901}))
    """Specific BCID values that are known artefacts and always discarded."""

    merge_delta: int = 3
    """Two consecutive BCIDs closer than this are merged into one time window."""

    min_slabs_hit: int = 10
    """A BCID window is kept only if it spans at least this many distinct slabs."""

    bcid_overflow: int = 4096
    """The BCID counter is 12-bit; it wraps (overflows) every 4096 counts."""

    bad_value: int = -999
    """Sentinel written by the converter for empty / invalid SCA cells."""

    # -------------------------------------------------------- hit selection ---
    adc_underflow_threshold: int = 11
    """High-gain ADC values <= this are treated as underflow and dropped."""

    # ------------------------------ calibration-only quality cuts (unused in --
    #                                event building, kept for pedestal/MIP) -----
    badbcid_max_good: int = 999
    """Upper bound on ``badbcid`` accepted during MIP calibration."""

    max_hits_per_sca: float = math.inf
    """Upper bound on ``nhits`` per SCA during MIP calibration (disabled)."""

    pedestal_fallback: float = 250.0
    """Pedestal assumed when a channel is missing from the calibration map."""

    default_mip_fallback: float = 20.0
    """MIP value assumed when no channel could be calibrated at all."""

    # ----------------------------------------------------------- processing ---
    max_hits_per_event: int = 15360
    """Hard cap on hits per event, sizing the writer's fixed per-hit buffers.

    Set to the total channel count (15 slabs x 16 chips x 64 channels = 15360):
    since a channel can fire at most once per event, no physical event can exceed
    it, so this cap never drops an event -- it only guards against buffer
    overflow. Lower it only if memory is a concern and you accept losing
    pathological high-multiplicity events.
    """

    default_workers: int = 5
    """Default number of parallel worker processes per run."""

    # ------------------------------------------------------ YAML overrides ----
    @classmethod
    def from_mapping(cls, overrides: Optional[Mapping] = None) -> "BuilderConfig":
        """Build a config from defaults, overriding only the keys provided.

        This is the bridge between the optional ``config.yml`` file and the
        immutable dataclass. Absent keys keep their default value, so the file
        only needs to list what the user wants to change. Passing ``None`` or an
        empty mapping returns the plain defaults.

        Parameters
        ----------
        overrides:
            Mapping whose keys must match :class:`BuilderConfig` field names
            (typically the ``builder:`` section of ``config.yml``). Two values
            receive light normalisation so the YAML stays natural to write:

            * ``drop_bcids`` -- any iterable (e.g. a YAML list ``[0, 901]``) is
              converted to a ``frozenset`` of ints.
            * ``max_hits_per_sca`` -- coerced to ``float`` so ``.inf`` works.

        Raises
        ------
        ValueError
            If a key does not name a configuration field; the message lists the
            valid field names so typos are caught early.
        """
        if not overrides:
            return cls()

        valid_names = {f.name for f in fields(cls)}
        unknown = [key for key in overrides if key not in valid_names]
        if unknown:
            raise ValueError(
                "Unknown BuilderConfig option(s) in config file: "
                f"{', '.join(sorted(unknown))}. "
                f"Valid options are: {', '.join(sorted(valid_names))}."
            )

        normalised = dict(overrides)
        if "drop_bcids" in normalised:
            normalised["drop_bcids"] = frozenset(
                int(value) for value in normalised["drop_bcids"])
        if "max_hits_per_sca" in normalised:
            normalised["max_hits_per_sca"] = float(normalised["max_hits_per_sca"])

        return replace(cls(), **normalised)
