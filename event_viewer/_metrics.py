"""
Per-event physics metrics for the SiW-ECAL validation.

Single home for *all* event-level quantities — both the original ones (number of
hits, summed energy, energy-weighted layer, per-layer hits, MIP-likeness) and the
particle-discrimination set ported from ``VarsJesus/analysis_template.C``
(per-layer energy profile, tungsten-weighted energy, transverse barycenter,
Molière radius, transverse RMS, longitudinal shower variables, ...).

Every function takes plain numpy arrays for a *single event* (the per-hit
``slab``/``energy``/``x``/``y`` vectors, or a per-layer profile) and returns
scalars or fixed-length per-layer arrays. Keeping the formulas here, isolated
from I/O and plotting, makes them directly unit-testable; :class:`EventData`
just calls them in its load loop.

Conventions
-----------
* Energies are in MIP units (our ``hit_energy``).
* The tungsten weight of a hit is ``energy * W[slab] / X0`` with ``X0 = 3.5 mm``
  (tungsten radiation length), i.e. the absorber depth in radiation lengths in
  front of the hit's layer — a sampling-corrected energy (see the C++ template).
* Empty events / non-positive total weight return zeros, mirroring the template.
* Shower-shape quantities are ``NaN`` when the event is not a shower, so they are
  naturally dropped from histograms and excluded by selection cuts.
"""

from dataclasses import dataclass

import numpy as np

NAN = float("nan")


# --------------------------------------------------------------- basic sums ---
def total_energy(energy: np.ndarray) -> float:
    """Summed hit energy of the event (``sume`` = Σ E_i), in MIP units."""
    return float(energy.sum())


def hits_per_layer(slab: np.ndarray, n_layers: int) -> np.ndarray:
    """Number of hit channels in each layer (length ``n_layers``)."""
    return np.bincount(slab, minlength=n_layers)[:n_layers].astype(float)


def energy_per_layer(slab: np.ndarray, energy: np.ndarray,
                     n_layers: int) -> np.ndarray:
    """Summed hit energy in each layer (``sume_layer``, length ``n_layers``)."""
    return np.bincount(slab, weights=energy, minlength=n_layers)[:n_layers]


def zbary(slab: np.ndarray, energy: np.ndarray) -> float:
    """Energy-weighted mean layer index: ``Σ(E_i · slab_i) / Σ E_i``.

    Returns ``NaN`` if the total energy is not positive.
    """
    total = energy.sum()
    if total <= 0:
        return NAN
    return float(np.dot(energy, slab.astype(float)) / total)


def mip_likeness(hits_layer: np.ndarray, n_layers: int) -> float:
    """MIP-likeness score: mean over hit layers of ``1 / hits_in_layer``.

    ``score = (1/n_layers) · Σ_layers [ 1/n_hit(layer) ]`` over layers with at
    least one hit. A through-going MIP (≈1 hit/layer) scores near 1; a dense
    shower (many hits/layer) scores near 0.
    """
    with np.errstate(divide="ignore"):
        inverse = np.where(hits_layer > 0, 1.0 / hits_layer, 0.0)
    return float(inverse.sum() / n_layers)


def layer_extent(hits_layer: np.ndarray):
    """``(first_layer, last_layer, n_layers_hit)`` from the per-layer hit counts.

    ``first``/``last`` are the lowest/highest layer index with any hit; both are
    ``NaN`` for an empty event. ``n_layers_hit`` is the count of layers hit.
    """
    hit = np.flatnonzero(hits_layer > 0)
    if hit.size == 0:
        return NAN, NAN, 0
    return float(hit[0]), float(hit[-1]), int(hit.size)


# --------------------------------------------------- tungsten-weighted energy -
def hit_weights(slab: np.ndarray, energy: np.ndarray,
                w_over_x0: np.ndarray) -> np.ndarray:
    """Per-hit tungsten-weighted energy: ``E_i · (W[slab_i] / X0)``.

    ``w_over_x0`` is the per-layer absorber depth in radiation lengths.
    """
    return energy * w_over_x0[slab]


def weighte_total(slab: np.ndarray, energy: np.ndarray,
                  w_over_x0: np.ndarray) -> float:
    """Sampling-corrected total energy ``Σ_i E_i · W[slab_i] / X0``."""
    return float(hit_weights(slab, energy, w_over_x0).sum())


def weighte_per_layer(slab: np.ndarray, energy: np.ndarray,
                      w_over_x0: np.ndarray, n_layers: int) -> np.ndarray:
    """Tungsten-weighted energy in each layer (length ``n_layers``)."""
    weights = hit_weights(slab, energy, w_over_x0)
    return np.bincount(slab, weights=weights, minlength=n_layers)[:n_layers]


# ---------------------------------------------------------- transverse shape --
def barycenter_xy(x: np.ndarray, y: np.ndarray, weights: np.ndarray):
    """Energy-weighted transverse barycenter ``(bar_x, bar_y, bar_r)``.

    ``bar_x = Σ w_i x_i / Σ w_i`` (idem y), ``bar_r = sqrt(bar_x² + bar_y²)``.
    Returns ``(0, 0, 0)`` when the total weight is not positive.
    """
    total = weights.sum()
    if total <= 0:
        return 0.0, 0.0, 0.0
    bx = float(np.dot(weights, x) / total)
    by = float(np.dot(weights, y) / total)
    return bx, by, float(np.hypot(bx, by))


def transverse_rms(x: np.ndarray, y: np.ndarray, weights: np.ndarray,
                   bar_x: float, bar_y: float) -> float:
    """Energy-weighted RMS hit radius about the shower axis (``bar_x, bar_y``):

    ``rms = sqrt( Σ w_i · r_i² / Σ w_i )`` with ``r_i² = (x_i-bx)² + (y_i-by)²``.
    A cheap compactness proxy. Returns 0 for non-positive total weight.
    """
    total = weights.sum()
    if total <= 0:
        return 0.0
    r2 = (x - bar_x) ** 2 + (y - bar_y) ** 2
    # max(0, .) guards against tiny negatives from any residual negative weight.
    return float(np.sqrt(max(0.0, np.dot(weights, r2) / total)))


def moliere_radius(x: np.ndarray, y: np.ndarray, weights: np.ndarray,
                   bar_x: float, bar_y: float, containment: float = 0.90) -> float:
    """Transverse radius containing ``containment`` of the (weighted) energy.

    Hits are ordered by their distance ``r_i`` to the shower axis; the radius
    returned is the smallest ``r`` whose cumulative energy reaches
    ``containment · Σ w``. With ``containment = 0.9`` this is the Molière radius
    estimate of the C++ template. Returns 0 for an empty / zero-weight event.
    """
    total = weights.sum()
    if x.size == 0 or total <= 0:
        return 0.0
    r = np.hypot(x - bar_x, y - bar_y)
    order = np.argsort(r, kind="stable")
    cumulative = np.cumsum(weights[order])
    reached = np.searchsorted(cumulative, containment * total)
    reached = min(reached, r.size - 1)
    return float(r[order][reached])


def core_hits_per_layer(slab: np.ndarray, x: np.ndarray, y: np.ndarray,
                        bar_x: float, bar_y: float, n_layers: int,
                        radius_mm: float = 30.0) -> np.ndarray:
    """Per-layer hit count restricted to the shower core.

    Only hits within ``radius_mm`` of the transverse barycentre ``(bar_x,
    bar_y)`` are counted. The spatial restriction is what keeps scattered noise
    from making a layer look dense: in real data the pedestal tail puts hits all
    over the detector, and a bare occupancy count would add them into a layer
    that has no core at all.

    Pass the *positive-energy* hits only, as the barycentre itself is computed
    from those (see :func:`barycenter_xy`).
    """
    if slab.size == 0 or not np.isfinite(bar_x) or not np.isfinite(bar_y):
        return np.zeros(n_layers, dtype=int)
    r2 = (x - bar_x) ** 2 + (y - bar_y) ** 2
    inside = slab[r2 < radius_mm * radius_mm]
    return np.bincount(inside, minlength=n_layers)[:n_layers].astype(int)


# ------------------------------------------------------ longitudinal (shower) -
@dataclass(frozen=True)
class ShowerFeatures:
    """Longitudinal shower descriptors derived from a per-layer profile."""

    is_shower: bool
    start_layer: float      # first layer above threshold (NaN if not a shower)
    max_layer: float        # layer of the profile peak
    end_layer: float        # last layer above threshold
    start_layer_10: float   # first layer above threshold AND > 10% of the peak
    end_layer_10: float     # last  layer above threshold AND > 10% of the peak
    length: float           # number of layers above threshold
    max_value: float        # height of the profile peak
    onset: float            # conversion layer: first of the dense run (NaN if not)
    n_layers_before_onset: float   # hit layers ahead of it (pre-shower track)


def shower_onset(core_hits: np.ndarray, profile: np.ndarray, max_min: float,
                 min_nhit: int = 4, min_consecutive: int = 2) -> int:
    """Layer where the cascade starts, or ``-1`` if the event does not show one.

    Two conditions, both required:

    * ``peak(profile) > max_min`` — cheap pre-filter;
    * ``min_consecutive`` consecutive layers each with at least ``min_nhit``
      hits in the core (``core_hits``, from :func:`core_hits_per_layer`); the
      onset is the first of them.

    Mirror of ``showerOnset()`` in ``k4SiWEcalReco/EcalShowerVars.h``, which is
    the source of truth — it is the criterion actually written into the files.
    It is also the definition the simulation's ``ShowerTagger`` uses to decide
    which hits stay out of the ACTS measurement pool, so simulation and analysis
    agree on where a shower begins.

    This replaced an edge detector ("three consecutive rising layers"). A rising
    requirement needs the profile to still be climbing, so a cascade starting
    near the back of the stack — a pion punching through and interacting late,
    the topology PID cares most about — can run out of layers before producing
    the required rise, and be missed. A dense run is still what separates a
    cascade from a MIP: a through-going muon lights one or two pads per layer, so
    it never reaches ``min_nhit``, let alone twice in a row.
    """
    if profile.size == 0 or min_consecutive < 1:
        return -1
    if float(profile.max()) <= max_min:
        return -1
    run = 0
    for layer in range(core_hits.size):
        if core_hits[layer] >= min_nhit:
            run += 1
            if run >= min_consecutive:
                return layer - (run - 1)
        else:
            run = 0
    return -1


def shower_features(profile: np.ndarray, threshold: float, max_min: float,
                    start_frac: float = 0.1, *,
                    core_hits: np.ndarray, hits_layer: np.ndarray,
                    onset_min_nhit: int = 4,
                    onset_min_consecutive: int = 2) -> ShowerFeatures:
    """Locate the shower along the layers from a per-layer ``profile``.

    The profile is the per-layer activity (hits, energy or weighted energy). The
    event is a shower when :func:`shower_onset` finds an onset; then:

    * ``onset`` = the conversion layer, i.e. the first of the dense run,
    * ``n_layers_before_onset`` = hit layers ahead of it, counted on
      ``hits_layer`` (the raw per-layer hits, so a MIP that scatters off the
      shower axis still counts as a traversed layer). This is the length of the
      incoming track segment: 0-1 for an electron converting immediately,
      several layers for a pion or radiative muon that traverses as a MIP first.
    * ``max_layer``  = argmax of the profile,
    * ``start_layer``/``end_layer`` = first/last layer above ``threshold``,
    * ``start_layer_10``/``end_layer_10`` additionally require the layer to
      exceed ``start_frac · peak`` (10% of the maximum by default),
    * ``length`` = number of layers above ``threshold``.

    ``onset`` and ``start_layer`` answer nearly the same question but are not
    interchangeable: ``start_layer`` searches only *before* the profile maximum,
    so it is ``NaN`` whenever the maximum sits in the first layer, whereas
    ``onset`` is defined whenever the flag is set.

    ``core_hits`` and ``hits_layer`` are keyword-only and required on purpose:
    the criterion needs the hit positions, and a profile-only signature with a
    silent fallback is exactly how this function drifted away from the C++ it is
    supposed to mirror.

    For non-showers everything but ``is_shower`` is ``NaN`` (the template uses
    a ``-1`` sentinel; ``NaN`` is cleaner for histograms and cuts here).
    """
    onset = -1
    if profile.size and profile.sum() > 0:
        onset = shower_onset(core_hits, profile, max_min,
                             onset_min_nhit, onset_min_consecutive)
    if onset < 0:
        return ShowerFeatures(False, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN, NAN)

    peak = float(profile.max())
    max_layer = int(np.argmax(profile))
    above = profile > threshold
    above_10 = above & (profile > start_frac * peak)

    def _first(mask, lo, hi):
        idx = np.flatnonzero(mask[lo:hi])
        return float(idx[0] + lo) if idx.size else NAN

    def _last(mask, lo, hi):
        idx = np.flatnonzero(mask[lo:hi])
        return float(idx[-1] + lo) if idx.size else NAN

    return ShowerFeatures(
        is_shower=True,
        start_layer=_first(above, 0, max_layer),
        max_layer=float(max_layer),
        end_layer=_last(above, max_layer + 1, profile.size),
        start_layer_10=_first(above_10, 0, max_layer),
        end_layer_10=_last(above_10, max_layer + 1, profile.size),
        length=float(int(above.sum())),
        max_value=peak,
        onset=float(onset),
        n_layers_before_onset=float(int(np.count_nonzero(hits_layer[:onset] > 0))),
    )
