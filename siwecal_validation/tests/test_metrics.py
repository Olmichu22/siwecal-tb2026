"""
Unit tests for :mod:`siwecal_validation.metrics`.

Synthetic single events with known answers. Runs under pytest, or standalone:

    python3 siwecal_validation/tests/test_metrics.py
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

from siwecal_validation import metrics as m  # noqa: E402


def test_basic_sums():
    slab = np.array([0, 0, 2], dtype=np.int32)
    energy = np.array([1.0, 2.0, 3.0], dtype=float)
    assert m.total_energy(energy) == 6.0
    hpl = m.hits_per_layer(slab, 4)
    assert list(hpl) == [2, 0, 1, 0]
    epl = m.energy_per_layer(slab, energy, 4)
    assert list(epl) == [3.0, 0.0, 3.0, 0.0]
    # zbary = (1*0 + 2*0 + 3*2) / 6 = 1.0
    assert m.zbary(slab, energy) == 1.0


def test_zbary_empty():
    assert math.isnan(m.zbary(np.array([], dtype=np.int32), np.array([])))


def test_mip_likeness():
    # layer 0 has 2 hits, layer 2 has 1 hit -> (1/2 + 1/1)/n_layers
    hpl = np.array([2.0, 0.0, 1.0, 0.0])
    assert m.mip_likeness(hpl, 4) == (0.5 + 1.0) / 4
    # a perfect MIP: 1 hit in each of n layers -> score 1.0
    assert m.mip_likeness(np.ones(5), 5) == 1.0


def test_layer_extent():
    hpl = np.array([0.0, 3.0, 0.0, 2.0, 0.0])
    first, last, n = m.layer_extent(hpl)
    assert (first, last, n) == (1.0, 3.0, 2)
    f, l, n0 = m.layer_extent(np.zeros(5))
    assert math.isnan(f) and math.isnan(l) and n0 == 0


def test_weighte():
    slab = np.array([0, 8], dtype=np.int32)        # W/X0 = 0.8 and 1.6
    energy = np.array([1.0, 1.0], dtype=float)
    w_over_x0 = np.array([0.8] + [1.2] * 7 + [1.6] * 7)
    # 1*0.8 + 1*1.6 = 2.4
    assert abs(m.weighte_total(slab, energy, w_over_x0) - 2.4) < 1e-9
    wpl = m.weighte_per_layer(slab, energy, w_over_x0, 15)
    assert abs(wpl[0] - 0.8) < 1e-9 and abs(wpl[8] - 1.6) < 1e-9


def test_barycenter_and_rms():
    # four unit-weight hits at (+-1, 0), (0, +-1): barycenter at origin, r=1
    x = np.array([1.0, -1.0, 0.0, 0.0])
    y = np.array([0.0, 0.0, 1.0, -1.0])
    w = np.ones(4)
    bx, by, br = m.barycenter_xy(x, y, w)
    assert abs(bx) < 1e-12 and abs(by) < 1e-12 and abs(br) < 1e-12
    assert abs(m.transverse_rms(x, y, w, bx, by) - 1.0) < 1e-12


def test_barycenter_zero_weight():
    x = np.array([1.0]); y = np.array([1.0]); w = np.array([0.0])
    assert m.barycenter_xy(x, y, w) == (0.0, 0.0, 0.0)


def test_moliere_ring():
    # ring of equal-energy hits at radius 5 -> 90% radius is 5; single hit -> 0
    n = 20
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    x = 5 * np.cos(ang); y = 5 * np.sin(ang)
    w = np.ones(n)
    assert abs(m.moliere_radius(x, y, w, 0.0, 0.0, 0.9) - 5.0) < 1e-9
    assert m.moliere_radius(np.array([0.0]), np.array([0.0]), np.array([1.0]),
                            0.0, 0.0) == 0.0


def test_moliere_two_radii():
    # 90 hits at r=1, 10 hits at r=10; 90% containment -> radius 1
    x = np.concatenate([np.ones(90), 10 * np.ones(10)])
    y = np.zeros(100)
    w = np.ones(100)
    assert abs(m.moliere_radius(x, y, w, 0.0, 0.0, 0.90) - 1.0) < 1e-9


def test_shower_features_shower():
    # rising profile peaking at layer 5, then falling
    profile = np.array([0, 1, 4, 8, 14, 20, 12, 6, 2, 0,
                        0, 0, 0, 0, 0], dtype=float)
    f = m.shower_features(profile, threshold=5.0, max_min=10.0, start_frac=0.1)
    assert f.is_shower is True
    assert f.max_layer == 5.0 and f.max_value == 20.0
    assert f.start_layer == 3.0           # first layer > 5 before the peak
    assert f.end_layer == 7.0             # last layer > 5 after the peak
    assert f.length == 5.0                # layers 3,4,5,6,7 are > 5


def test_shower_features_mip():
    # flat MIP-like profile (a couple of hits per layer): not a shower
    profile = np.full(15, 2.0)
    f = m.shower_features(profile, threshold=5.0, max_min=10.0)
    assert f.is_shower is False
    assert math.isnan(f.start_layer) and math.isnan(f.max_layer)


def _run_standalone():
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for fn in fns:
        fn()
        print(f"  ok  {fn.__name__}")
    print(f"\nAll {len(fns)} metric tests passed.")


if __name__ == "__main__":
    _run_standalone()
