"""The tungsten stack is defined in several places. They must all agree.

In July 2026 they did not: mappings/Tungsten_thickness.yml (and, copied from it,
siwecal_validation) put the first 5.6 mm plate in front of slab 8, while the C++
reconstruction had it in front of slab 9. Nothing caught it, because the two
wrong copies were the ones no result depended on. These tests would have.

The stack: 2.8 mm of W in front of slab 0, 4.2 mm for slabs 1-8, 5.6 mm for
slabs 9-14 -- 70.0 mm, 20.00 X0 with X0(W) = 3.5 mm.
"""
import os
import re

import pytest
import yaml

from siwecal_validation import config

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))

EXPECTED_MM = (2.8,) + (4.2,) * 8 + (5.6,) * 6
EXPECTED_X0_MM = 3.5


def _read_thicknesses_from_source(path, pattern):
    """Pull a flat list of floats out of a literal in a source file."""
    with open(os.path.join(_REPO, path)) as handle:
        match = re.search(pattern, handle.read(), re.S)
    assert match, f"{path}: could not find the thickness literal ({pattern})"
    return tuple(float(tok) for tok in match.group(1).replace("\n", " ").split(",") if tok.strip())


def test_validation_default_is_the_real_stack():
    assert config.W_THICKNESSES_DEFAULT == EXPECTED_MM
    assert config.W_X0_MM == EXPECTED_X0_MM


def test_tungsten_yaml_resolves_to_the_real_stack():
    """mappings/Tungsten_thickness.yml, walked in beam order."""
    resolved = config.load_w_thicknesses(os.path.join(_REPO, "mappings", "Tungsten_thickness.yml"))
    assert resolved == EXPECTED_MM


def test_slab8_sits_behind_4_2_mm():
    """The exact bug, pinned: slab 8 is behind 4.2 mm, not 5.6."""
    resolved = config.load_w_thicknesses(os.path.join(_REPO, "mappings", "Tungsten_thickness.yml"))
    assert resolved[8] == pytest.approx(4.2)
    assert resolved[9] == pytest.approx(5.6)
    assert sum(resolved) == pytest.approx(70.0)
    assert sum(resolved) / EXPECTED_X0_MM == pytest.approx(20.0)


def test_cpp_geometry_defaults_agree():
    """gaudi_source's compiled-in geometry -- what the event builder uses for
    hit_w_energy and hit_X0, and (since the tungsten geometry was unified) what
    EcalPidTransformer derives `weighte` from."""
    cpp = _read_thicknesses_from_source(
        "gaudi_source/include/k4SiWEcalReco/PadMapGeometry.h",
        r"kDefaultSlabWThicknessMm\s*=\s*\{([^}]*)\}")
    assert cpp == EXPECTED_MM

    with open(os.path.join(_REPO, "gaudi_source/include/k4SiWEcalReco/PadMapGeometry.h")) as handle:
        x0 = re.search(r"kWX0Mm\s*=\s*([\d.]+)", handle.read())
    assert x0 and float(x0.group(1)) == EXPECTED_X0_MM


def test_slab_z_yaml_agrees():
    """event_display/conversion/slab_z_positions.yml -- the file the event
    builder is actually pointed at in production (EVBLD_SLAB_Z_FILE)."""
    with open(os.path.join(_REPO, "event_display", "conversion", "slab_z_positions.yml")) as handle:
        doc = yaml.safe_load(handle)
    assert tuple(float(v) for v in doc["w_thickness_mm"]) == EXPECTED_MM


def test_event_viewer_default_agrees():
    viewer = _read_thicknesses_from_source(
        "event_viewer/_geometry.py", r"DEFAULT_SLAB_W_THICKNESS_MM\s*=\s*\(([^)]*)\)")
    assert viewer == EXPECTED_MM
