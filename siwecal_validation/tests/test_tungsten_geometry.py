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


# --------------------------------------------------------------------------- #
# The per-slab technology block of slab_z_positions.yml
# --------------------------------------------------------------------------- #

def _slab_yaml():
    return os.path.join(_REPO, "mappings", "slab_z_positions.yml")


def test_slab_12_is_the_chip_on_board():
    """Slab 12 is the FEV11 COB: its own pad map, its own threshold DAC (243,
    run book), read from the YAML and not from a hard-coded 12 anywhere."""
    from siwecal_eventbuilder.geometry import load_slab_technology

    t = load_slab_technology(_slab_yaml())
    assert t.has_technology_block
    assert t.technology[12] == "FEV11_COB"
    assert all(v == "FEV10" for s, v in enumerate(t.technology) if s != 12)
    assert list(t.pad_map_overrides()) == [12]
    assert os.path.basename(t.pad_map_overrides()[12]).startswith("fev11_cob")
    assert os.path.isfile(t.pad_map_overrides()[12])
    assert t.pad_map_overrides_env().startswith("12:")
    assert t.slabs_with_threshold_override() == {12: 243}
    assert t.sensor_thickness_um[14] == 650


def test_cli_and_viewer_pad_maps_come_from_the_yaml():
    from siwecal_eventbuilder.cli import PAD_MAP_FILES_DEFAULT
    from siwecal_eventbuilder.geometry import load_slab_technology

    files = load_slab_technology(_slab_yaml()).pad_map_files()
    assert PAD_MAP_FILES_DEFAULT == {k: os.path.basename(v) for k, v in files.items()}
    assert set(PAD_MAP_FILES_DEFAULT) == {"default", 12}


def test_conversion_copy_of_the_slab_yaml_is_identical():
    """event_display/conversion/slab_z_positions.yml is a copy the reco.sh jobs
    read; it must not drift from mappings/."""
    twin = os.path.join(_REPO, "event_display", "conversion", "slab_z_positions.yml")
    assert open(twin).read() == open(_slab_yaml()).read()


def test_cpp_parser_reads_the_technology_block_without_corrupting_the_lists():
    """SlabGeometry::fromYamlFile used to append any '- value' line to the last
    list it had opened, so the technology block would have landed in
    w_thickness_mm.  The source must reset on every top-level key and know the
    three new ones."""
    src = open(os.path.join(_REPO, "gaudi_source", "include", "k4SiWEcalReco", "PadMapGeometry.h")).read()
    for key in ("slab_technology:", "sensor_thickness_um:", "threshold_dac:"):
        assert key in src, key
    assert "thresholdDacOverride" in src and "sensorThicknessUm" in src and "technology(" in src
    assert "current = nullptr;" in src
