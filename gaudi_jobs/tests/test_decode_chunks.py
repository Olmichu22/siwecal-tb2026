"""The decode health check must actually fire.

A single-process decode of a whole run silently drops acquisitions on some runs:
run_000012 was decoded, reconstructed, plotted and believed with 75% of its data
missing. Nothing failed. The only thing that would have caught it is the ratio
between tree entries and the largest acqNumber -- one entry per acquisition means
a healthy ratio of ~1, and run_000012's was 4.

decode_chunks.assert_healthy() is that check. These tests pin it, because a guard
that does not fire is worse than no guard: it makes the data look verified.

Run:  python3 -m pytest gaudi_jobs/tests/ -q
"""
import array
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from decode_chunks import assert_healthy, decode_health  # noqa: E402

ROOT = pytest.importorskip("ROOT")


def _tree(path, entries, acq_step):
    """A siwecaldecoded-like tree: `entries` rows whose acqNumber climbs by
    `acq_step`. acq_step=1 is a healthy decode (one entry per acquisition);
    acq_step=4 is the run_000012 signature (three of every four gone)."""
    handle = ROOT.TFile(str(path), "RECREATE")
    tree = ROOT.TTree("siwecaldecoded", "")
    acq = array.array("i", [0])
    tree.Branch("acqNumber", acq, "acqNumber/I")
    for i in range(entries):
        acq[0] = (i + 1) * acq_step
        tree.Fill()
    tree.Write()
    handle.Close()
    return str(path)


def test_healthy_decode_passes(tmp_path):
    path = _tree(tmp_path / "good.root", 100, 1)
    entries, max_acq, ratio = decode_health([path])
    assert (entries, max_acq) == (100, 100)
    assert ratio == pytest.approx(1.0)
    assert assert_healthy([path], "RUN_OK") == 100


def test_lossy_decode_is_rejected(tmp_path):
    """The exact run_000012 failure: a valid ROOT file with a quarter of the data."""
    path = _tree(tmp_path / "lossy.root", 100, 4)
    assert decode_health([path])[2] == pytest.approx(4.0)
    with pytest.raises(SystemExit) as excinfo:
        assert_healthy([path], "RUN_LOSSY")
    message = str(excinfo.value)
    assert "lossy" in message
    assert "75%" in message          # 1 - 1/4 of the acquisitions missing


def test_empty_decode_is_rejected(tmp_path):
    path = _tree(tmp_path / "empty.root", 0, 1)
    with pytest.raises(SystemExit):
        assert_healthy([path], "RUN_EMPTY")


def test_chunks_are_chained(tmp_path):
    """The check runs over a whole run's chunks, not one file: loss shows up in
    the total, and a per-chunk view would miss it."""
    paths = [_tree(tmp_path / f"chunk_{i:04d}.root", 50, 1) for i in range(3)]
    entries, _, _ = decode_health(paths)
    assert entries == 150
