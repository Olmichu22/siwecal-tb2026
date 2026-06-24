#!/usr/bin/env python3
"""
add_xy_branches.py

Reads an ecal ROOT file (TTree 'ecal') produced by the SiW-ECAL event builder
and writes a new ROOT file that is identical except for three additional branches:

    hit_x[nhit_chan]/F   -- transverse x position in mm
    hit_y[nhit_chan]/F   -- transverse y position in mm
    hit_z[nhit_chan]/F   -- longitudinal z position in mm

x/y come from the pad geometry mapping files:

    fev10_chip_channel_x_y_mapping.txt        (default, slabs 0-11, 13-14)
    fev11_cob_rotate_chip_channel_x_y_mapping.txt  (slab 12 override)

z comes from slab_z_positions.yml (one entry per slab, same directory as this
script), indexed by hit_slab.  Edit that file to match your detector survey.

Mapping files live two levels above this script (i.e. in the project root),
but can be overridden with --map-dir.

Usage
-----
    python3 add_xy_branches.py ecal_P5_74GeV_runs_000047.root
    python3 add_xy_branches.py ecal_P5_74GeV_runs_000047.root -o ecal_P5_74GeV_xyz.root
    python3 add_xy_branches.py ecal_P5_74GeV_runs_000047.root --map-dir /path/to/maps/
    python3 add_xy_branches.py ecal_P5_74GeV_runs_000047.root --z-file /path/to/slab_z_positions.yml
"""

import argparse
import array
import os
import sys

import ROOT
import yaml

ROOT.gROOT.SetBatch(True)

# ---------------------------------------------------------------------------
# Mapping configuration  (mirrors siwecal_eventbuilder/cli.py PAD_MAP_FILES_DEFAULT)
# ---------------------------------------------------------------------------
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Script is at  event_display/conversion/  -> project root is two levels up
_DEFAULT_MAP_DIR = os.path.normpath(os.path.join(_SCRIPT_DIR, "..", ".."))

# Rotated templates (standard maps with x0,y0,x,y negated = 180 deg about z).
# Fixes the hit-position bug: the un-rotated maps placed pads/hits in the wrong
# detector orientation (chip 0 in the bottom-left instead of its real top-right),
# so hit_x/hit_y came out mirrored through the origin. Keep in sync with
# siwecal_eventbuilder/cli.py PAD_MAP_FILES_DEFAULT.
MAP_FILES = {
    "default": "fev10_rotate_chip_channel_x_y_mapping.txt",
    12: "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt",
}

TREE_NAME = "ecal"
NAN_SENTINEL = -999.0   # written for hits with no mapping entry
MAX_HITS = 20000        # upper bound on nhit_chan per event

_DEFAULT_Z_FILE = os.path.join(_SCRIPT_DIR, "slab_z_positions.yml")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_mapping(path: str) -> dict:
    """Parse 'chip x0 y0 channel x y' text file.

    Returns {(chip, channel): (x, y)}.  Header and blank/comment lines are
    skipped automatically.
    """
    mapping = {}
    with open(path) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) < 6:
                continue
            try:
                chip = int(parts[0])
                channel = int(parts[3])
                x = float(parts[4])
                y = float(parts[5])
            except ValueError:
                continue  # header row or garbled line
            mapping[(chip, channel)] = (x, y)
    return mapping


def _build_pad_maps(map_dir: str) -> dict:
    pad_maps = {}
    for key, fname in MAP_FILES.items():
        path = os.path.join(map_dir, fname)
        if not os.path.exists(path):
            raise FileNotFoundError(f"Mapping file not found: {path}")
        pad_maps[key] = _load_mapping(path)
        print(f"  [map] {fname}: {len(pad_maps[key])} (chip,channel) entries")
    return pad_maps


def _get_xy(pad_maps: dict, slab: int, chip: int, channel: int):
    """Return (x, y) in mm, or (NAN_SENTINEL, NAN_SENTINEL) if unmapped."""
    m = pad_maps.get(slab, pad_maps["default"])
    xy = m.get((chip, channel))
    return xy if xy is not None else (NAN_SENTINEL, NAN_SENTINEL)


def _load_z_positions(z_file: str) -> list:
    """Load slab z positions (mm) from a YAML file.

    Returns a list indexed by slab number.  The file must contain a top-level
    key ``slab_z_mm`` whose value is a sequence with one entry per slab.
    """
    with open(z_file) as fh:
        doc = yaml.safe_load(fh)
    z_list = doc.get("slab_z_mm")
    if z_list is None:
        raise ValueError(f"Key 'slab_z_mm' not found in {z_file}")
    print(f"  [z]   {z_file}: {len(z_list)} slab positions")
    return [float(v) for v in z_list]


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def add_xyz_branches(input_path: str, output_path: str,
                     map_dir: str, z_file: str) -> None:
    in_file = ROOT.TFile.Open(input_path, "READ")
    if not in_file or in_file.IsZombie():
        raise RuntimeError(f"Cannot open input file: {input_path}")

    in_tree = in_file.Get(TREE_NAME)
    if not in_tree:
        raise RuntimeError(f"Tree '{TREE_NAME}' not found in {input_path}")

    n_entries = in_tree.GetEntries()
    print(f"  Input : {input_path}  ({n_entries} events)")

    pad_maps = _build_pad_maps(map_dir)
    z_positions = _load_z_positions(z_file)

    # ------------------------------------------------------------------
    # Create output file and clone the input tree structure (0 entries).
    # CloneTree(0) shares the leaf objects between in_tree and out_tree:
    # reading in_tree.GetEntry(i) automatically populates out_tree's
    # existing branches, so a single out_tree.Fill() copies the event.
    # ------------------------------------------------------------------
    out_file = ROOT.TFile(output_path, "RECREATE")
    out_tree = in_tree.CloneTree(0)

    # Allocate C-compatible float arrays for the new branches.
    hit_x = array.array("f", [NAN_SENTINEL] * MAX_HITS)
    hit_y = array.array("f", [NAN_SENTINEL] * MAX_HITS)
    hit_z = array.array("f", [NAN_SENTINEL] * MAX_HITS)

    # nhit_chan is already in the cloned tree; ROOT uses it as the
    # array-length counter for the variable-size branch declarations below.
    out_tree.Branch("hit_x", hit_x, "hit_x[nhit_chan]/F")
    out_tree.Branch("hit_y", hit_y, "hit_y[nhit_chan]/F")
    out_tree.Branch("hit_z", hit_z, "hit_z[nhit_chan]/F")

    n_unmapped = 0
    n_slabs = len(z_positions)

    for i in range(n_entries):
        in_tree.GetEntry(i)
        n = int(in_tree.nhit_chan)

        slabs = list(in_tree.hit_slab)[:n]
        chips = list(in_tree.hit_chip)[:n]
        chans = list(in_tree.hit_chan)[:n]

        for j in range(n):
            slab = slabs[j]
            x, y = _get_xy(pad_maps, slab, chips[j], chans[j])
            hit_x[j] = x
            hit_y[j] = y
            hit_z[j] = z_positions[slab] if slab < n_slabs else NAN_SENTINEL
            if x == NAN_SENTINEL:
                n_unmapped += 1

        out_tree.Fill()

        if (i + 1) % 10000 == 0:
            print(f"  {i + 1}/{n_entries} events processed", flush=True)

    out_file.Write("", ROOT.TObject.kOverwrite)
    out_file.Close()
    in_file.Close()

    if n_unmapped:
        print(f"  WARNING: {n_unmapped} hits had no (chip,channel) mapping entry "
              f"(written as {NAN_SENTINEL})")
    print(f"  Output: {output_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Add hit_x / hit_y / hit_z position branches to an ecal ROOT file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  %(prog)s ecal_P5_74GeV_runs_000047.root\n"
            "  %(prog)s ecal_P5_74GeV_runs_000047.root -o ecal_P5_74GeV_xyz.root\n"
        ),
    )
    parser.add_argument("input", help="Path to the input ecal ROOT file")
    parser.add_argument(
        "-o", "--output", default=None,
        help="Output path (default: <input stem>_xyz.root next to the input file)",
    )
    parser.add_argument(
        "--map-dir", default=None,
        help=f"Directory containing the .txt mapping files "
             f"(default: {_DEFAULT_MAP_DIR})",
    )
    parser.add_argument(
        "--z-file", default=None,
        help=f"YAML file with slab z positions "
             f"(default: {_DEFAULT_Z_FILE})",
    )
    return parser.parse_args(argv)


def main(argv=None) -> None:
    args = _parse_args(argv)

    map_dir = args.map_dir or _DEFAULT_MAP_DIR
    if not os.path.isdir(map_dir):
        sys.exit(f"ERROR: mapping directory not found: {map_dir}")

    z_file = args.z_file or _DEFAULT_Z_FILE
    if not os.path.exists(z_file):
        sys.exit(f"ERROR: z-positions file not found: {z_file}")

    if args.output:
        output_path = args.output
    else:
        base, ext = os.path.splitext(args.input)
        output_path = base + "_xyz" + (ext or ".root")

    print("[add_xyz_branches]")
    print(f"  input   : {args.input}")
    print(f"  output  : {output_path}")
    print(f"  map_dir : {map_dir}")
    print(f"  z_file  : {z_file}")

    add_xyz_branches(args.input, output_path, map_dir, z_file)
    print("[Done]")


if __name__ == "__main__":
    main()
