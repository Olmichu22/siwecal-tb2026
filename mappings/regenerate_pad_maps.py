#!/usr/bin/env python3
"""Rewrite the pad maps onto the measured sensor grid from the spec sheet.

The maps are the ``(chip, channel) -> (x, y)`` tables the event builder uses to
turn a raw hit's ``(chip, channel)`` into a transverse position. They were
written on the *nominal* design grid: a round 5.5 mm pitch, 7.6 mm across the
boundary between two sensors. The sensor spec sheet gives the real numbers -- a
5.52 mm diode on a 5.53 mm pixel pitch, with a 0.61 mm inactive rim and no guard
ring -- which put the pad centres up to 0.43 mm away from where the nominal
table has them. Small, but it is a real position error on every reconstructed
hit, so the map should describe the sensor as built, not as designed.

Unlike the simulation project, there is no DD4hep compact file here to read the
grid from: this is real test-beam data. The target grid is therefore built
straight from the spec-sheet constants below.

This script maps the old columns onto the new ones by rank -- the i-th sorted
column of the old map becomes the i-th column of the new grid -- and rewrites
the files in place. Only coordinates change: the chip/channel assignment, the
row order and the header are copied through untouched, so anything keyed on
(chip, channel) -- masking lists, MIP calibration -- is unaffected.

IMPORTANT: this changes where every real hit lands. Runs already built with the
nominal map keep their old positions until they are re-run through the event
builder; regenerating the map only affects builds done afterwards.

Usage:
    python mappings/regenerate_pad_maps.py            # rewrite in place
    python mappings/regenerate_pad_maps.py --dry-run  # just report the shift
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys

import numpy as np

MAPPINGS = os.path.abspath(os.path.dirname(__file__))

MAP_FILES = [
    "fev10_rotate_chip_channel_x_y_mapping.txt",
    "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt",
]

# --- Sensor geometry, from the spec sheet (all in mm) ------------------------
PAD_PITCH = 5.53          # pixel pitch: 5.52 mm diode + 0.01 mm inter-pad margin
N_PADS_PER_WAFER = 16     # per side
WAFER_MARGIN = 0.61       # inactive rim of a sensor, per side (no guard ring)
N_WAFERS = 2              # per side (2x2 sensors per layer)
WAFER_GAP = 0.0           # extra clearance between two butted sensors

WAFER_ACTIVE = N_PADS_PER_WAFER * PAD_PITCH          # 88.48 mm sensitive area
WAFER_SIZE = WAFER_ACTIVE + 2 * WAFER_MARGIN         # 89.70 mm physical sensor


def geometry_columns() -> np.ndarray:
    """Pad centres of one layer along one axis [mm], from the spec sheet.

    A pad's centre is its position inside its sensor's pad array, plus the
    position of that sensor. The array is centred on the sensor, so the rim
    shifts nothing; it only sets how far apart two sensors sit.
    """
    columns = []
    for s in range(N_WAFERS):
        centre = (s - (N_WAFERS - 1) / 2.0) * (WAFER_SIZE + WAFER_GAP)
        for p in range(N_PADS_PER_WAFER):
            columns.append(centre - WAFER_ACTIVE / 2.0 + (p + 0.5) * PAD_PITCH)
    return np.asarray(sorted(columns))


def read_rows(path: str):
    """``(header, [(chip, x0, y0, channel, x, y)])`` of a pad map file."""
    header, rows = None, []
    for line in open(path):
        parts = line.split()
        if not parts:
            continue
        try:
            chip = int(parts[0])
        except ValueError:
            if header is None:
                header = line.rstrip("\n")
            continue
        rows.append((chip, float(parts[1]), float(parts[2]),
                     int(parts[3]), float(parts[4]), float(parts[5])))
    return header, rows


def fmt(value: float) -> str:
    """Match the style of the existing files: plain decimals, no trailing zeros."""
    text = f"{value:.4f}".rstrip("0").rstrip(".")
    return "0" if text in ("", "-0") else text


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--dry-run", action="store_true",
                    help="report the shift without touching the files")
    ap.add_argument("--no-backup", action="store_true",
                    help="do not leave a .nominal copy of the old file")
    args = ap.parse_args()

    new_cols = geometry_columns()
    print(f"spec-sheet grid: {new_cols.size} columns, "
          f"{new_cols[0]:.4f} .. {new_cols[-1]:.4f} mm "
          f"(pitch {PAD_PITCH}, boundary {PAD_PITCH + 2 * WAFER_MARGIN:.2f} mm)")

    for name in MAP_FILES:
        path = os.path.join(MAPPINGS, name)
        header, rows = read_rows(path)
        old_cols = np.asarray(sorted({round(r[4], 4) for r in rows}
                                     | {round(r[5], 4) for r in rows}))
        if old_cols.size != new_cols.size:
            sys.exit(f"{name}: map has {old_cols.size} columns, "
                     f"spec-sheet grid has {new_cols.size}")

        def remap(v: float) -> float:
            return float(np.interp(v, old_cols, new_cols))

        shift = max(abs(remap(c) - c) for c in old_cols)
        print(f"{name}: {len(rows)} pads, largest move {shift:.4f} mm")
        if args.dry_run:
            continue

        if not args.no_backup:
            shutil.copyfile(path, path + ".nominal")
        with open(path, "w") as out:
            if header:
                out.write(header + "\n")
            for chip, x0, y0, channel, x, y in rows:
                out.write(f"{chip} {fmt(remap(x0))} {fmt(remap(y0))} "
                          f"{channel} {fmt(remap(x))} {fmt(remap(y))}\n")

    if args.dry_run:
        print("\ndry run, nothing written")


if __name__ == "__main__":
    main()
