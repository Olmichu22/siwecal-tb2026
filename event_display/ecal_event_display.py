#!/usr/bin/env python3
"""
ecal_event_display.py

TEve-based 3D event display for the SiW-ECAL TB2026CERN prototype.
Draws the detector geometry (Si planes, W absorber plates, Al box) and
overlays reconstructed hits for a single event.

Usage:
    python3 ecal_event_display.py
    python3 ecal_event_display.py -i data/run.root -e 3
    python3 ecal_event_display.py --z-file /path/to/slab_z_positions.yml
"""

import argparse
import os
import yaml
import ROOT

# Keep TColor references alive — Python GC would destroy them otherwise.
_color_refs = []

_DEFAULT_Z_FILE = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "conversion", "slab_z_positions.yml",
)

# ---------------------------------------------------------------------------
# Geometry constants  (all in cm)
# ---------------------------------------------------------------------------
_SI_HALF_XY   = 9.0    # 180 mm active area / 2
_SI_HALF_Z    = 0.025  # 0.5 mm Si plane half-thickness
_W_GAP        = 0.5    # gap from W downstream face to Si plane centre
_AL_HALF_XY   = 9.5    # Al box is slightly wider than active area
_AL_HALF_Z    = 0.3    # 6 mm Al plate half-thickness
_AL_FRONT_OFF = -3.0    # Al front plate centre: this far upstream of slab 0
_AL_REAR_OFF  = -3.0    # Al rear plate centre: this far downstream of last slab

_HIT_HALF_XY  = 0.275  # 5.5 mm cell / 2, in cm
_HIT_HALF_Z   = _SI_HALF_Z  # same as Si plane thickness

_N_PALETTE    = 64     # number of energy colour steps
_hit_palette  = []     # populated by _init_hit_palette()


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _rgb(r, g, b):
    """Register an RGB colour (0..1 floats) and return its ROOT index."""
    idx = ROOT.TColor.GetFreeColorIndex()
    c = ROOT.TColor(idx, r, g, b, "")
    _color_refs.append(c)
    return idx


def _box(name, dx, dy, dz, x, y, z, color, transparency=0):
    """Create a TEveGeoShape axis-aligned box.  All values in cm."""
    gs = ROOT.TEveGeoShape(name)
    gs.SetShape(ROOT.TGeoBBox(dx, dy, dz))
    gs.SetMainColor(color)
    gs.SetFillColor(color)
    gs.SetLineColor(color)
    gs.SetMainTransparency(transparency)
    gs.RefMainTrans().SetPos(x, y, z)
    gs.ResetBBox()
    return gs


def _init_hit_palette():
    """Register _N_PALETTE colours from blue → green → yellow → red."""
    global _hit_palette
    if _hit_palette:
        return
    n = _N_PALETTE
    for i in range(n):
        t = i / (n - 1)
        if t < 0.33:
            s = t / 0.33
            r, g, b = 0.0, s, 1.0 - s
        elif t < 0.66:
            s = (t - 0.33) / 0.33
            r, g, b = s, 1.0, 0.0
        else:
            s = (t - 0.66) / 0.34
            r, g, b = 1.0, 1.0 - s, 0.0
        _hit_palette.append(_rgb(r, g, b))


def _energy_color(energy, e_max):
    """Map energy → palette index.  Clamps negatives to index 0."""
    if e_max <= 0:
        return _hit_palette[0]
    t = max(0.0, min(1.0, energy / e_max))
    return _hit_palette[int(t * (_N_PALETTE - 1))]


# ---------------------------------------------------------------------------
# Hit reader
# ---------------------------------------------------------------------------

def load_hits(root_file, event_idx):
    """Return list of (x_cm, y_cm, z_index, energy) for one event."""
    f = ROOT.TFile.Open(root_file)
    if not f or f.IsZombie():
        raise SystemExit(f"ERROR: cannot open {root_file}")
    tree = f.Get("ecal")
    if not tree:
        raise SystemExit("ERROR: TTree 'ecal' not found in file")
    n_entries = tree.GetEntries()
    if event_idx >= n_entries:
        raise SystemExit(f"ERROR: event {event_idx} does not exist (tree has {n_entries} entries)")

    tree.GetEntry(event_idx)
    nhit = tree.nhit_chan
    hits = []
    for i in range(nhit):
        hits.append((
            tree.hit_x[i] / 10.0,       # mm → cm
            tree.hit_y[i] / 10.0,
            int(tree.hit_slab[i]),
            float(tree.hit_energy[i]),
        ))
    print(f"[Event {event_idx}] {nhit} hits  (run={tree.run} spill={tree.spill} bcid={tree.bcid})")
    f.Close()
    return hits


# ---------------------------------------------------------------------------
# Hit builder
# ---------------------------------------------------------------------------

def build_hits(eve, hits, slab_z_cm):
    """Draw hit cells in the event scene, coloured by energy."""
    _init_hit_palette()

    energies = [e for *_, e in hits]
    e_max = max((e for e in energies if e > 0), default=1.0)
    print(f"[Hits] energy range: {min(energies):.2f} .. {e_max:.2f}")

    hit_grp = ROOT.TEveElementList("Hits")
    n_slabs = len(slab_z_cm)
    for i, (x, y, slab, energy) in enumerate(hits):
        if slab < 0 or slab >= n_slabs:
            continue
        z = slab_z_cm[slab]
        col = _energy_color(energy, e_max)
        hit_grp.AddElement(_box(
            f"hit_{i}",
            _HIT_HALF_XY, _HIT_HALF_XY, _HIT_HALF_Z,
            x, y, z,
            col, transparency=20,
        ))

    eve.GetEventScene().AddElement(hit_grp)
    print(f"[Hits] {hit_grp.NumChildren()} hit boxes added")


# ---------------------------------------------------------------------------
# Geometry builder
# ---------------------------------------------------------------------------

def build_geometry(eve, slab_z_mm, w_thickness_mm):
    col_si = _rgb(0.20, 0.80, 0.30)   # green
    col_w  = _rgb(0.50, 0.30, 0.60)   # purple
    col_al = _rgb(0.65, 0.65, 0.70)   # silver-grey

    geo = ROOT.TEveElementList("Geometry")

    # --- Si planes ---
    si_grp = ROOT.TEveElementList("Si planes")
    for i, z_mm in enumerate(slab_z_mm):
        z_cm = z_mm / 10.0
        si_grp.AddElement(_box(
            f"Si_{i}",
            _SI_HALF_XY, _SI_HALF_XY, _SI_HALF_Z,
            0.0, 0.0, z_cm,
            col_si, transparency=80,
        ))
    geo.AddElement(si_grp)
    print(f"[Geometry] {len(slab_z_mm)} Si planes")

    # --- W absorber plates ---
    w_grp = ROOT.TEveElementList("W plates")
    for i, (z_mm, w_mm) in enumerate(zip(slab_z_mm, w_thickness_mm)):
        z_cm   = z_mm / 10.0
        w_half = w_mm / 20.0                        # mm → cm → half
        w_z    = z_cm - _W_GAP - w_half             # centre of W plate
        w_grp.AddElement(_box(
            f"W_{i}",
            _SI_HALF_XY, _SI_HALF_XY, w_half,
            0.0, 0.0, w_z,
            col_w, transparency=60,
        ))
    geo.AddElement(w_grp)
    print(f"[Geometry] {len(w_thickness_mm)} W plates")

    # --- Al box (front and rear plates) ---
    first_z = slab_z_mm[0]  / 10.0
    last_z  = slab_z_mm[-1] / 10.0
    al_grp  = ROOT.TEveElementList("Al box")
    for name, z_cm in (("Al_front", first_z - _AL_FRONT_OFF),
                        ("Al_rear",  last_z  + _AL_REAR_OFF)):
        al_grp.AddElement(_box(
            name,
            _AL_HALF_XY, _AL_HALF_XY, _AL_HALF_Z,
            0.0, 0.0, z_cm,
            col_al, transparency=70,
        ))
    geo.AddElement(al_grp)
    print(f"[Geometry] Al box: front z={first_z - _AL_FRONT_OFF:.1f} cm, "
          f"rear z={last_z + _AL_REAR_OFF:.1f} cm")

    eve.GetGlobalScene().AddElement(geo)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="SiW-ECAL TEve event display (TB2026CERN)",
    )
    parser.add_argument(
        "--z-file", default=_DEFAULT_Z_FILE,
        help=f"slab_z_positions.yml  (default: {_DEFAULT_Z_FILE})",
    )
    parser.add_argument(
        "-i", "--input", default=None,
        help="ROOT file with the 'ecal' TTree (omit to show geometry only)",
    )
    parser.add_argument(
        "-e", "--event", type=int, default=0,
        help="Event index to display (default: 0)",
    )
    args = parser.parse_args(argv)

    if not os.path.exists(args.z_file):
        raise SystemExit(f"ERROR: z-file not found: {args.z_file}")
    with open(args.z_file) as fh:
        doc = yaml.safe_load(fh)
    slab_z_mm      = [float(v) for v in doc["slab_z_mm"]]
    w_thickness_mm = [float(v) for v in doc["w_thickness_mm"]]
    slab_z_cm      = [z / 10.0 for z in slab_z_mm]
    print(f"[Config] {len(slab_z_mm)} slabs from {args.z_file}")

    ROOT.gErrorIgnoreLevel = ROOT.kWarning
    ROOT.gSystem.Load("libEve")
    ROOT.gSystem.Load("libEG")

    eve = ROOT.TEveManager.Create(True, "FI")

    # Explicitly wire both scenes to a viewer — without this the GL
    # canvas stays grey (scenes exist but are not attached to any viewer).
    viewer = eve.SpawnNewViewer("3D View", "")
    viewer.AddScene(eve.GetGlobalScene())
    viewer.AddScene(eve.GetEventScene())

    # Let the initial GL draw pass complete before modifying scenes;
    # without this the GL scenes are DrawLocked when we call AddElement.
    ROOT.gSystem.ProcessEvents()

    build_geometry(eve, slab_z_mm, w_thickness_mm)

    if args.input is not None:
        hits = load_hits(args.input, args.event)
        build_hits(eve, hits, slab_z_cm)

    eve.Redraw3D(True)
    ROOT.gSystem.ProcessEvents()

    # Quit the application when the TEve browser window is closed; without
    # this the GUI event loop keeps the Python process alive forever.
    browser = eve.GetBrowser()
    browser.Connect("CloseWindow()", "TApplication", ROOT.gApplication,
                    "Terminate()")

    # Restore default signal handling so Ctrl-C in the terminal also exits
    # (ROOT installs its own SIGINT handler that otherwise swallows it).
    ROOT.gSystem.ResetSignals()

    print("[Ready] Close the TEve window (or press Ctrl-C) to exit.")
    ROOT.gApplication.Run(True)


if __name__ == "__main__":
    main()
