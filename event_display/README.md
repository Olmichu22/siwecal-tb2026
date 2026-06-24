# SiW-ECAL Event Display — TB2026CERN

3D event display for the SiW-ECAL prototype beam test, built on ROOT TEve.
Hits are drawn as energy-coloured pad boxes (blue → red), with detector planes
and W absorbers as transparent overlays from the geometry config. A **single
event** is rendered per launch; the TEve window is interactive (rotate / zoom /
pan the camera). To view a different event, relaunch with `-e N`.

---

## Prerequisites

- ROOT 6 with PyROOT available (`python3 -c "import ROOT"` must succeed)
- `PyYAML` (`pip install pyyaml` or available through the environment)
- A `.valcache.root` file produced by the validation pipeline
- An **X11 display** — the TEve window is a GUI. Over SSH connect with
  `ssh -Y` (or run on a machine with a local display); a headless session
  with an empty `$DISPLAY` cannot open the window.

The software dependencies are satisfied by sourcing `setup.sh` in the
repository root.

---

## Quick start

No conversion step is needed. Point the display directly at a `.valcache.root`
file with `-i` — z positions are assigned per hit from `slab_z_positions.yml`
using the `hit_slab` branch.

Use the launcher (file as first argument, optional start event as second):

```bash
cd event_display/
./launch.sh ../data/TB2026CERN_run_000013/ecal_TB2026CERN_run_000013.valcache.root
./launch.sh ../data/TB2026CERN_run_000013/ecal_TB2026CERN_run_000013.valcache.root 42
```

Or call the script directly. The ROOT file is passed with `-i`, and the event
index with `-e`:

```bash
python3 ecal_event_display.py -i ../data/TB2026CERN_run_000013/ecal_TB2026CERN_run_000013.valcache.root -e 42
```

To show the detector geometry only (no hits), omit `-i`:

```bash
python3 ecal_event_display.py
```

Each launch shows one event. Close the TEve window (or press Ctrl-C) to exit;
relaunch with a different `-e` (or `launch.sh ... N`) to view another event.

---

## Command-line options (`ecal_event_display.py`)

```
ecal_event_display.py  [-i INPUT]  [-e N]  [--z-file FILE]
```

| Option | Default | Description |
|--------|---------|-------------|
| `-i / --input INPUT` | *(none — geometry only)* | `.valcache.root` file (must have `hit_x`, `hit_y`, `hit_slab`, `hit_energy`) |
| `-e / --event N` | `0` | Event index to display |
| `--z-file FILE` | `conversion/slab_z_positions.yml` | Slab z-position and W absorber geometry config |

---

## Geometry configuration (`conversion/slab_z_positions.yml`)

Controls the z position of each slab and the thickness of its upstream
W absorber. Edit this file to match the actual survey for your run.

```yaml
slab_z_mm:
  - 0.0     # slab 0
  - 11.0    # slab 1
  # ... 15 entries total (slabs 0–14)

w_thickness_mm:
  - 2.8     # slab 0  (2.8 mm)
  - 4.2     # slab 1  (4.2 mm, slabs 1–8)
  # ... 15 entries total
```

Z is assigned to each hit as `slab_z_mm[hit_slab]`.
The W absorber for slab *i* is drawn with its downstream face 5 mm upstream
of the Si plane, i.e. centred at `slab_z[i] - 5 - w_thickness[i]/2` (mm).

---

## Coordinate orientation (why it looks rotated 90° vs. the 4×4 plot)

The display reads the **same** `hit_x` / `hit_y` branches as the per-slab
matplotlib scripts (`event_xy_4x4.py`, `heatmap_xy_4x4.py`). It does **not**
transpose, swap or negate them — see
[`load_hits()`](ecal_event_display.py#L118-L124): each hit is placed at world
position `(hit_x, hit_y, slab_z)`, hit_x → world X, hit_y → world Y, slab → Z.
So the data is identical in both tools; there is no mapping bug here.

The apparent 90° rotation is a **viewing convention**, not a data difference:

- **`event_xy_4x4.py`** uses matplotlib with a *fixed 2D frame*: X grows to the
  right, Y grows upward — the standard math orientation.
- **This display** builds a **TEve 3D GL scene** and attaches it to a freely
  rotatable GL viewer
  ([`SpawnNewViewer("3D View")`](ecal_event_display.py#L257-L259)). There is no
  fixed 2D projection. What you see is whatever orientation the **default GL
  camera** happens to give — and ROOT's default perspective camera does not
  align world-X with screen-right and world-Y with screen-up. The detector is
  long in Z (15 slabs), so the camera frames the Z axis prominently and the
  X/Y plane ends up tilted/rotated relative to matplotlib's flat view.

Because the camera is interactive, the orientation is **not locked**: rotate
with the mouse (or pan/zoom) until X is horizontal and Y vertical, looking
straight down the beam (Z) axis, and the two views match.

Two terms worth keeping distinct:

| Effect | Transform | Meaning |
|--------|-----------|---------|
| **Reflection / transpose** | swap `x ↔ y` | a real pad-mapping bug — what `debug_pad_reflection*.py` and `update_hit_xy.py` exist to catch/fix |
| **90° rotation between views** | `(x, y) → (−y, x)` | *not* a bug — purely the 3D camera vs. matplotlib axis convention |

If you want the matplotlib PNG to match the on-screen camera instead, apply a
true rotation in `event_xy_4x4.py` (`scatter(-hy[m], hx[m], ...)`), not an
axis swap.

---

## Hit colour scale

Hits are coloured by energy (MIP units) on a linear blue → red scale
recomputed per event:

| Colour | Energy |
|--------|--------|
| Blue | minimum hit energy in event |
| Cyan → Green → Yellow | intermediate |
| Red | maximum hit energy in event |

Hits with no pad-mapping entry (sentinel value −999) are filtered out
before display.

---

## Files

```
event_display/
├── ecal_event_display.py       SiW-ECAL TEve single-event display
├── launch.sh                   Convenience launcher
└── conversion/
    └── slab_z_positions.yml    Slab z positions + W absorber thicknesses
```
