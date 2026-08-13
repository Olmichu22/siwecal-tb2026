# ACTS Track Reconstruction Integration

## Overview

This repo reconstructs tracks through the **15 silicon slabs of the SiW-ECAL
prototype**, using them as tracking surfaces — real test-beam hits, not
simulated ones. There is one detector, "SiPad" (the same name
`../siwecal_k4sim` uses for the equivalent silicon pad detector in the
simulation).

This is a direct port of `../siwecal_k4sim`'s tracking stage: **same class
names, same property names, same address convention, same tuned Hough/CKF
defaults** — see that repo's `docs/acts_integration.md` and
`docs/gaudi_pipeline.md`, which this document mirrors section-for-section. The
intent is that a run's reconstructed tracks and a simulated sample's are
directly comparable without translating variable names or conventions. The
**one real difference** is where the detector geometry comes from (see
"ACTSGeoSvc" below) — everything downstream of the geometry service is
unmodified between the two repos, down to internal C++ type names.

ACTS version: 44.3.0 (from the key4hep stack pinned in `.key4hep-release`,
shared with `siwecal_k4sim`).

Four components, all in `gaudi_source/src/components/`:

| Component | Kind | Role |
|---|---|---|
| `ACTSGeoSvc` | Gaudi service | Builds the `Acts::TrackingGeometry`: one plane surface per slab, with that slab's tungsten-absorber material |
| `ShowerTagger` | Gaudi algorithm | Identifies EM cascades and flags their hits so they never become measurements; writes the showers |
| `SiPadMeasConverter` | Gaudi algorithm | Pad hits → `edm4hep::TrackerHit3D` measurements bound to those surfaces, minus the flagged ones |
| `ACTSProtoTracker` | Gaudi algorithm | Hough seeding → CKF → KalmanFitter refit → `edm4hep::TrackCollection` |

Configured by `gaudi_source/options/run_tracking.py`; see
`gaudi_source/README.md`'s "Tracking (ACTS)" section for the job's
environment variables and how it fits the rest of the pipeline.

> **Heritage note.** `ACTSProtoTracker` is copied unmodified from
> `siwecal_k4sim` (a fork of `key4ship_PoC`, which tracked three subdetectors).
> Anything mentioning stations, planes or U/V pairing is vestigial: here
> `detID` is always 1, `station`/`plane` are always −1, and every measurement
> is a 2D pad hit. Internal type names (`SND*`) are kept as-is on purpose —
> matching names make the two repos' sources diffable against each other.

---

## ACTSGeoSvc

**Files:** `gaudi_source/include/k4SiWEcalReco/ACTSGeoSvc.h`,
`gaudi_source/src/components/ACTSGeoSvc.cpp`, interface
`gaudi_source/include/k4SiWEcalReco/ISNDGeoSvc.h`.

### Where the geometry comes from — the one real difference vs. siwecal_k4sim

`siwecal_k4sim`'s `ACTSGeoSvc` walks a DD4hep compact-XML detector
description to find its 15 sensitive slices. This repo has **no DD4hep
detector description** — the real geometry is `mappings/slab_z_positions.yml`
(per-slab z position and tungsten thickness), already ported to C++ as
`k4siwecal::SlabGeometry` in `gaudi_source/include/k4SiWEcalReco/PadMapGeometry.h`
(the same header `EcalEventBuilder` uses for `hit_z`/`hit_X0`). So this
`ACTSGeoSvc` loops `hit_slab` 0..14 directly over `SlabGeometry` instead of a
`DetElement` tree, and every downstream piece (the `PlaneLayer`/
`CuboidVolumeBounds`/`LayerArrayCreator` scaffolding, the `rot90Y` convention,
the address map) is copied unchanged.

### Configuration

| Property | Description |
|---|---|
| `SlabZFile` | Path to `mappings/slab_z_positions.yml` (absolute — the job may run from any directory). Replaces `siwecal_k4sim`'s `CompactFile` |
| `PadPitchMm` | Pad pitch [mm], default `5.53` — same physical hardware as the simulation |
| `HalfSizeX` / `HalfSizeY` | Transverse half-size [mm] of each plane surface, default `90.0` |

### Interface methods (identical to `siwecal_k4sim`)

| Method | Returns |
|---|---|
| `trackingGeometry()` | The `Acts::TrackingGeometry` |
| `allSurfaces()` | The 15 sensitive surfaces, sorted by ACTS X (the beam axis) |
| `surfaceByAddress(detID, station, layer, plane)` | Surface lookup by address; SiPad uses `(1, -1, hit_slab, -1)` |
| `geometryContext()` | The `Acts::GeometryContext` |

### Surface material

Each surface carries **only its slab's tungsten absorber** (`SlabGeometry::
slabWOverX0(slab) * kWX0Mm` recovers the raw thickness in mm), folded into an
`Acts::HomogeneousSurfaceMaterial` via `Acts::Material::fromMassDensity` with
the same `X0 = 3.5 mm` constant `EcalEventBuilder` uses for `hit_X0`. This is
simpler than `siwecal_k4sim`'s full DD4hep slice-stack query (which also sums
Si/PCB/Cu/carbon-fibre), because W dominates the material budget and this repo
has no per-slice geometry to sum from — the physical quantity (radiation
lengths of tungsten crossed) is the same one either way.

### Unit notes

`mappings/slab_z_positions.yml` is already in mm (unlike DD4hep/TGeo, which
work in cm) — no unit conversion needed on the way in.

---

## Coordinates — the one thing to get right

Identical to `siwecal_k4sim` (see that repo's `docs/acts_integration.md` for
the full derivation). Two frames are in play:

* **Detector convention** (this repo's own, matching `mappings/
  slab_z_positions.yml` and `EcalEventBuilder`'s `hit_x/hit_y/hit_z`): z is the
  beam, (x, y) transverse.
* **ACTS global:** X is the beam (after `rot90Y`), so the transverse plane is
  (Y, Z).

`rot90Y = R_Y(pi/2)` sends the local axes to

```
local x -> (0, 0, -1)     local y -> (0, 1, 0)     local z -> (1, 0, 0)
```

so local z is the beam as intended, but **local x maps to minus global Z**.
`SiPadMeasConverter`/`SNDCalibrator` define the measurement as `(loc0, loc1) =
(hit_x, hit_y)`, so consistency requires `global Y = +hit_y`, `global Z =
-hit_x`, for both the seed position and direction (see `ACTSProtoTracker`'s
seed-construction comments — copied unchanged from the simulation, where the
consequence of getting this sign wrong is spelled out in detail).

On output the convention is inverted back, so `TrackState::referencePoint`
carries this repo's own `(x, y, z_beam)`.

---

## ShowerTagger — why an EM event yields no track

**File:** `gaudi_source/src/components/ShowerTagger.cpp`

Same logic as `siwecal_k4sim`'s `ShowerTagger` (see that repo's doc for the
full rationale — every slab here is 1.2-2.0 X0, so a showering particle
leaves no clean incoming segment). Two adaptations for real data:

* **Input**: `edm4hep::CalorimeterHitCollection` (this repo's `ECalHits`,
  real reconstructed hits from `EcalToEDM4hep`) instead of
  `edm4hep::SimCalorimeterHitCollection`.
* **`BitField`**: defaults to `"system:8,slab:8,chip:16,channel:8,sca:8"` —
  the exact cellID convention `EcalToEDM4hep.cpp` writes — decoding a `"slab"`
  field (the simulation's DD4hep bitfield calls the equivalent field
  `"layer"`).

Properties, outputs and the onset/veto algorithm are otherwise identical:
`ShowerNHitsThreshold` (4), `ShowerMinConsecutive` (2), `MinTrackLayers` (4),
`Enabled`, `OutputFlags` (`SiPadShowerFlags`), `OutputShowers` (`EMShowers`,
`shapeParameters = [start slab, slab of maximum, slabs spanned, transverse RMS
in mm, number of hits]`).

---

## SiPadMeasConverter

**File:** `gaudi_source/src/components/SiPadMeasConverter.cpp`

Same adaptation as `ShowerTagger`: input `edm4hep::CalorimeterHitCollection`
(`ECalHits`) instead of `edm4hep::SimCalorimeterHitCollection`, and hit time
comes straight from `hit.getTime()` (a real field on `CalorimeterHit`) instead
of the simulation's MC-truth-only `hit.getContributions(0).getTime()`.

Everything else is unchanged: `layer`/`slab` is decoded from the cellID with a
`dd4hep::DDSegmentation::BitFieldCoder` and written into `quality` (how
`ACTSProtoTracker` finds the surface, via `surfaceByAddress(1, -1, slab, -1)`);
position is taken from the hit, not recomputed; covariance is `pitch²/12` on
each transverse axis (`PixelSizeX`/`PixelSizeY`, both `5.53` mm); `InputFlags`
names the per-hit veto collection (`ShowerTagger`'s `OutputFlags`).

---

## ACTSProtoTracker

**File:** `gaudi_source/src/components/ACTSProtoTracker.cpp` — an **unmodified
copy** of `siwecal_k4sim/gaudi_source/ACTSProtoTracker.cpp` (only the
`ISNDGeoSvc.h` include path changes). It only ever talks to `ISNDGeoSvc` and
`SiPadMeasConverter`'s `edm4hep::TrackerHit3DCollection` output, neither of
which depends on DD4hep, so nothing about real vs. simulated data reaches this
file. See `siwecal_k4sim/docs/acts_integration.md`'s "ACTSProtoTracker"
section for the full internal-types table (`SNDMeasurement`, `SNDSourceLink`,
`SNDFixedNavigator`, `SNDCalibrator`, ...), the two non-obvious obligations
(`endOfWorldReached()` must report `navigationBreak`; `SNDCalibrator` must call
`setUncalibratedSourceLink()`), and the full algorithm flow (Hough seeding →
shower-hit purge → CKF → final KF refit → chi2/ndf acceptance → seed cleaning
→ event-level deduplication) — all identical here.

`IronFieldRanges` stays empty by default, i.e. zero magnetic field, correct
for a test beam (same as the simulation's beam pipelines).

### Neighbour windows are in pad pitches

`IsolationWindow` and `HitPurgeWindow` must be expressed in units of the pad
pitch (5.53 mm) — any window at or below one pitch counts zero neighbours by
construction, silently disabling the filter. `run_tracking.py` uses 1.5
pitches, same as the simulation.

---

## Interpreting chi2/ndf on this detector

Same geometry, same conclusion as `siwecal_k4sim`: with a 5.53 mm pitch and
~11-30 mm between slabs, a track must be tilted by a sizeable angle to change
pad, so a beam muon deposits in the *same* pad on most slabs and the residual
is small by construction — chi2/ndf sits near zero for straight tracks, and
`MaxChi2PerNdf` is a weak discriminant; the number of measurements, holes and
per-hit residual are the more useful quality handles.

---

## See also

* `gaudi_source/README.md` — "Tracking (ACTS)" section: the job's environment
  variables and how it fits the rest of the pipeline
* `../siwecal_k4sim/docs/acts_integration.md` — the simulation's version of
  this document (DD4hep-based `ACTSGeoSvc`, otherwise the same components)
* `../siwecal_k4sim/docs/gaudi_pipeline.md` — the simulation's job
  configuration and measured tracking performance (muon efficiency, shower
  rejection rates) as a reference point for what "the same physics story" on
  real data should look like
