# Plan — Particle-discrimination metrics & plots

**Goal:** add shower-shape / topological metrics that help separate particle
types (MIP-like muons vs EM showers vs hadronic showers) in the validation
output, porting the physics from `VarsJesus/analysis_template.C` but adapted to
**the variables we actually have**.

Status: PLAN ONLY — nothing is implemented until approved.

---

## 1. What we have vs what the C++ template assumes

The template reads an `ecal` tree with per-hit vectors
`hit_energy, hit_x, hit_y, hit_z, hit_slab, hit_isMasked` plus a per-layer
tungsten-thickness vector `W_thicknesses`.

Our `ecal` tree (verified on `ecal_TB2026CERN_run_000007.root`) exposes:

- **Per-event scalars:** `run, event, spill, bcid, nhit_slab, nhit_chip,
  nhit_chan, sum_hg, sum_energy`.
- **Per-hit arrays (length `nhit_chan`):** `hit_slab, hit_chip, hit_chan,
  hit_sca, hit_hg, hit_lg, hit_energy, hit_x, hit_y`.

Today `EventData.from_root` ([event_data.py](event_data.py)) reads only
`nhit_chan, hit_slab, hit_energy` and stores per event: `nhit`, `energy`
(Σ MIP), `zbary` (E-weighted layer), `hits_per_layer`, `mip_likeness`.

### Gaps / differences (important)
- **No `hit_z`.** Not needed for transverse metrics; where a real z is wanted,
  derive `z = hit_slab * layer_z_mm` (we already keep `layer_z_mm = 10`).
- **No `hit_isMasked`.** The template's `masked` path is simply unused
  (`masked = False`); no functionality lost for us.
- **Energy units differ.** Our `hit_energy` is in MIP units; the template's
  absolute thresholds (`3.`, `5.` in `is_Shower`/`shower_variables`) were tuned
  to its simulation scale and **must be re-tuned** for us → made configurable.

### Tungsten map (provided — enables Tier B)
From `Tungsten_thickness.yml`. The beam enters and crosses one absorber plate
*in front of* each silicon slab; that local thickness is what weights the layer.
Thicknesses: `t1=2.8`, `t2=4.2`, `t3=5.6` mm. Resolving the `structure:` list to
the absorber directly preceding each slab gives the **per-slab W [mm]**:

| slab | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 |
|------|---|---|---|---|---|---|---|---|---|---|----|----|----|----|----|
| W mm |2.8|4.2|4.2|4.2|4.2|4.2|4.2|4.2|4.2|5.6|5.6|5.6|5.6|5.6|5.6|
| X0   |0.8|1.2|1.2|1.2|1.2|1.2|1.2|1.2|1.2|1.6|1.6|1.6|1.6|1.6|1.6|

The template's weight is `E_hit * W[slab] / 3.5`, where **3.5 mm ≈ X0 of
tungsten**, so the factor is the absorber depth in radiation lengths (bottom
row). `weighte = Σ_i E_i · (W[slab_i]/X0)` is therefore a sampling-corrected
energy estimator (thicker plates sample more shower per layer). We will read the
map from `Tungsten_thickness.yml` so it stays the single source of truth, with
the resolved per-slab array baked into `PlotConfig` as the default.

### Thresholds — what they are used for (re-tuned, configurable)
Applied to the per-layer profile (nhit-per-layer or energy-per-layer):
- **`shower_e_threshold` (template `3.`)** — minimum per-layer activity to count
  as shower material; used both to flag a shower (rising edge) and to locate the
  shower **start/end** layers.
- **`shower_max_min` (template `5.`)** — minimum height of the profile peak for
  the event to qualify as a shower at all.
- **`shower_start_frac` (template `0.1`)** — the `start_10`/`end_10` variants
  also require the layer to exceed 10% of the peak.

**Profile choice & noise.** We anchor the shower flag / start-end to the
**nhit-per-layer** profile (geometric, ~scale-independent). But a MIP layer is
**not** a single hit here: detector **noise adds extra hits per layer**, so the
per-layer floor for a MIP sits a few hits above 1. The threshold must therefore
clear that noise+MIP floor. Recommended defaults (to be confirmed against the
muon run, which *is* the noise+MIP floor):
- `shower_e_threshold = 5` hits/layer (was template `3`; raised above noise).
- `shower_max_min = 10` (peak must sit clearly above the per-layer threshold;
  raised in step with it so the "rising profile with a real peak" logic stays
  meaningful).
- `shower_start_frac = 0.1`.

Energy-based profiles (`sume`, `weighte`) keep MIP-scale thresholds; for
`weighte` the threshold is scaled by `W[slab]/X0` so it means the same per layer.

---

## 2. Variable mapping (C++ → ours) and viability

| Template quantity | Inputs needed | Viable for us? |
|---|---|---|
| `nhit`, `sume` | hit_energy | ✅ already `nhit`, `energy` |
| `MIP_Likeness` | nhit_layer | ✅ already `mip_likeness` (same formula) |
| `nhit_layer[]` | hit_slab | ✅ already `hits_per_layer` |
| `sume_layer[]` (E profile) | hit_slab, hit_energy | ➕ add (one pass) |
| barycenter `bar_x/y/r` (E-weighted) | hit_x, hit_y, hit_energy | ✅ add |
| `bar_z` | slab×layer_z_mm | ✅ have `zbary` (layer units) |
| `is_Shower` + shower start/max/end/length | sume_layer / nhit_layer | ✅ add |
| **Molière radius** (90% transverse) | hit_x, hit_y, hit_energy | ✅ add (key) |
| transverse E-weighted RMS radius | hit_x, hit_y, hit_energy | ✅ add (cheap) |
| first/last hit layer, `n_layers_hit` | hit_slab | ✅ add (cheap) |
| `weighte`, `weighte_layer[]`, W-weighted shower vars | + W map (provided) | ✅ Tier B |
| `radius90_layer[]`, `bary_layer[]` (15-dim) | hit_x, hit_y, hit_energy | ⚠️ Tier C (heavy) |

---

## 3. Proposed metrics (tiered)

**Tier A — high discrimination, only uses what we already read.** Recommended
first deliverable.
- `energy_per_layer[n_layers]` (longitudinal E profile) — enables the rest.
- `moliere_radius` — transverse 90% containment radius. Strong e/μ/hadron
  separator (μ≈point, e compact, hadron wide).
- `transverse_rms` — E-weighted RMS of hit radius about the shower axis (cheap
  proxy for compactness, no sorting).
- `bar_x`, `bar_y`, `bar_r` — E-weighted transverse barycenter (also QA).
- `is_shower` (bool) + `shower_start_layer`, `shower_max_layer`,
  `shower_end_layer`, `shower_length` (= layers with E above threshold).
- `first_layer`, `last_layer`, `n_layers_hit` — longitudinal extent.
- `e_over_nhit` = energy / nhit — cheap density discriminator.

**Tier B — IN SCOPE (tungsten map provided).**
- `weighte = Σ E_i·W[slab_i]/X0` (sampling-corrected total energy) and
  `weighte_per_layer[]`.
- W-weighted shower variables (radiation-length-correct longitudinal
  development), and the `weighte` profile feeding `is_shower`/shower start-max-
  end. These are the most physically meaningful longitudinal discriminators
  because they account for the non-uniform absorber.

**Tier C — optional, heavier (15-dim per event).**
- `radius90_per_layer[]`, `bary_per_layer[]` — per-layer transverse profiles.

---

## 4. Design / where each change goes

Keep the existing architecture: compute per-event scalars (and fixed-width
per-layer arrays) **once** in `EventData.from_root`, then each Plotter only
reads `EventData` fields. This mirrors both the current code and the C++ (which
computes per event, then fills).

1. **`event_data.py` — single-pass computation.**
   - Read `hit_x`, `hit_y` alongside `hit_slab`, `hit_energy` per event.
   - Add a small pure-numpy helpers module (`metrics.py`) with vectorised
     equivalents of: `energy_per_layer`, `barycenter_xy`, `moliere_radius`,
     `transverse_rms`, `shower_variables` (on a 1-D layer profile). Keeping them
     in their own file makes them unit-testable in isolation.
   - Extend the `EventData` dataclass with the new fields and update `select()`
     to slice them (all aligned per event).

2. **`metrics.py` (new) — single home for ALL metric logic.** One function per
   metric, each taking plain arrays (`slab, energy, x, y`) for a single event (or
   a layer profile). Each carries the **formula in the docstring** (e.g. Molière
   90% containment, RMS radius). Edge cases mirrored from the C++: empty event →
   0; `sume<=0` → barycenter 0; `n_layers_hit<3` ignored for per-layer radius;
   Molière only when `is_shower`.
   - **The existing/old metrics are moved here too (if not already):** the
     current per-event computations living inline in `EventData.from_root` —
     `nhit` (Σ channels), `energy`/`sume` (Σ hit_energy), `zbary` (E-weighted
     layer), `hits_per_layer` (bincount), and `mip_likeness` (mean over hit
     layers of `1/hits_in_layer`) — are extracted into `metrics.py` as named,
     documented, unit-tested functions, so `from_root` becomes a thin loop that
     just calls `metrics.*`. This keeps one consistent, testable home for both
     old and new variables and avoids duplicated formulas. Behaviour is
     preserved exactly (regression-checked against current outputs).

3. **`config.py` `PlotConfig` — new tunables (all defaulted, documented):**
   - `w_thicknesses: tuple` — per-slab W [mm], default = the resolved map above
     `(2.8, 4.2×8, 5.6×6)`; `w_x0_mm = 3.5` (tungsten X0). A small loader parses
     `Tungsten_thickness.yml` (`W_thick` + `structure`) into the per-slab array
     so the map stays authoritative; the baked tuple is the fallback default.
   - `shower_profile = "nhit"` (base profile for the shower flag / start-end),
     `shower_e_threshold = 5` (hits/layer, raised above the noise+MIP floor —
     noise gives a MIP layer several hits, so `3` is too low),
     `shower_max_min = 10`, `shower_start_frac = 0.1` — all MIP-/hit-scale and
     configurable. Energy/`weighte` profile thresholds expressed in MIP and, for
     `weighte`, scaled by `W[slab]/X0`.
   - `moliere_containment = 0.9`, `moliere_bins`, `radius_bins`, profile bins.

4. **`plots.py` — new Plotters** (each auto-appears as PNG + grid + per-energy
   grid; no other wiring needed):
   - `LongitudinalProfilePlotter` — mean E per layer (line/step).
   - `MoliereRadiusPlotter` — 1-D hist.
   - `TransverseRmsPlotter` — 1-D hist.
   - `ShowerStartLayerPlotter`, `ShowerMaxLayerPlotter` — 1-D hist.
   - `NLayersHitPlotter` — 1-D hist.
   - `Discriminator2DPlotter` — 2-D hist of two separating vars (default
     **Molière radius vs shower-start layer**), the single most useful view for
     particle-type separation.
   - Register in `DEFAULT_PLOTTERS`.

5. **`selection.py` `CutSet` (IN SCOPE — every new variable is selectable).**
   Extend `_CUT_SPEC` with min/max bounds for **all** new per-event scalars so
   they can be cut both from the CLI and from the per-energy YAML `cuts:` blocks
   already supported: `weighte`, `moliere`, `transverse_rms`, `bar_r`,
   `shower_start` (= `shower_start_layer`), `shower_max` (layer), `shower_end`,
   `shower_length`, `n_layers_hit`, `first_layer`, `last_layer`, `e_over_nhit`,
   plus `is_shower` as a boolean filter. Each needs: a field pair in `CutSet`, a
   `_CUT_SPEC` row (short name + `EventData` attr), and a short token for the
   file-name suffix/title. (Per-layer arrays in Tier C are not cut variables.)

6. **`results.py` (IN SCOPE).** Add summary columns: mean Molière radius, mean
   `weighte`, shower fraction (`is_shower` rate), mean shower-start layer.

7. **`cli.py` (IN SCOPE — follows from cuts).** Add `--<var>-min/--<var>-max`
   flags mirroring every new `CutSet` field, consistent with the existing
   `--nhit-min` etc., wired through `cutset_from_args`.

No change to `runner.py` is required unless we add new cross-energy summary
plots; if we do, it follows the existing `mu_fits`/`sigma_fits` collection
pattern.

---

## 5. Phased execution

Scope for this implementation: **Tier A + Tier B**.

- **Phase 0 — scaffolding:** `metrics.py` (Tier A + B functions, incl. `weighte`
  and W-weighted profiles) + unit tests; extend `PlotConfig` with the W map
  loader and thresholds. No behaviour change yet.
- **Phase 1 — data:** extend `EventData.from_root`/`select` to read `hit_x/y`,
  apply the per-slab W weights, and compute & carry all Tier A+B fields.
  Validate numbers on a known event.
- **Phase 2 — plots:** add the Plotters (including the `weighte` longitudinal
  profile and W-weighted shower vars) + register them.
- **Phase 3 — cuts/results integration:** expose **every** new per-event
  variable as a `CutSet` field + CLI flag + YAML `cuts:` key, and add the
  summary columns. Validate a cut end-to-end (suffix, title, mask).
- **Phase 4 — Tier C** (only if per-layer profiles are later wanted).

Each phase is independently runnable and reviewable.

---

## 6. Validation strategy (rigour)

- **Unit tests** (`metrics.py`): synthetic events with known answers — single
  hit → Molière 0 & RMS 0; symmetric ring → known radius; flat vs peaked layer
  profile → expected shower start/max; barycenter of a known geometry.
- **Cross-check** one real event end-to-end against a hand/`numpy` computation.
- **Physics sanity (the real proof):** run on the **muon** run (`run_000004`)
  vs an **electron** energy point and confirm the expected separation — muons:
  low Molière, `is_shower` mostly False, MIP-like; electrons: showering, larger
  Molière, shower start in early layers. This is the acceptance criterion.
- **Regression:** existing plots/metrics must be unchanged (the new code only
  adds fields/plotters).

---

## 7. Resolved decisions

All scope questions are settled; defaults below, tunable at review.

- **Tungsten map:** provided (`Tungsten_thickness.yml`) → **Tier A + B in scope**.
- **Threshold profile & defaults:** anchor the shower flag / start-end to the
  **nhit-per-layer** profile, with `shower_e_threshold = 5` (raised above the
  noise+MIP floor — with noise a MIP layer has several hits, so `3` is too low),
  `shower_max_min = 10`, `shower_start_frac = 0.1`. Energy/`weighte` thresholds
  in MIP (the latter scaled by `W[slab]/X0`). All in `PlotConfig`; confirmed
  against the muon run (the noise+MIP floor) before freezing.
- **Cuts/results:** **all** new per-event variables are exposed as selectable
  cuts (CLI + YAML) and the summary columns are added.
- **2-D discriminator axes:** **Molière radius vs shower-start layer.**
