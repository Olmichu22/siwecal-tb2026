# Structured execution plan — from claude_orders.txt

Rewritten, dependency-ordered orders to myself. Threshold points: **th230, th220, th210**.
Core correction driving everything: **MIPs and pedestals shift with threshold, so every
threshold needs its OWN pedestal + MIP fit.** Never reuse another threshold's ADC MIP value.

## Phase A — Data readiness (decode th210 muons via LCIO)
- **A1. Finish the LCIODecoder path.** The `EcalLcioDecoder` component + `run_lcio_decode.py`
  steering are committed and validated (exact-match vs raw on run_000254). What is missing is a
  *batch launcher* that feeds the th210 muon runs through it, since the raw decoder cannot parse
  their EUDAQ `.bin` (runs 146-178 are LCIO-only in `Data/eudaq/ROC_run_<N>_tp.slcio`).
  → Add a small batch/DAG generator alongside the existing decode jobs.
- **A2. Decode th210 muon runs.** Confirm 152-155 and 158-160 are muon runs, then decode each
  from its `ROC_run_000NNN_tp.slcio` into a `siwecaldecoded` ROOT tree.

## Phase B — Per-threshold pedestal + MIP calibration
GOVERNING DIRECTIVE (user, 2026-07-15): **every event-builder pipeline auto-calibrates from its
OWN threshold — both pedestal AND MIP. The th220-for-all-thresholds MIP bypass is retired.**
- **B0. [DONE] Retire the bypass in the resolver.** `MIP_CALIB_TH` is now `None` = "run's own
  threshold"; `resolve_gaudi_calib_files` defaults MIP to `th`. `--mip-th` stays only as an
  explicit cross-threshold override. Touched `siwecal_eventbuilder/cli.py`,
  `gaudi_jobs/run_full_pipeline_batch.py`, `gaudi_jobs/condor/generate_reco_dag.py`.
- **B1. [CODE DONE] Falling-edge MIP fit — OPTIONAL, th230 only.** th230's old MIP table
  (run_000004, MPV 32.40) is the trigger-cut-biased fit. Added an OPT-IN high-gain falling-edge
  fit (`fitLangau(..., fallingEdgeHighGain)`, property `MipFallingEdgeHighGain`, env
  `CALIB_MIP_FALLING_EDGE`) that fits only the peak-upward Landau tail. Default OFF; both the
  calibration DAG (`generate_dag.py`) and interactive batch (`run_calibration_batch.py`) enable it
  ONLY when `th == 230`. th220/th210 keep the full-range fit (their cut sits below the peak and
  already matches the reference tool). REMAINING: regenerate th230's deployed MIP table with the
  flag on (validation of MPV shift in progress).
- **B2.** Produce the th210 pedestal + MIP tables (HG and LG) from the A2 decode, using the same
  MuonCalib batch. This is the missing third threshold.

## Phase C — Diagnostic plots
- **C1.** Pedestals + MIP full-mapping plots (complete channel map) for th230, th220, th210.
- **C2.** HG & LG fit gallery — one full chip per case (many individual fits shown).
- **C3.** MIP stability plot for th230 — NOT POSSIBLE: only run_004 is a th230 muon run
  (confirmed by user 2026-07-15), so there is no run-to-run time series to plot. th220 stability
  already exists and stands.
- **C4.** Redo `plot_mip_snr` dropping the redundant 4th panel (noise-scaling); keep panels 1-3.
- **C5.** LG-vs-HG fit plot, one per threshold (defends the per-threshold-fit correction).

## Phase D — Energy calibration
- **D1.** Per threshold: use one electron run to FIT the calibration, a second to TEST it
  (extrapolation + HG→LG switch). **Switch at 1500 ADC of high-gain.**
- Electron runs (fit / test):
  - th230: 12 / 13
  - th220: 91 / 92
  - th210: 166 / 167

## Execution order & status
1. C4 (mip_snr 3-panel, th220/230/210) .... DONE (+ user text tweaks: legend size, "pedestal vs MIP", no "noise", 2-dec %)
2. B0 (retire th220 MIP bypass) ........... DONE (resolver verified per-th)
3. A1 LCIO batch launcher ................. DONE (decode_lcio_runs.py)
4. A2 decode th210 muons (152-5,158-60) ... DONE (7/7 valid, 511,794 acqs)
5. B1 falling-edge fit ................... ABANDONED/REVERTED per user; th230 uses full-range like all
6. B2 th210 muon calib .................... DONE (tables written; HG MPV 18.56, 211 masked)
7. C1 full ped+MIP maps ................... DONE (th220, th230, th210) plot_calib_maps.py
8. C2 chip-shapes HG+LG ................... DONE (th220, th230, th210; s13/c5 & s2/c5) plot_mip_chip_shapes.py
9. C3 th230 stability ..................... NOT POSSIBLE (only run_004 is th230 muon)
10. Switch on ped-subtracted ADC .......... DONE (EventBuilder.h, compiled, 21 tests pass)
11. gain_anchor edits ..................... DONE (1 panel, no least-squares, "fiducial region")
12. Trend found: MIP HG th210=18.56 < th220=21.69 < th230=32.40 (truncation with threshold)

## Remaining
- D energy calibration (switch HG->LG at 1500 ADC ped-sub hg). Electron runs fit/test:
  th230 12/13 (decoded), th220 91/92 (raw, 630+254 chunks -- need decode),
  th210 166/167 (LCIO -- need decode). Then reconstruct each & fit energy scale + validate.
- Save muon-run-by-threshold list to a .txt (waiting on user).
- C5 "LG vs HG fit per threshold": effectively covered by the C2 chip-shapes (HG & LG per th).
