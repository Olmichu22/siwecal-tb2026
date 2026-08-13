# trains_out_of_spill

Event-builds the isolated, high-occupancy acquisitions that sit **between** spill
trains, so they can be opened in `event_viewer`.

`diagnostics/plot_beam_spill_structure.py` finds that on `run_000060` (th220) 3 of
the 29 above-cut blocks are single acquisitions outside any train, and
`diagnostics/plot_isolated_bursts.py` shows they are three different things — none
of them a cosmic shower. This pipeline puts those three in front of the viewer.

```bash
./gaudi_jobs/trains_out_of_spill/run_pipeline.sh
python3 -m event_viewer --file gaudi_jobs/trains_out_of_spill/data/ecal_trains_out_of_spill.root
```

A few minutes on one core; no Condor. Override `RUN` and `ACQS` to point it at
other acquisitions or another run.

## Stages

1. **`extract_acquisitions.py`** — scans the run's 59 chunks on `acqNumber` alone
   and copies just the wanted entries into one small `siwecaldecoded` file.
2. **`steer_trains_out_of_spill.py`** — the ordinary `EcalEventBuilder` over that
   file with the full th220 calibration from `calibration/MuonCalib_gaudi/`.

## Two things that are easy to get wrong here

**`MinSlabsHit` must be relaxed, and this pipeline does it deliberately.** The
default is 10 of 15 slabs — right for physics, wrong here: two of the three
acquisitions are confined to 2 and 3 slabs. With the default, building these
three yields 7 events and *all* of them come from the one detector-wide
acquisition; the other two disappear with no warning. This pipeline sets
`MinSlabsHit=1` (`TOS_MIN_SLABS_HIT`). The events it produces are therefore **not**
the physics selection and are not comparable with a normal `ecal_*.root`; many are
single-channel noise by construction.

**Copy through a `TChain`, never `TTree::CloneTree` per file.** `CloneTree(0)`
binds the output tree's branch addresses to the *source* tree's buffers, so on
the next chunk `Fill()` silently writes stale/zero entries. The first version of
stage 1 did that and produced a file whose 2nd and 3rd entries were empty
(`acqNumber=0`, no hits) — which the builder then reported as a chunk-boundary
splice, because two consecutive entries shared `acqNumber` 0. `TChain::LoadTree`
re-copies addresses to registered clones on every file switch.

Related: an acquisition split across a chunk boundary is written twice as two
complementary halves sharing an `acqNumber`, and the builder merges them by
comparing with the previous entry. Both halves must be copied and must stay
consecutive — copying in chain order preserves that.

## Reading the output

`spill` in the `ecal` tree is the builder's **entry index in its input**, not
`acqNumber`. Stage 1 prints the mapping; for the default three:

| spill | acquisition | what it is |
|---|---|---|
| 0 | 2631 | detector-wide activity at many times — beam halo/leakage |
| 1 | 16038 | slab 12 (the COB) firing coincidentally in one instant |
| 2 | 26053 | a chip in slab 1 retriggering across the acquisition |

Two things worth knowing before drawing conclusions from what the viewer shows:

* **Energies are meaningless for these events.** Spill 1's `bcid 2998` comes out
  at 726 channels and ~69,000 in energy; spill 0 has an event at −220. Negative
  totals are what a coherent baseline excursion looks like after per-channel
  pedestal subtraction, which is the same common-mode effect studied in section 7
  of `diagnostics/pedestal_multipeak_tables.txt`. The calibration is not built for
  this regime.
* **Channel counts here do not match the `nhits` numbers in the diagnostics
  plots.** `plot_isolated_bursts.py` reports 236 hits for acq 16038 from the
  chip-reported `nhits` field; the builder finds 726 channels for the same
  acquisition from the hit bits. The two count different things and the gap has
  not been chased down — do not quote them interchangeably.
