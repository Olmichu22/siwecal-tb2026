# SiW-ECAL Event Viewer (TB2026CERN)

Interactive Plotly Dash app to inspect reconstructed `ecal` events: a 3-D
detector scene (silicon planes + tungsten absorbers) with the per-pad signals
coloured by energy, a 2-D per-layer view, a lateral table of per-event metrics,
plus a distributions tab with dynamic cuts and simple clustering.

## Load the environment

```bash
source setup.sh                      # from the repo root: key4hep + .venv-viewer
```

The `.venv-viewer` is a `--system-site-packages` venv on top of key4hep (which
already provides uproot, pandas, scikit-learn, ...) with only `dash` + `plotly`
added. It is **not** versioned (147 MB, with absolute CVMFS paths baked in);
create it once per machine:

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2026-04-08
python -m venv --system-site-packages .venv-viewer
source .venv-viewer/bin/activate
pip install -r requirements.txt      # dash + plotly
```

After that, `source setup.sh` activates it automatically.

## Run

From the **repo root** (so `siwecal_eventbuilder` / `siwecal_validation` /
`siwecal_common` are importable):

```bash
python -m event_viewer                                   # pick a file from the dropdown
python -m event_viewer --file /path/to/ecal_run.valcache.root   # open one directly
python -m event_viewer --data-dir /some/other/dir        # scan a custom directory
```

Then open the UI from your laptop through an SSH tunnel:

```bash
ssh -L 8050:localhost:8050 you@lxplus.cern.ch
# browse to http://localhost:8050
```

Without `--file` the viewer scans every `data_dir` configured in `settings.yml`
(plus `cache_dir` if set) for `*.valcache.root` and `ecal_*.root`, and you pick a
run from the dropdown or change file at runtime.

## Tabs

- **Event** — 3-D scene + 2-D layer grid + metrics table. Navigate with the
  *Previous*/*Next* buttons or by typing an index. Only events passing the active
  cuts are reachable.
- **Distributions** — file-level histograms; add **dynamic cuts** by picking
  variables (a `RangeSlider` appears per variable) which also filter the Event
  tab; run **clustering** (K-Means / DBSCAN / GMM / Spectral) on chosen features
  and view the labels on any 2-variable scatter.

## Files without metrics

Opening a plain `ecal_*.root` (no `.valcache.root`) still works: the viewer
computes only the cheap derivable quantities (`nhit`, `sum_energy`,
`n_layers_hit`, `first/last_layer`, `zbary`, `mip_likeness`, `e_over_nhit`) and
the status bar flags the file as "SIN métricas". The event and layer views are
unaffected.

## Architecture (OOP)

```
config.py        ViewerConfig (paths, display constants)
io/              EventFileReader            — uproot, per-event vs per-hit split
model/           DetectorModel, Event, EventDataset
analysis/        CutModel, ClusteringService, compute_cheap_metrics (reuses metrics.py)
viz/             DetectorScene3D, LayerGrid2D, DistributionPlots — Plotly figures
controller.py    ViewerController           — geometry + builders + per-file cache
ui/              build_layout, register_callbacks — Dash
app.py           build_app factory; __main__.py CLI
```

Reuses the existing geometry (`siwecal_eventbuilder.geometry` / `pad_map`) and
physics (`siwecal_validation.metrics`); reads with uproot (no PyROOT needed).

## Troubleshooting

- **Distributions / clustering "stuck" on a previous selection.** A clustering
  run is tied to the cuts it was computed under; changing the cuts invalidates it
  (the histogram, scatter and cluster examples then follow the current selection
  again). Just re-run clustering after changing cuts.
- **The UI occasionally gets into an inconsistent state** after a particular
  clustering ↔ cuts ↔ event-navigation sequence (a Dash client/server state race
  we could not reproduce deterministically). If a panel looks out of sync,
  **refresh the browser tab** (or restart `python -m event_viewer`) — that always
  resolves it. State is rebuilt from the files, so nothing is lost.