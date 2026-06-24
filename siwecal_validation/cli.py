"""
Command-line interface for the SiW-ECAL validation plots.

Preserves the original ``plot_val_plots.py`` flags:

    --run RUN            single run (resolves to <base>/<run>/ecal_<run>.root)
    --file FILE          explicit path to an events ROOT file
    --all                process every entry of the YAML 'event_data' map
    --cfg FILE           data_reference YAML
    --out DIR            output base (external; default: validation_output)
    --grid-only          produce only the combined grid PNG (no individuals)
    --no-grid            produce only the individual PNGs (no combined grid)
    --nhit-min/max, --energy-min/max, --rate-min/max   general selection cuts

Cache / augmented-tree flags:

    --create-tree        force (re)generate the per-file metrics cache
                         (``<stem>.valcache.root`` next to the input); useful
                         after a config change or to pre-warm a fresh dataset.
    --save-tree          also write a cut-applied augmented TTree in the output
                         directory (``<out>/<label>/trees/<label><suffix>.root``).

General cuts apply to every sample; in ``--all`` mode a per-energy ``cuts:``
block in the YAML overrides them, so selections can be optimised per energy.
"""

import argparse
import os
from dataclasses import fields

import yaml

from siwecal_common import paths

from .config import (BASE_PATH, DEFAULT_CONFIG, DEFAULT_OUTPUT_DIR,
                     DEFAULT_RUN, PlotConfig)
from .output import OutputLayout
from .runner import ValidationRunner
from .selection import CutSet, _CUT_SPEC


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="SiW-ECAL validation plots (object-oriented).")
    p.add_argument("--run", default=None,
                   help="Single run name. With --file it only sets the output "
                        "label (default: derived from the file name).")
    p.add_argument("--file", default=None,
                   help="Explicit path to events ROOT file (ecal tree)")
    p.add_argument("--cfg", default=DEFAULT_CONFIG,
                   help="data_reference YAML (run/energy -> events mapping)")
    bulk = p.add_mutually_exclusive_group()
    bulk.add_argument("--all", action="store_true",
                      help="Process all entries of the YAML 'event_data' map")
    bulk.add_argument("--point", type=int, choices=range(1, 6), default=None,
                      help="Process only one beam point (P<N>_... entries), like "
                           "--all restricted to that point so the summary plots "
                           "compare energies within it.")
    p.add_argument("--out", default=DEFAULT_OUTPUT_DIR,
                   help="External output base directory")
    # Combined grid vs individual plots (default: produce both).
    grid = p.add_mutually_exclusive_group()
    grid.add_argument("--grid-only", action="store_true",
                      help="Only the combined grid PNG (skip individual plots)")
    grid.add_argument("--no-grid", action="store_true",
                      help="Only the individual PNGs (skip the combined grid)")
    # General selection cuts (applied to all samples; per-energy YAML overrides).
    # One --<var>-min/--<var>-max float flag per CutSet variable, generated from
    # the spec so the CLI never drifts out of sync with selection.py.
    cuts = p.add_argument_group(
        "selection cuts",
        "min/max bounds applied to all samples (per-energy YAML cuts override)")
    for field_min, field_max, _name, _attr in _CUT_SPEC:
        cuts.add_argument(f"--{field_min.replace('_', '-')}",
                          type=float, default=None)
        cuts.add_argument(f"--{field_max.replace('_', '-')}",
                          type=float, default=None)
    shower = cuts.add_mutually_exclusive_group()
    shower.add_argument("--is-shower", dest="is_shower", action="store_const",
                        const=True, default=None, help="keep only shower events")
    shower.add_argument("--no-shower", dest="is_shower", action="store_const",
                        const=False, help="keep only non-shower events")
    # Augmented-tree / cache flags.
    p.add_argument("--create-tree", action="store_true", dest="create_tree",
                   help="Force (re)generate the per-file metrics cache next to "
                        "the input (*.valcache.root). Use after a config change "
                        "or to pre-warm a fresh dataset.")
    p.add_argument("--save-tree", action="store_true", dest="save_tree",
                   help="Write a cut-applied augmented TTree to the output "
                        "directory (<out>/<label>/trees/<label><suffix>.root).")
    p.add_argument("--cache-dir", default=paths.cache_dir(), dest="cache_dir",
                   help="Directory for the *.valcache.root metric caches "
                        "(default: settings.yml cache_dir, else next to each "
                        "input). Use this when the data directory is read-only.")
    return p.parse_args(argv)


def label_from_events_path(events_path: str) -> str:
    """Derive an output label from an events file name.

    ``.../ecal_TB2026CERN_run_000048.root`` -> ``TB2026CERN_run_000048``. The
    ``ecal_`` prefix and ``.root`` suffix are stripped; otherwise the bare stem
    is used.
    """
    stem = os.path.splitext(os.path.basename(events_path))[0]
    if stem.startswith("ecal_"):
        stem = stem[len("ecal_"):]
    return stem


def cutset_from_args(args) -> CutSet:
    """Build the general CutSet from the CLI flags (one per CutSet field)."""
    return CutSet(**{f.name: getattr(args, f.name)
                     for f in fields(CutSet) if hasattr(args, f.name)})


def main(argv=None):
    args = parse_args(argv)
    general_cut = cutset_from_args(args)
    layout = OutputLayout(args.out)
    make_individual = not args.grid_only
    make_grid = not args.no_grid
    runner = ValidationRunner(layout, config=PlotConfig(),
                              make_individual=make_individual,
                              make_grid=make_grid,
                              create_tree=args.create_tree,
                              save_tree=args.save_tree,
                              cache_dir=args.cache_dir)

    if args.all or args.point is not None:
        with open(args.cfg) as handle:
            cfg = yaml.safe_load(handle) or {}
        event_data_map = cfg.get("event_data", {})
        if not event_data_map:
            raise SystemExit(f"ERROR: no 'event_data' mapping in {args.cfg}")
        if args.point is not None:
            prefix = f"P{args.point}_"
            event_data_map = {label: entry
                              for label, entry in event_data_map.items()
                              if label.startswith(prefix)}
            if not event_data_map:
                raise SystemExit(
                    f"ERROR: no entries for point P{args.point} in {args.cfg}")
            print(f"[Point] P{args.point}: {len(event_data_map)} energy point(s): "
                  f"{', '.join(event_data_map)}")
        base_path = cfg.get("main_path", BASE_PATH)
        runner.run_all(event_data_map, base_path, general_cut=general_cut)
    else:
        if args.file:
            # Read the given file; label from --run, else derived from the name.
            events_path = args.file
            label = args.run or label_from_events_path(events_path)
        else:
            run = args.run or DEFAULT_RUN
            # Search every data root (settings.yml data_dir) for the run file.
            events_path = paths.resolve_input(
                os.path.join(run, f"ecal_{run}.root"))
            label = run
        runner.run_sample(events_path, label=label, cutset=general_cut)
        runner.write_results()

    print(f"\n[Done] output -> {args.out}")


if __name__ == "__main__":
    main()
