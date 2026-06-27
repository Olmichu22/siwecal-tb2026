"""
Command-line interface for the SiW-ECAL event builder.

Mirrors the original ``build_events.py`` CLI: select a single run or a whole
energy point (a group of runs listed in ``data_reference.yml``), choose how to
obtain the calibration, and process each run with the parallel pipeline.
"""

import argparse
import glob
import os
import sys

import yaml

from dataclasses import replace

from siwecal_common import paths

from .calibration import Calibration
from .config import BuilderConfig
from .run_settings import read_threshold_dac, run_settings_path
from .geometry import DetectorGeometry, load_slab_z_mm
from .pad_map import PadMap
from .pipeline import EventBuildingPipeline
from .settings import AppSettings

# ---------------------------------------------------------------- defaults ---
# Filesystem locations come from the shared settings.yml (siwecal_common.paths),
# so no absolute /eos path is hard-coded here. CLI flags and config.yml still
# override these at runtime (see main()).
BASE_PATH = paths.data_dir()
CALIB_DIR_DEFAULT = paths.calib_dir()
CALIB_PEDESTAL_NAME = "dummy_pedestal_15_highgain.txt"
CALIB_MIP_NAME = "dummy_mip_map_15_highgain.txt"
MUON_CALIB_DIR = os.path.join(paths.calib_dir(), "MuonCalib_it2_corrected")
MIP_RUN_DEFAULT = "TB2026CERN_run_000004"
DEFAULT_RUN = "TB2026CERN_run_000013"

_PROJECT_ROOT = paths.REPO_ROOT
# Run list + base paths (energy -> [runs] mapping). Kept separate, unchanged.
CONFIG_FILE = paths.config_file("data_reference_base.yml")
# Optional tunables file (cuts / geometry / calibration). Loaded only if present.
CONFIG_YML_DEFAULT = os.path.join(_PROJECT_ROOT, "config.yml")

# Pad -> (x, y) position maps. A single dict holds the mandatory "default" map
# plus any per-slab override (integer slab key). Overridable via config.yml
# ``mapping:`` (dir + files). Relative paths resolve against PAD_MAP_DIR_DEFAULT.
PAD_MAP_DIR_DEFAULT = paths.geometry_dir()
# Rotated templates (standard maps with x0,y0,x,y negated, i.e. rotated 180 deg
# about the beam axis z). This fixes the hit-position bug: the un-rotated maps
# placed pads/hits in the wrong detector orientation (chip 0 ended up in the
# bottom-left corner instead of its real top-right position), so the
# reconstructed hit_x/hit_y were mirrored through the origin w.r.t. the real
# detector. The rotated maps put every layer back in the true orientation.
PAD_MAP_FILES_DEFAULT = {
    "default": "fev10_rotate_chip_channel_x_y_mapping.txt",
    12: "fev11_cob_good_rotate_chip_channel_x_y_mapping.txt",
}

# Per-slab z positions for hit_z. This file is the live source of truth; if it
# exists it overrides the DetectorGeometry default.
SLAB_Z_FILE_DEFAULT = paths.geometry_file("slab_z_positions.yml")

# Raw DAQ data directory: contains per-run subdirs with Run_Settings.txt.
RAW_BASE_DEFAULT = "/eos/experiment/drdcalo/siw-ecal/TB2026-06/rundata"


# ------------------------------------------------- MuonCalib_it2 resolution ---
def resolve_muon_calib_files(th: str):
    """Return (pedestal_path, mip_path) for the given threshold label (e.g. '220').

    MIP file: the cumulative *_run_000th{th}_highgain.txt in the mips/th{th} folder.
    Pedestal file: the lexicographically latest Pedestal_*_highgain.txt in
    pedestals/th{th} (run numbers are zero-padded so sort order = numeric order).
    """
    mip_path = os.path.join(
        MUON_CALIB_DIR, "mips", f"th{th}",
        f"MIP_pedestalsubmode1_TB2026CERN_run_000th{th}_highgain.txt",
    )
    if not os.path.exists(mip_path):
        print(f"ERROR: cumulative MIP file not found for --th {th}: {mip_path}",
              file=sys.stderr)
        sys.exit(1)

    ped_pattern = os.path.join(
        MUON_CALIB_DIR, "pedestals", f"th{th}", "Pedestal_*_highgain.txt"
    )
    ped_candidates = sorted(glob.glob(ped_pattern))
    if not ped_candidates:
        print(f"ERROR: no pedestal files found for --th {th} in "
              f"{os.path.dirname(ped_pattern)}", file=sys.stderr)
        sys.exit(1)
    pedestal_path = ped_candidates[-1]
    return pedestal_path, mip_path


# ----------------------------------------------------- config / path helpers --
def load_reference_config(path=CONFIG_FILE) -> dict:
    """Read ``data_reference.yml`` (energy -> [runs] mapping and base paths)."""
    if not os.path.exists(path):
        print(f"ERROR: data-reference file not found: {path}\n"
              f"  Pass --data-reference PATH pointing to your "
              f"data_reference YAML.", file=sys.stderr)
        sys.exit(1)
    with open(path) as handle:
        return yaml.safe_load(handle) or {}


def input_path_for(run_name: str, base_path: str) -> str:
    """Absolute path of a run's converted input file."""
    return os.path.join(base_path, run_name, run_name + ".root")


def default_output_for(run_name: str, base_path: str, out_dir=None) -> str:
    """Default output path: ``ecal_<run>.root`` in the run dir or ``out_dir``."""
    file_name = f"ecal_{run_name}.root"
    if out_dir:
        return os.path.join(out_dir, file_name)
    return os.path.join(base_path, run_name, file_name)


def combined_output_for(energy: str, run_names: list, base_path: str,
                        out_dir=None) -> str:
    """Single output path for an energy point, concatenating its run numbers.

    Example: ``P1_20GeV`` with runs 16..20 ->
    ``ecal_P1_20GeV_runs_000016_000017_000018_000019_000020.root``. The file is
    written to ``out_dir`` if given, otherwise to the data ``base_path`` root.
    """
    run_numbers = [name.rsplit("run_", 1)[-1] for name in run_names]
    file_name = f"ecal_{energy}_runs_{'_'.join(run_numbers)}.root"
    return os.path.join(out_dir or base_path, file_name)


def output_for(label: str, run_names: list, base_path: str, out_dir=None,
               explicit=None, is_energy=True) -> str:
    """Resolve the output path for a job (single run or energy point).

    An explicit ``--output`` always wins. Otherwise an energy point with several
    runs gets one combined file named after its runs; a single run (energy point
    or ``--run``) is written as ``ecal_<run>.root`` inside its run directory, so
    the layout matches the ``event_data:`` paths the event display expects.
    """
    if explicit:
        return explicit
    if is_energy and len(run_names) > 1:
        return combined_output_for(label, run_names, base_path, out_dir)
    return default_output_for(run_names[0], base_path, out_dir)


def resolve_available(run_names: list, base_path: str) -> list:
    """Keep only ``(run, path)`` pairs whose converted input file exists."""
    available = []
    for run_name in run_names:
        data_path = input_path_for(run_name, base_path)
        if os.path.exists(data_path):
            available.append((run_name, data_path))
        else:
            print(f"WARNING: input not found, skipping: {data_path}",
                  file=sys.stderr)
    return available


def write_event_data_yaml(event_data: dict, out_path: str) -> None:
    """Write the ``event_data: {label: {path: ...}}`` map, preserving order."""
    with open(out_path, "w") as handle:
        yaml.safe_dump({"event_data": event_data}, handle,
                       default_flow_style=False, sort_keys=False)


# --------------------------------------------------------------- arguments ----
def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="TB2026CERN SiW-ECAL event builder (object-oriented)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  %(prog)s --run TB2026CERN_run_000007\n"
            "  %(prog)s --energy P1_74GeV\n"
            "  %(prog)s --energy P1_20GeV --outdir /tmp/ecal_20gev\n"
            "  %(prog)s --all --data-reference configs/data/data_reference_base.yml\n"
        ),
    )

    target = parser.add_mutually_exclusive_group()
    target.add_argument("--run", default=None, metavar="RUN",
                        help="Single run name (e.g. TB2026CERN_run_000007)")
    target.add_argument("--energy", default=None, metavar="LABEL",
                        help="Energy label from the data-reference YAML; all its "
                             "runs are processed in sequence.")
    target.add_argument("--all", action="store_true",
                        help="Process every entry of the data-reference 'data:' "
                             "map, one after another (each run in parallel).")

    parser.add_argument("--data-reference", default=CONFIG_FILE, metavar="PATH",
                        help="YAML mapping energy label -> [runs] and base paths "
                             f"(default: {CONFIG_FILE}).")
    parser.add_argument("--exclude", nargs="*", default=[], metavar="LABEL",
                        help="Labels to skip in --all mode (e.g. --exclude muons).")
    parser.add_argument("--ref-out", default=None, metavar="PATH",
                        help="Where to write the produced 'event_data:' YAML. "
                             "Defaults in --all mode to <data-reference>_event_data.yml.")

    parser.add_argument("--config", default=CONFIG_YML_DEFAULT, metavar="PATH",
                        help="Optional YAML with cut/geometry/calibration "
                             "overrides. Missing file => dataclass defaults.")

    parser.add_argument("--output", default=None,
                        help="Output file path. In energy mode it sets the name "
                             "of the single combined file.")
    parser.add_argument("--outdir", default=None,
                        help="Output directory (energy mode writes one combined "
                             "file named after its runs here)")

    parser.add_argument("--max-entries", type=int, default=None,
                        help="Max acquisitions per run (debug)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Parallel workers per run")

    parser.add_argument("--th", default=None, metavar="N",
                        help="Threshold label (e.g. 220): loads cumulative "
                             "MuonCalib_it2 calibration for th<N>. "
                             "Overridden by --ped-file / --mip-file.")
    parser.add_argument("--calib-dir", default=None,
                        help="Directory with pedestal/MIP text files")
    parser.add_argument("--ped-file", default=None)
    parser.add_argument("--mip-file", default=None)
    parser.add_argument("--no-calib", action="store_true",
                        help="Raw-ADC mode: no pedestal subtraction, no MIP scale")
    parser.add_argument("--compute-calib", action="store_true",
                        help="Compute pedestals/MIPs from data (slow)")
    parser.add_argument("--ped-run", default=None)
    parser.add_argument("--mip-run", default=None)
    parser.add_argument("--ped-max-entries", type=int, default=None)
    parser.add_argument("--mip-max-entries", type=int, default=None)

    parser.add_argument("--no-mapping", action="store_true",
                        help="Do not assign pad (x,y) positions "
                             "(hit_x/hit_y are written as NaN)")

    parser.add_argument("--raw-dir", default=None, metavar="PATH",
                        help="Base directory containing raw run subdirs with "
                             f"Run_Settings.txt (default: {RAW_BASE_DEFAULT}). "
                             "Used to read ThresholdDAC and store it in the output "
                             "tree.")
    return parser.parse_args(argv)


# --------------------------------------------------------- calibration set-up -
def build_calibration(args, config: BuilderConfig, geometry: DetectorGeometry,
                      base_path: str, calib_dir: str, runs: list,
                      calib_settings: dict = None) -> Calibration:
    """Create the :class:`Calibration` selected by CLI flags and ``config.yml``.

    Precedence per value: CLI flag > ``config.yml`` ``calibration:`` section >
    built-in default. The mode is ``files`` unless ``--no-calib`` (``none``) or
    ``--compute-calib`` (``data``) is given, or ``calibration.mode`` says so.
    """
    calib_settings = calib_settings or {}

    if args.no_calib:
        mode = "none"
    elif args.compute_calib:
        mode = "data"
    else:
        mode = calib_settings.get("mode", "files")
    if mode not in {"files", "data", "none"}:
        print(f"ERROR: invalid calibration mode '{mode}' (files|data|none)",
              file=sys.stderr)
        sys.exit(1)

    if mode == "none":
        print("[Calib] raw ADC mode (no pedestal subtraction, no MIP)")
        return Calibration.disabled(config)

    if mode == "data":
        pedestal_run = args.ped_run or calib_settings.get("pedestal_run") or runs[0]
        mip_run = args.mip_run or calib_settings.get("mip_run") or MIP_RUN_DEFAULT
        pedestal_path = input_path_for(pedestal_run, base_path)
        mip_path = input_path_for(mip_run, base_path)
        print(f"[Calib] Pedestals from {pedestal_run} (computed from data)")
        print(f"[Calib] MIP from       {mip_run} (computed from data)")
        return Calibration.from_data(
            config, geometry, pedestal_path, mip_path,
            pedestal_max_entries=args.ped_max_entries,
            mip_max_entries=args.mip_max_entries)

    if args.th and not (args.ped_file or args.mip_file):
        pedestal_file, mip_file = resolve_muon_calib_files(args.th)
    else:
        pedestal_file = (args.ped_file or calib_settings.get("pedestal_file")
                         or os.path.join(calib_dir, CALIB_PEDESTAL_NAME))
        mip_file = (args.mip_file or calib_settings.get("mip_file")
                    or os.path.join(calib_dir, CALIB_MIP_NAME))
    for path in (pedestal_file, mip_file):
        if not os.path.exists(path):
            print(f"ERROR: calibration file not found: {path}", file=sys.stderr)
            sys.exit(1)
    print(f"[Calib] Pedestals : {pedestal_file}")
    print(f"[Calib] MIP       : {mip_file}")
    return Calibration.from_files(config, pedestal_file, mip_file)


# ----------------------------------------------------------- pad-map set-up --
def build_pad_map(args, mapping_settings: dict = None):
    """Create the :class:`PadMap`, or ``None`` when ``--no-mapping`` is given.

    Precedence per value: ``config.yml`` ``mapping:`` section > built-in default.
    ``files`` from config is merged on top of :data:`PAD_MAP_FILES_DEFAULT`, so a
    config may override just the ``default`` map or add/replace per-slab entries.
    """
    if args.no_mapping:
        print("[PadMap] disabled (--no-mapping): hit_x/hit_y = NaN")
        return None

    mapping_settings = mapping_settings or {}
    files = {**PAD_MAP_FILES_DEFAULT, **(mapping_settings.get("files") or {})}
    base_dir = mapping_settings.get("dir") or PAD_MAP_DIR_DEFAULT
    return PadMap.from_files(files, base_dir)


# ------------------------------------------------------------------- main -----
def resolve_jobs(args, reference: dict) -> list:
    """Build the ordered list of ``(label, run_names, is_energy)`` jobs.

    ``--all`` expands to every entry of the ``data:`` map (minus ``--exclude``);
    ``--energy`` to a single energy point; otherwise a single run.
    """
    data_map = reference.get("data", {})

    if args.all:
        if not data_map:
            print(f"ERROR: no 'data:' mapping in {args.data_reference}",
                  file=sys.stderr)
            sys.exit(1)
        excluded = set(args.exclude or [])
        return [(label, list(runs), True)
                for label, runs in data_map.items() if label not in excluded]

    if args.energy:
        if args.energy not in data_map:
            print(f"ERROR: energy '{args.energy}' not in {args.data_reference}",
                  file=sys.stderr)
            print(f"  Available: {', '.join(data_map.keys())}", file=sys.stderr)
            sys.exit(1)
        return [(args.energy, list(data_map[args.energy]), True)]

    run = args.run or DEFAULT_RUN
    return [(run, [run], False)]


def main(argv=None) -> None:
    args = parse_args(argv)
    reference = load_reference_config(args.data_reference)

    # Optional config.yml layer: defaults <- config.yml (CLI is applied below).
    try:
        settings = AppSettings.from_yaml(args.config)
    except ValueError as error:
        print(f"ERROR in config file '{args.config}': {error}", file=sys.stderr)
        sys.exit(1)
    if args.config and os.path.exists(args.config):
        print(f"[Config] Loaded overrides from {args.config}")

    # base_path: data_reference.yml -> config.yml[paths] -> default.
    base_path = (settings.paths.get("main_path")
                 or reference.get("main_path", BASE_PATH)).rstrip("/")

    jobs = resolve_jobs(args, reference)

    # Shared objects: one config, geometry, calibration and pad-map for all jobs.
    # config/geometry already carry the config.yml overrides (or pure defaults).
    config = settings.config
    geometry = settings.geometry
    # hit_z: slab_z_positions.yml is the live source of truth; fall back to the
    # geometry default if it is absent or empty.
    if os.path.exists(SLAB_Z_FILE_DEFAULT):
        slab_z = load_slab_z_mm(SLAB_Z_FILE_DEFAULT)
        if slab_z:
            geometry = replace(geometry, slab_z_mm=slab_z)
            print(f"[Geometry] hit_z from {SLAB_Z_FILE_DEFAULT} "
                  f"({len(slab_z)} slabs)")
    # calib_dir: CLI -> config.yml[paths] -> data_reference.yml -> default.
    calib_dir = (args.calib_dir
                 or settings.paths.get("calibration_dir")
                 or reference.get("pedestrial_and_mip", CALIB_DIR_DEFAULT)).rstrip("/")
    all_runs = [run for _label, runs, _is_energy in jobs for run in runs]

    # Read ThresholdDAC for every run upfront so we can auto-select calibration.
    raw_dir = args.raw_dir or RAW_BASE_DEFAULT
    all_runs_dacs = {}
    for run_name in all_runs:
        all_runs_dacs[run_name] = read_threshold_dac(
            run_settings_path(raw_dir, run_name))

    # Auto-select --th from Run_Settings when no explicit calibration is forced.
    # If --th / --ped-file / --mip-file / --no-calib / --compute-calib are given
    # the user's choice takes full precedence; otherwise we pick the threshold
    # that appears in the majority of runs and check that the MuonCalib files
    # exist for it before committing.
    if not (args.th or args.ped_file or args.mip_file
            or args.no_calib or args.compute_calib):
        valid_dacs = [v for v in all_runs_dacs.values() if v != -1]
        if valid_dacs:
            auto_th = str(max(set(valid_dacs), key=valid_dacs.count))
            th_ped_dir = os.path.join(MUON_CALIB_DIR, "pedestals", f"th{auto_th}")
            if os.path.isdir(th_ped_dir):
                print(f"[Calib] Auto-detected ThresholdDAC={auto_th} from Run_Settings "
                      f"-> MuonCalib_it2_corrected/th{auto_th} "
                      f"(override with --th or --ped-file/--mip-file)")
                args.th = auto_th
            else:
                print(f"ERROR: ThresholdDAC={auto_th} detected from Run_Settings but no "
                      f"MuonCalib pedestal directory found for that threshold.\n"
                      f"  Looked for: {th_ped_dir}\n"
                      f"  To use a different threshold: --th N\n"
                      f"  To use explicit calibration files: --ped-file / --mip-file\n"
                      f"  To skip calibration entirely: --no-calib",
                      file=sys.stderr)
                sys.exit(1)

    calibration = build_calibration(args, config, geometry, base_path, calib_dir,
                                    all_runs, settings.calibration)
    pad_map = build_pad_map(args, settings.mapping)
    pipeline = EventBuildingPipeline(config, geometry, calibration, pad_map)

    if args.all and args.output:
        print("WARNING: --output is ignored in --all mode (one file per energy)",
              file=sys.stderr)

    event_data = {}
    n_ok = n_fail = 0
    for label, run_names_all, is_energy in jobs:
        available = resolve_available(run_names_all, base_path)
        if not available:
            print(f"SKIP {label}: no input files found.", file=sys.stderr)
            n_fail += 1
            continue

        run_names = [name for name, _path in available]
        input_paths = [path for _name, path in available]
        explicit = None if args.all else args.output
        output_path = output_for(label, run_names, base_path, out_dir=args.outdir,
                                 explicit=explicit, is_energy=is_energy)

        # Reuse the ThresholdDAC values already read above (stored in output tree).
        threshold_dacs = []
        for run_name in run_names:
            dac = all_runs_dacs.get(run_name, -1)
            if dac == -1:
                print(f"  [Settings] WARNING: ThresholdDAC not found for {run_name} "
                      f"(stored as -1). Check --raw-dir.", file=sys.stderr)
            else:
                print(f"  [Settings] {run_name}: ThresholdDAC = {dac}")
            threshold_dacs.append(dac)

        print(f"\n{'=' * 60}")
        if is_energy:
            print(f"  Energy  : {label}  ({len(run_names)} run(s) -> 1 file)")
            print(f"  Runs    : {', '.join(run_names)}")
        else:
            print(f"  Run     : {run_names[0]}")
        print(f"  Output  : {output_path}")
        print(f"  Workers : {args.workers or config.default_workers}")
        print(f"{'=' * 60}")

        os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
        try:
            pipeline.build_runs(input_paths, output_path,
                                max_entries=args.max_entries, n_workers=args.workers,
                                threshold_dacs=threshold_dacs)
        except Exception as error:                       # keep going in --all mode
            print(f"FAIL {label}: {error}", file=sys.stderr)
            n_fail += 1
            if not args.all:
                raise
            continue

        event_data[label] = {"path": os.path.relpath(output_path, base_path)}
        n_ok += 1
        print(f"[Done] {label} -> {output_path}")

    # Emit the event_data YAML for the event display (always in --all mode).
    if event_data and (args.all or args.ref_out):
        ref_out = (args.ref_out
                   or os.path.splitext(args.data_reference)[0] + "_event_data.yml")
        write_event_data_yaml(event_data, ref_out)
        print(f"[Reference] Wrote event_data for {len(event_data)} entr(ies) "
              f"-> {ref_out}")

    print(f"\n[Summary] ok={n_ok} failed={n_fail} of {len(jobs)} job(s)")
    if n_fail and n_ok == 0:
        sys.exit(1)
