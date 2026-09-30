#!/usr/bin/env python
"""
Write the valtree (valcache-schema ``ecal`` TTree) of already-produced PID EDM4hep files, without re-running the
PID stage: every event of ``ecal_<label>.edm4hep.root`` -> ``ecal_<label>.valtree.root`` next to it (or in
--outdir). Same writer as run_pid_batch.py --format valtree (siwecal_common.edm4hep_pid.write_valtree).

    python gaudi_jobs/pid_to_valtree.py <file.edm4hep.root> [...] [--outdir DIR] [--force]
    python gaudi_jobs/pid_to_valtree.py --campaign /eos/.../Reconstructed_final      # every run in it

The MIP-variant blocks (mip05_/mip1_) are written only if the file carries them (validation mode).
"""
import argparse
import glob
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(_HERE))

from siwecal_common import edm4hep_pid  # noqa: E402
from siwecal_validation.config import PlotConfig  # noqa: E402


def convert(path, outdir=None, force=False, block=40000):
    label = os.path.basename(path)[len("ecal_"):-len(".edm4hep.root")]
    out = os.path.join(outdir or os.path.dirname(path), f"ecal_{label}.valtree.root")
    if os.path.exists(out) and not force:
        print(f"[skip] {out} exists")
        return out
    config = PlotConfig()
    reader = edm4hep_pid.PidFileReader(path, n_layers=config.n_layers, mip_thresholds=(0.5, 1.0))
    mip = (0.5, 1.0) if any(n.startswith("mip05_") for n in reader.shape_names) else ()
    if not mip:
        reader = edm4hep_pid.PidFileReader(path, n_layers=config.n_layers, mip_thresholds=())
    n = len(reader.identifiers()["event"])
    # Blocks of events, then hadd: write_valtree holds its events' hits in memory (~45 kB/event), which for the
    # 475k-event runs would be 20 GB.
    if n <= block:
        written = edm4hep_pid.write_valtree(reader, out, config, range(n), mip_thresholds=mip)
    else:
        import subprocess, tempfile
        tmp = tempfile.mkdtemp(prefix="valtree_", dir=os.environ.get("TMPDIR", "/tmp"))
        parts, written = [], 0
        for i, start in enumerate(range(0, n, block)):
            part = os.path.join(tmp, f"part_{i:03d}.root")
            written += edm4hep_pid.write_valtree(reader, part, config, range(start, min(start + block, n)),
                                                 mip_thresholds=mip)
            parts.append(part)
        local = os.path.join(tmp, "merged.root")
        subprocess.run(["hadd", "-f", "-v", "0", local] + parts, check=True)
        subprocess.run(["cp", "-f", local, out], check=True)
        subprocess.run(["rm", "-rf", tmp], check=True)
    print(f"[valtree] {written} event(s) -> {out}")
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="*")
    p.add_argument("--campaign", default=None, help="Convert <campaign>/<run>/ecal_<run>.edm4hep.root for every run")
    p.add_argument("--outdir", default=None)
    p.add_argument("--force", action="store_true")
    a = p.parse_args(argv)
    files = list(a.files)
    if a.campaign:
        files += sorted(glob.glob(os.path.join(a.campaign, "*", "ecal_*.edm4hep.root")))
    files = [f for f in files if f.endswith(".edm4hep.root") and "_moliere_" not in f and "_weighte_" not in f]
    for f in files:
        convert(f, a.outdir, a.force)


if __name__ == "__main__":
    main()
