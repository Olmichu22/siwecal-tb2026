"""
Parser for the DAQ ``Run_Settings.txt`` files produced alongside each raw run.

The settings file records per-chip hardware configuration, including the
discriminator threshold ``ThresholdDAC`` that determines which channels fire.
This value is needed later for pedestal calibration (calibration tables are
organised per threshold, e.g. ``MuonCalib_it2_corrected/pedestals/th230/``).
"""

import os
import re
from collections import Counter

_THRESHOLD_RE = re.compile(r"ThresholdDAC:\s*(\d+)")


def read_threshold_dac(settings_path: str) -> int:
    """Return the nominal ThresholdDAC from a Run_Settings.txt file.

    Reads every ``ThresholdDAC: N`` token in the file and returns the mode
    (most common value). Per-chip fine-tuning entries that appear during
    calibration scans are thus outvoted by the bulk of data-taking lines.

    Returns -1 when the file cannot be opened or contains no ThresholdDAC token.
    """
    values = []
    try:
        with open(settings_path) as fh:
            for line in fh:
                match = _THRESHOLD_RE.search(line)
                if match:
                    values.append(int(match.group(1)))
    except OSError:
        return -1
    if not values:
        return -1
    return Counter(values).most_common(1)[0][0]


def run_settings_path(raw_dir: str, run_name: str) -> str:
    """Absolute path to the Run_Settings.txt for *run_name* inside *raw_dir*."""
    return os.path.join(raw_dir, run_name, "Run_Settings.txt")
