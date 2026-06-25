#!/usr/bin/env python3
"""Clean SiWECAL pedestal tables so they satisfy the calibration policy.

Each data row is ``layer chip channel`` followed by 15 SCA triplets
``mean err width``. A ``-nan`` only ever appears in the ``mean``/``width``
sub-columns (``err`` carries a finite sentinel ``-5``/``-10``). A mean is a
real pedestal only if it is finite and positive; ``0``, ``-nan`` and negative
values all mark an SCA without a usable fit. This tool rewrites every row so
the event-builder loader accepts it (each pad is either all-zero -> masked, or
all finite/in-range -> calibrated):

1. **No usable mean** (all SCA means NaN, zero or negative) -> the pad is dead;
   the row becomes 15 ``0 -10 0`` triplets (all-zero means -> the loader masks
   it).
2. **At least one usable mean** -> every non-usable mean (NaN/0/negative) is
   replaced by the *robust* mean of the usable ones (a 2-sigma clip drops
   outliers from that average, but each usable SCA keeps its own value),
   truncated to the same number of decimals as the usable means. NaN widths are
   likewise filled from the row's valid widths.
3. **Final range check** -> any resulting mean below MIN_MEAN or above MAX_MEAN
   is replaced by DEFAULT_MEAN. Dead pads (all-zero) are exempt.

Already-clean rows (no NaN, and means either all-zero or all in-range) are
passed through byte-for-byte. Cleaned files keep their original basename and
are written into a mandatory output directory (no default), leaving the inputs
untouched.

Usage:
    python clean_pedestals.py --out-dir DIR file1.txt [file2.txt ...]
"""

import argparse
import math
import os
import statistics


N_SCA = 15
MIN_MEAN = 100.0          # below this a mean is non-physical
MAX_MEAN = 300.0          # above this a mean is non-physical
DEFAULT_MEAN = 240.0      # fallback for out-of-range means


def _is_nan(token: str) -> bool:
    """True for ``nan``/``-nan`` tokens (case-insensitive)."""
    return token.lower().lstrip("+-") == "nan"


def _is_usable_mean(token: str) -> bool:
    """True for a mean that is a real pedestal: finite and strictly positive."""
    return not _is_nan(token) and float(token) > 0.0


def _decimals(token: str) -> int:
    """Number of decimal places written in ``token`` (0 if integer)."""
    return len(token.split(".")[1]) if "." in token else 0


def _truncate(value: float, ndec: int) -> str:
    """Format ``value`` truncated (not rounded) to ``ndec`` decimals."""
    factor = 10 ** ndec
    truncated = math.trunc(value * factor) / factor
    return f"{truncated:.{ndec}f}"


def _robust_mean(values: list) -> float:
    """Mean of ``values`` after dropping points beyond 2 sigma (outliers)."""
    if len(values) < 3:
        return statistics.fmean(values)
    mean, std = statistics.fmean(values), statistics.pstdev(values)
    if std == 0.0:
        return mean
    kept = [v for v in values if abs(v - mean) <= 2 * std]
    return statistics.fmean(kept) if kept else mean


def _fill_nan(tokens: list) -> list:
    """Replace NaN tokens by the truncated mean of the non-NaN ones."""
    valid = [t for t in tokens if not _is_nan(t)]
    if not valid or not any(_is_nan(t) for t in tokens):
        return tokens
    mean = sum(float(t) for t in valid) / len(valid)
    ndec = max(_decimals(t) for t in valid)
    filled = _truncate(mean, ndec)
    return [filled if _is_nan(t) else t for t in tokens]


def _clean_mean(token: str, fill: str, ndec: int) -> str:
    """Resolve one mean to a usable, in-range value (final range check applied).

    A usable, in-range mean keeps its original token; an out-of-range usable
    mean becomes DEFAULT_MEAN; a non-usable mean takes ``fill`` (already clamped
    by the caller).
    """
    if _is_usable_mean(token):
        value = float(token)
        if MIN_MEAN <= value <= MAX_MEAN:
            return token
        return _truncate(DEFAULT_MEAN, ndec)
    return fill


def clean_line(line: str) -> tuple:
    """Return ``(cleaned_line, category)`` for one input line.

    ``category`` is ``""`` (passed through), ``"dead"`` (masked as 0) or
    ``"filled"`` (means/widths repaired).
    """
    if line.startswith("#") or not line.strip():
        return line, ""
    items = line.split()
    header, data = items[:3], items[3:3 + 3 * N_SCA]
    means = data[0::3]
    usable = [m for m in means if _is_usable_mean(m)]
    has_nan = any(_is_nan(t) for t in data)

    # Already clean: no NaN and means either all-zero (dead) or all in-range.
    all_zero = not has_nan and all(float(m) == 0.0 for m in means)
    all_in_range = bool(usable) and all(
        _is_usable_mean(m) and MIN_MEAN <= float(m) <= MAX_MEAN for m in means)
    if not has_nan and (all_zero or all_in_range):
        return line, ""

    # Dead pad: nothing usable to anchor a calibration -> mask with zeros.
    if not usable:
        triplets = " ".join(["0 -10 0"] * N_SCA)
        return f"{' '.join(header)} {triplets}\n", "dead"

    # Calibrated pad: fill non-usable means from the robust mean, range-check.
    ndec = max(_decimals(m) for m in usable)
    reference = _robust_mean([float(m) for m in usable])
    if MIN_MEAN <= reference <= MAX_MEAN:
        fill = _truncate(reference, ndec)
    else:
        fill = _truncate(DEFAULT_MEAN, ndec)
    new_means = [_clean_mean(m, fill, ndec) for m in means]
    new_errs = _fill_nan(data[1::3])
    new_widths = _fill_nan(data[2::3])
    rebuilt = [tok for triplet in zip(new_means, new_errs, new_widths)
               for tok in triplet]
    return f"{' '.join(header + rebuilt)}\n", "filled"


def clean_file(in_path: str, out_dir: str) -> tuple:
    """Clean ``in_path`` into ``out_dir`` (same basename). Returns counts."""
    n_dead = n_filled = 0
    out_path = os.path.join(out_dir, os.path.basename(in_path))
    with open(in_path) as src, open(out_path, "w") as dst:
        for line in src:
            cleaned, category = clean_line(line)
            n_dead += category == "dead"
            n_filled += category == "filled"
            dst.write(cleaned)
    return out_path, n_dead, n_filled


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("files", nargs="+", help="pedestal .txt files to clean")
    parser.add_argument("--out-dir", required=True,
                        help="output directory (mandatory, no default)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    for in_path in args.files:
        out_path, n_dead, n_filled = clean_file(in_path, args.out_dir)
        print(f"{in_path} -> {out_path}  "
              f"(dead pads masked: {n_dead}, pads repaired: {n_filled})")


if __name__ == "__main__":
    main()
