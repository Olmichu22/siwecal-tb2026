"""
Results table writer for the validation run.

Collects one row per processed sample and writes both a machine-readable
``results.csv`` and a human-readable ``results.txt``. The key quantity is the
**signal rate surviving the cuts** (``signal_rate = n_selected / n_total``).
"""

import csv

# Column order for both the CSV and the aligned text table.
COLUMNS = [
    ("label", "label"),
    ("energy_gev", "E[GeV]"),
    ("n_total", "N_total"),
    ("n_selected", "N_kept"),
    ("signal_rate", "kept_frac"),
    ("mu_fit", "mu"),
    ("sigma_fit", "sigma"),
    ("resolution", "res"),
    ("shower_frac", "shower_frac"),
    ("mean_weighte", "wE_mean"),
    ("mean_moliere", "mol_mean"),
    ("mean_shower_start", "shStart_mean"),
    ("cuts", "cuts"),
]


class ResultsWriter:
    """Accumulates per-sample result rows and writes them out."""

    def __init__(self):
        self._rows = []

    def add(self, **row):
        """Record one sample. Missing columns are filled with ``None``."""
        if "resolution" not in row or row.get("resolution") is None:
            mu, sigma = row.get("mu_fit"), row.get("sigma_fit")
            row["resolution"] = (sigma / mu) if (mu not in (None, 0)
                                                 and sigma is not None) else None
        self._rows.append({key: row.get(key) for key, _ in COLUMNS})

    @property
    def rows(self):
        return list(self._rows)

    # ------------------------------------------------------------- writers --
    @staticmethod
    def _fmt(value):
        if value is None:
            return ""
        if isinstance(value, float):
            return f"{value:.5g}"
        return str(value)

    def write_csv(self, path):
        with open(path, "w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([key for key, _ in COLUMNS])
            for row in self._rows:
                writer.writerow([self._fmt(row[key]) for key, _ in COLUMNS])
        print(f"Saved results: {path}")

    def write_txt(self, path):
        headers = [head for _, head in COLUMNS]
        table = [headers]
        table += [[self._fmt(row[key]) for key, _ in COLUMNS] for row in self._rows]
        widths = [max(len(r[i]) for r in table) for i in range(len(headers))]
        with open(path, "w") as handle:
            for r, row in enumerate(table):
                handle.write("  ".join(cell.ljust(widths[i])
                                       for i, cell in enumerate(row)).rstrip() + "\n")
                if r == 0:
                    handle.write("  ".join("-" * widths[i]
                                           for i in range(len(headers))) + "\n")
        print(f"Saved results: {path}")
