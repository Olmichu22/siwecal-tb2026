"""
Structured output layout for the validation plots and results.

Everything is written under a single external base directory (never inside the
event directories). Layout::

    <base>/
      <label>/<category>/<file>.png        # per-run / per-energy plots
      summary/<file>_id<NN>.png            # aggregate (all-energies) plots
      results_id<NN>.csv / results_id<NN>.txt   # selection / fit summary table

The aggregate plots and the results table carry a per-invocation run id
(``id01``, ``id02``, ...) so successive runs (e.g. with different cuts) do not
overwrite each other. The id is the highest one already present in ``<base>``
plus one; the matching ``results_id<NN>.txt`` documents the cut used per energy.
"""

import os
import re


class OutputLayout:
    """Builds and creates the output directory structure on demand."""

    #: matches the ``_id<NN>`` token used by aggregate plots / results files.
    _ID_RE = re.compile(r"_id(\d+)\b")

    def __init__(self, base: str):
        self.base = base

    def label_dir(self, label: str, category: str) -> str:
        """``<base>/<label>/<category>`` (created if missing)."""
        directory = os.path.join(self.base, label, category)
        os.makedirs(directory, exist_ok=True)
        return directory

    def summary_dir(self) -> str:
        """``<base>/summary`` (created if missing)."""
        directory = os.path.join(self.base, "summary")
        os.makedirs(directory, exist_ok=True)
        return directory

    # -------------------------------------------------------------- run id ---
    @staticmethod
    def id_token(run_id: int) -> str:
        """File-name token for a run id, e.g. ``5 -> 'id05'``."""
        return f"id{run_id:02d}"

    def allocate_run_id(self) -> int:
        """Return the next free run id: highest ``_id<NN>`` in ``<base>`` + 1.

        Scans the results files directly under ``<base>`` and the aggregate
        plots under ``<base>/summary``. Returns ``1`` when nothing is present.
        """
        max_id = 0
        for directory in (self.base, os.path.join(self.base, "summary")):
            if not os.path.isdir(directory):
                continue
            for name in os.listdir(directory):
                match = self._ID_RE.search(name)
                if match:
                    max_id = max(max_id, int(match.group(1)))
        return max_id + 1

    def summary_path(self, stem: str, run_id: int) -> str:
        """``<base>/summary/<stem>_id<NN>.png`` (dir created if missing)."""
        return os.path.join(self.summary_dir(),
                            f"{stem}_{self.id_token(run_id)}.png")

    def results_path(self, name: str, run_id: int = None) -> str:
        """Path of a results file under ``<base>`` (created if missing).

        When ``run_id`` is given, the id token is inserted before the extension
        (``results.csv -> results_id<NN>.csv``).
        """
        os.makedirs(self.base, exist_ok=True)
        if run_id is not None:
            root, ext = os.path.splitext(name)
            name = f"{root}_{self.id_token(run_id)}{ext}"
        return os.path.join(self.base, name)
