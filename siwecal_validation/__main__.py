"""
Entry point so the validation runs as ``python -m siwecal_validation``.

All the logic lives in :mod:`siwecal_validation.cli`; this module only wires it
to the ``-m`` launcher. Examples::

    python -m siwecal_validation --run TB2026CERN_run_000007
    python -m siwecal_validation --all
    python -m siwecal_validation --file FILE.root --nhit-min 20 --cache-dir /tmp/cache
"""

from .cli import main

if __name__ == "__main__":
    main()
