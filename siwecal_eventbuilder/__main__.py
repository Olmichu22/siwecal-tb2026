"""
Entry point so the event builder runs as ``python -m siwecal_eventbuilder``.

All the logic lives in :mod:`siwecal_eventbuilder.cli`; this module only wires it
to the ``-m`` launcher. Examples::

    python -m siwecal_eventbuilder --run TB2026CERN_run_000007
    python -m siwecal_eventbuilder --energy P1_20GeV --outdir /tmp/ecal_20gev
    python -m siwecal_eventbuilder --all --data-reference configs/data/data_reference_base.yml
"""

from .cli import main

if __name__ == "__main__":
    main()
