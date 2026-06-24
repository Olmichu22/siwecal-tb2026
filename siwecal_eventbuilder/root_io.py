"""
ROOT input/output for the event builder.

This module is the *only* place that talks to PyROOT directly. Everything else
in the package works with plain Python objects (:class:`Acquisition`,
:class:`~siwecal_eventbuilder.models.ReconstructedEvent`, ...). Isolating ROOT
here keeps the physics logic testable and easy to read.

Two responsibilities live here:

* :class:`AcquisitionReader` / :class:`Acquisition` -- read the converted input
  tree one acquisition (a.k.a. *spill* / *cycle*) at a time.
* :class:`EcalWriter` -- write reconstructed events into the output ``ecal`` tree.
"""

import numpy as np
import ROOT

from .geometry import DetectorGeometry


class Acquisition:
    """A read-only view of a single acquisition currently loaded in the TTree.

    A SKIROC2 acquisition (also called a *spill* or *readout cycle*) holds the
    buffered data of every chip since the previous readout. The TTree is a
    stateful cursor: only the most recently loaded acquisition is valid, so an
    ``Acquisition`` is a thin, short-lived accessor created by
    :meth:`AcquisitionReader.load`.

    All accessors take structured ``(slab, chip, sca[, channel])`` coordinates
    and translate them to flat array indices via the shared geometry.
    """

    def __init__(self, tree, leaves: dict, geometry: DetectorGeometry, index: int):
        self._tree = tree            # the underlying ROOT TTree (already positioned)
        self._leaves = leaves        # name -> TLeaf, for the large flat arrays
        self._geometry = geometry
        self.index = index           # global acquisition index in the file

    # -- small per-(slab[,chip]) quantities: read via direct branch access ----

    @property
    def n_slab_positions(self) -> int:
        """Number of slab slots present in this acquisition (always 15 here)."""
        return int(self._tree.n_slboards)

    def slboard_id(self, slab: int) -> int:
        """Hardware ID of the slab in slot ``slab`` (``-1`` if the slot is empty)."""
        return int(self._tree.slboard_id[slab])

    def chip_id(self, slab: int, chip: int) -> int:
        """Hardware chip ID (``-999`` for chips that are not present/configured)."""
        return int(self._tree.chipid[slab][chip])

    # -- per-SCA quantities: read via leaves and flat indices for speed -------

    def n_hits(self, slab: int, chip: int, sca: int) -> int:
        """Number of channels that fired in this SCA (0 marks an empty cell)."""
        return int(self._leaves["nhits"].GetValue(
            self._geometry.sca_index(slab, chip, sca)))

    def badbcid(self, slab: int, chip: int, sca: int) -> int:
        """DAQ quality flag for this SCA (0 = good, >0 = retrigger/bad)."""
        return int(self._leaves["badbcid"].GetValue(
            self._geometry.sca_index(slab, chip, sca)))

    def raw_bcid(self, slab: int, chip: int, sca: int) -> int:
        """Raw (uncorrected) bunch-crossing ID stamped on this SCA."""
        return int(self._leaves["bcid"].GetValue(
            self._geometry.sca_index(slab, chip, sca)))

    # -- per-channel quantities ---------------------------------------------

    def hitbit_high(self, slab: int, chip: int, sca: int, channel: int) -> int:
        """High-gain discriminator bit: >0 means the channel triggered."""
        return int(self._leaves["hitbit_high"].GetValue(
            self._geometry.channel_index(slab, chip, sca, channel)))

    def adc_high(self, slab: int, chip: int, sca: int, channel: int) -> int:
        """High-gain ADC value of this channel."""
        return int(self._leaves["adc_high"].GetValue(
            self._geometry.channel_index(slab, chip, sca, channel)))

    def adc_low(self, slab: int, chip: int, sca: int, channel: int) -> int:
        """Low-gain ADC value of this channel."""
        return int(self._leaves["adc_low"].GetValue(
            self._geometry.channel_index(slab, chip, sca, channel)))

    # -- bulk read used by the overflow correction ---------------------------

    def raw_bcid_matrix(self) -> np.ndarray:
        """Return all raw BCIDs as an ``(n_chip_rows, n_scas)`` integer matrix.

        Row ``slab * n_chips_per_slab + chip`` holds the 15 SCA BCIDs of one
        chip. This mirrors the ``reshape(-1, n_scas)`` the reference performs and
        is the input to the per-chip overflow-cycle detection.
        """
        n_boards = self.n_slab_positions
        n_chips = self._geometry.n_chips_per_slab
        n_scas = self._geometry.n_scas_per_chip
        leaf = self._leaves["bcid"]
        sca_index = self._geometry.sca_index

        matrix = np.empty((n_boards * n_chips, n_scas), dtype=np.int64)
        for slab in range(n_boards):
            for chip in range(n_chips):
                row = slab * n_chips + chip
                for sca in range(n_scas):
                    matrix[row, sca] = int(leaf.GetValue(sca_index(slab, chip, sca)))
        return matrix


class AcquisitionReader:
    """Opens a converted input file and yields its acquisitions on demand."""

    #: Leaf names of the large flat arrays read through ``TLeaf.GetValue``.
    _LEAF_NAMES = ("nhits", "bcid", "badbcid", "hitbit_high", "adc_high", "adc_low")

    def __init__(self, input_path: str, geometry: DetectorGeometry, tree_name: str):
        self._geometry = geometry
        self._file = ROOT.TFile(input_path, "READ")
        if self._file.IsZombie():
            raise IOError(f"Could not open input file: {input_path}")
        self._tree = self._file.Get(tree_name)
        if not self._tree:
            raise IOError(f"Tree '{tree_name}' not found in {input_path}")
        self._leaves = {name: self._tree.GetLeaf(name) for name in self._LEAF_NAMES}

    @property
    def n_acquisitions(self) -> int:
        """Total number of acquisitions stored in the file."""
        return int(self._tree.GetEntries())

    def load(self, index: int) -> Acquisition:
        """Position the cursor on acquisition ``index`` and return a view of it."""
        self._tree.GetEntry(index)
        return Acquisition(self._tree, self._leaves, self._geometry, index)

    def close(self) -> None:
        self._file.Close()

    def __enter__(self) -> "AcquisitionReader":
        return self

    def __exit__(self, *_exc) -> None:
        self.close()


class EcalWriter:
    """Creates the output ``ecal`` tree and appends reconstructed events.

    The numpy scalar/array buffers must stay alive for as long as the tree is
    open because ROOT stores their *addresses* as branch buffers; they are kept
    as instance attributes for exactly that reason.
    """

    def __init__(self, output_path: str, max_hits_per_event: int, run_id: int = -1):
        self._file = ROOT.TFile(output_path, "RECREATE")
        self._tree = ROOT.TTree("ecal", "Reconstructed SiW-ECAL events")
        self._max_hits = max_hits_per_event
        self._run_id = int(run_id)

        # Scalar (per-event) branch buffers.
        self._run = np.zeros(1, dtype=np.int32)
        self._event = np.zeros(1, dtype=np.int32)
        self._spill = np.zeros(1, dtype=np.int32)
        self._bcid = np.zeros(1, dtype=np.int32)
        self._n_slab = np.zeros(1, dtype=np.int32)
        self._n_chip = np.zeros(1, dtype=np.int32)
        self._n_chan = np.zeros(1, dtype=np.int32)
        self._sum_hg = np.zeros(1, dtype=np.float32)
        self._sum_energy = np.zeros(1, dtype=np.float32)

        # Variable-length (per-hit) branch buffers, sized by ``nhit_chan``.
        self._hit_slab = np.zeros(max_hits_per_event, dtype=np.int32)
        self._hit_chip = np.zeros(max_hits_per_event, dtype=np.int32)
        self._hit_chan = np.zeros(max_hits_per_event, dtype=np.int32)
        self._hit_sca = np.zeros(max_hits_per_event, dtype=np.int32)
        self._hit_hg = np.zeros(max_hits_per_event, dtype=np.float32)
        self._hit_lg = np.zeros(max_hits_per_event, dtype=np.float32)
        self._hit_energy = np.zeros(max_hits_per_event, dtype=np.float32)
        self._hit_x = np.zeros(max_hits_per_event, dtype=np.float32)
        self._hit_y = np.zeros(max_hits_per_event, dtype=np.float32)
        self._hit_z = np.zeros(max_hits_per_event, dtype=np.float32)

        self._tree.Branch("run", self._run, "run/I")
        self._tree.Branch("event", self._event, "event/I")
        self._tree.Branch("spill", self._spill, "spill/I")
        self._tree.Branch("bcid", self._bcid, "bcid/I")
        self._tree.Branch("nhit_slab", self._n_slab, "nhit_slab/I")
        self._tree.Branch("nhit_chip", self._n_chip, "nhit_chip/I")
        self._tree.Branch("nhit_chan", self._n_chan, "nhit_chan/I")
        self._tree.Branch("sum_hg", self._sum_hg, "sum_hg/F")
        self._tree.Branch("sum_energy", self._sum_energy, "sum_energy/F")
        self._tree.Branch("hit_slab", self._hit_slab, "hit_slab[nhit_chan]/I")
        self._tree.Branch("hit_chip", self._hit_chip, "hit_chip[nhit_chan]/I")
        self._tree.Branch("hit_chan", self._hit_chan, "hit_chan[nhit_chan]/I")
        self._tree.Branch("hit_sca", self._hit_sca, "hit_sca[nhit_chan]/I")
        self._tree.Branch("hit_hg", self._hit_hg, "hit_hg[nhit_chan]/F")
        self._tree.Branch("hit_lg", self._hit_lg, "hit_lg[nhit_chan]/F")
        self._tree.Branch("hit_energy", self._hit_energy, "hit_energy[nhit_chan]/F")
        self._tree.Branch("hit_x", self._hit_x, "hit_x[nhit_chan]/F")
        self._tree.Branch("hit_y", self._hit_y, "hit_y[nhit_chan]/F")
        self._tree.Branch("hit_z", self._hit_z, "hit_z[nhit_chan]/F")

    def write(self, event, spill_index: int, event_index: int) -> bool:
        """Fill one :class:`ReconstructedEvent`. Returns ``False`` if skipped.

        ``spill_index`` is the global acquisition index; ``event_index`` is the
        running index of this event within its acquisition. They are combined
        into the ``event`` number (unique within a run). The ``run`` branch
        records which run the event came from, so ``(run, event)`` is globally
        traceable even after several runs are merged into one file.
        """
        if not event.hits or event.n_channels > self._max_hits:
            return False

        self._run[0] = self._run_id
        self._spill[0] = spill_index
        self._event[0] = spill_index * 1000 + event_index
        self._bcid[0] = event.bcid
        self._n_chan[0] = event.n_channels
        self._n_slab[0] = event.n_slabs
        self._n_chip[0] = event.n_chips
        self._sum_hg[0] = event.sum_adc_high
        self._sum_energy[0] = event.sum_energy

        for i, hit in enumerate(event.hits):
            self._hit_slab[i] = hit.slab_position
            self._hit_chip[i] = hit.chip_id
            self._hit_chan[i] = hit.channel
            self._hit_sca[i] = hit.sca
            self._hit_hg[i] = hit.adc_high_pedsub
            self._hit_lg[i] = hit.adc_low_pedsub
            self._hit_energy[i] = hit.energy_mip
            self._hit_x[i] = hit.x
            self._hit_y[i] = hit.y
            self._hit_z[i] = hit.z

        self._tree.Fill()
        return True

    def close(self) -> None:
        self._file.Write()
        self._file.Close()
