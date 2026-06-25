/*
 * EcalToEDM4hep: read the (non-podio) reconstructed `ecal` TTree one event at a
 * time and emit EDM4hep collections. The `ecal` tree is already one row per
 * physics event, so this is a strict 1->1 source: EvtMax is set to the number of
 * tree entries in the steering file and each call produces one event's hits.
 *
 * Output per event:
 *   - CalorimeterHitCollection ECalHits : one hit per fired channel, energy =
 *     hit_energy (MIP units), position = (hit_x, hit_y, hit_z) [mm], cellID
 *     encoding (slab, chip, channel, sca), and `type` = slab (layer) for native
 *     access without a bitfield decoder.
 *   - EventHeaderCollection EventHeader : run/event (+ bcid in timeStamp, spill
 *     in weight).
 *   - UserDataCollections (parallel to ECalHits, same order): ECalHitChip,
 *     ECalHitChan, ECalHitSca (int) and ECalHitHG, ECalHitLG (float). These
 *     carry the per-hit quantities EDM4hep CalorimeterHit has no native field
 *     for, kept as podio-native collections (no cellID decode needed downstream).
 *
 * NOTE: ROOT TTree reading is not thread-safe; the sequential entry cursor is
 * guarded by a mutex. Run single-threaded for now (one source).
 */
#include "k4FWCore/Producer.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "podio/UserDataCollection.h"

#include "DDSegmentation/BitFieldCoder.h"

#include "Gaudi/Property.h"

#include "TFile.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"

#include <atomic>
#include <memory>
#include <mutex>
#include <string>
#include <tuple>

using OutputType = std::tuple<edm4hep::CalorimeterHitCollection, edm4hep::EventHeaderCollection,
                              podio::UserDataCollection<std::int32_t>,   // chip
                              podio::UserDataCollection<std::int32_t>,   // chan
                              podio::UserDataCollection<std::int32_t>,   // sca
                              podio::UserDataCollection<float>,          // hg
                              podio::UserDataCollection<float>>;         // lg

struct EcalToEDM4hep final : k4FWCore::Producer<OutputType()> {
  EcalToEDM4hep(const std::string& name, ISvcLocator* svcLoc)
      : Producer(name, svcLoc, {},
                 {KeyValues("OutputCaloHits", {"ECalHits"}),
                  KeyValues("OutputEventHeader", {"EventHeader"}),
                  KeyValues("OutputHitChip", {"ECalHitChip"}),
                  KeyValues("OutputHitChan", {"ECalHitChan"}),
                  KeyValues("OutputHitSca", {"ECalHitSca"}),
                  KeyValues("OutputHitHG", {"ECalHitHG"}),
                  KeyValues("OutputHitLG", {"ECalHitLG"})}) {}

  StatusCode initialize() override {
    m_file.reset(TFile::Open(m_inputFile.value().c_str(), "READ"));
    if (!m_file || m_file->IsZombie())
      return error() << "Cannot open input file: " << m_inputFile.value() << endmsg, StatusCode::FAILURE;
    auto* tree = dynamic_cast<TTree*>(m_file->Get(m_treeName.value().c_str()));
    if (!tree)
      return error() << "Tree '" << m_treeName.value() << "' not found" << endmsg, StatusCode::FAILURE;

    m_reader = std::make_unique<TTreeReader>(tree);
    m_nhit = std::make_unique<TTreeReaderValue<int>>(*m_reader, "nhit_chan");
    m_run = std::make_unique<TTreeReaderValue<int>>(*m_reader, "run");
    m_event = std::make_unique<TTreeReaderValue<int>>(*m_reader, "event");
    m_spill = std::make_unique<TTreeReaderValue<int>>(*m_reader, "spill");
    m_bcid = std::make_unique<TTreeReaderValue<int>>(*m_reader, "bcid");
    m_hitSlab = std::make_unique<TTreeReaderArray<int>>(*m_reader, "hit_slab");
    m_hitChip = std::make_unique<TTreeReaderArray<int>>(*m_reader, "hit_chip");
    m_hitChan = std::make_unique<TTreeReaderArray<int>>(*m_reader, "hit_chan");
    m_hitSca = std::make_unique<TTreeReaderArray<int>>(*m_reader, "hit_sca");
    m_hitEnergy = std::make_unique<TTreeReaderArray<float>>(*m_reader, "hit_energy");
    m_hitHG = std::make_unique<TTreeReaderArray<float>>(*m_reader, "hit_hg");
    m_hitLG = std::make_unique<TTreeReaderArray<float>>(*m_reader, "hit_lg");
    m_hitX = std::make_unique<TTreeReaderArray<float>>(*m_reader, "hit_x");
    m_hitY = std::make_unique<TTreeReaderArray<float>>(*m_reader, "hit_y");
    m_hitZ = std::make_unique<TTreeReaderArray<float>>(*m_reader, "hit_z");
    // Masked (uncalibrated) channels are flagged by the event builder and
    // dropped here so the downstream PID sees only calibrated hits. Older ecal
    // files predate the branch; bind it only when present (else keep all hits).
    if (tree->GetBranch("hit_ismasked"))
      m_hitMasked = std::make_unique<TTreeReaderArray<int>>(*m_reader, "hit_ismasked");

    m_coder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_cellIDEncoding.value());
    info() << "EcalToEDM4hep: " << m_reader->GetEntries(false) << " entries in '"
           << m_treeName.value() << "'" << endmsg;
    return StatusCode::SUCCESS;
  }

  OutputType operator()() const override {
    edm4hep::CalorimeterHitCollection hits;
    edm4hep::EventHeaderCollection headers;
    podio::UserDataCollection<std::int32_t> chip, chan, sca;
    podio::UserDataCollection<float> hg, lg;

    std::scoped_lock lock(m_mutex);
    const Long64_t entry = m_cursor.fetch_add(1);
    if (m_reader->SetEntry(entry) != TTreeReader::kEntryValid) {
      warning() << "Failed to read entry " << entry << endmsg;
      return std::make_tuple(std::move(hits), std::move(headers), std::move(chip), std::move(chan),
                             std::move(sca), std::move(hg), std::move(lg));
    }

    auto header = headers.create();
    header.setRunNumber(static_cast<std::uint32_t>(**m_run));
    header.setEventNumber(static_cast<std::uint64_t>(**m_event));
    header.setTimeStamp(static_cast<std::uint64_t>(**m_bcid));
    header.setWeight(static_cast<double>(**m_spill));

    const int n = **m_nhit;
    for (int i = 0; i < n; ++i) {
      // Skip masked channels: dropped from ECalHits and all parallel collections
      // at once, so they stay index-aligned and the PID recomputes on the rest.
      if (m_hitMasked && (*m_hitMasked)[i])
        continue;
      auto hit = hits.create();
      std::uint64_t cellID = 0;
      m_coder->set(cellID, "slab", (*m_hitSlab)[i]);
      m_coder->set(cellID, "chip", (*m_hitChip)[i]);
      m_coder->set(cellID, "channel", (*m_hitChan)[i]);
      m_coder->set(cellID, "sca", (*m_hitSca)[i]);
      hit.setCellID(cellID);
      hit.setEnergy((*m_hitEnergy)[i]);
      hit.setPosition({(*m_hitX)[i], (*m_hitY)[i], (*m_hitZ)[i]});
      // The full coordinate lives in cellID; the layer is also stored in `type`
      // so consumers can read it natively (hit.getType()) without a decoder.
      hit.setType((*m_hitSlab)[i]);
      // Per-hit quantities without a native CalorimeterHit field, as parallel
      // podio collections (same order as the hits).
      chip.push_back((*m_hitChip)[i]);
      chan.push_back((*m_hitChan)[i]);
      sca.push_back((*m_hitSca)[i]);
      hg.push_back((*m_hitHG)[i]);
      lg.push_back((*m_hitLG)[i]);
    }
    return std::make_tuple(std::move(hits), std::move(headers), std::move(chip), std::move(chan),
                           std::move(sca), std::move(hg), std::move(lg));
  }

  Gaudi::Property<std::string> m_inputFile{this, "InputFile", "", "Path to the ecal ROOT file"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "ecal", "Input TTree name"};
  Gaudi::Property<std::string> m_cellIDEncoding{
      this, "CellIDEncoding", "system:8,slab:8,chip:16,channel:8,sca:8",
      "DD4hep bitfield descriptor for the CalorimeterHit cellID"};

private:
  std::unique_ptr<TFile> m_file;
  std::unique_ptr<TTreeReader> m_reader;
  std::unique_ptr<TTreeReaderValue<int>> m_nhit, m_run, m_event, m_spill, m_bcid;
  std::unique_ptr<TTreeReaderArray<int>> m_hitSlab, m_hitChip, m_hitChan, m_hitSca;
  std::unique_ptr<TTreeReaderArray<int>> m_hitMasked;   // null for pre-mask ecal files
  std::unique_ptr<TTreeReaderArray<float>> m_hitEnergy, m_hitHG, m_hitLG, m_hitX, m_hitY, m_hitZ;
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_coder;
  mutable std::mutex m_mutex;
  mutable std::atomic<Long64_t> m_cursor{0};
};

DECLARE_COMPONENT(EcalToEDM4hep)
