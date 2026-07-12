/*
 * EcalEventBuilder: event building as a Gaudi component, porting
 * siwecal_eventbuilder (BCID clustering, calibration application, geometry
 * mapping) to C++. Reads the `siwecaldecoded` tree written by EcalRawDecoder
 * and pre-computed pedestal/MIP calibration text tables written by
 * PedestalMipCalibrator (from a DIFFERENT data-taking run -- pedestal and
 * MIP calibration are generated upstream, not as part of this component),
 * and writes the `ecal` tree in the exact schema
 * siwecal_eventbuilder/root_io.py::EcalWriter already produces, so the
 * existing EcalToEDM4hep component reads it unchanged.
 *
 * `threshold_dac` is read directly from the `thresholdDac` branch EcalRawDecoder
 * writes into `siwecaldecoded` (see RunSettings.h / EcalRawDecoder.cpp) --
 * this component does not re-parse Run_Settings.txt itself.
 *
 * Plain Gaudi::Algorithm (not a k4FWCore Producer/Transformer), same
 * reasoning as EcalRawDecoder/PedestalMipCalibrator: this is a global
 * reduction over an entire run (the number of reconstructed events isn't
 * known ahead of time -- BCID clustering can yield zero or several events
 * per acquisition), and the output is a legacy ROOT tree, not an EDM4hep
 * collection. All work happens in initialize(); execute() is a no-op.
 */
#include "k4SiWEcalReco/CalibrationTables.h"
#include "k4SiWEcalReco/EventBuilder.h"
#include "k4SiWEcalReco/PadMapGeometry.h"
#include "k4SiWEcalReco/SlbFrameDecoder.h"  // kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include "TFile.h"
#include "TTree.h"

#include <limits>
#include <map>
#include <memory>
#include <regex>
#include <string>
#include <vector>

using k4siwecal::kChannelsInSkiroc;
using k4siwecal::kScasInSkiroc;
using k4siwecal::kSkirocsPerAsu;
using k4siwecal::kSlbDepth;

namespace {
// Port of siwecal_eventbuilder/pipeline.py::_run_id_from_path: extracts the
// numeric run id from a "..._run_<N>..." filename. Only used as a fallback
// when the RunId property is left at its -1 default.
int parseRunIdFromPath(const std::string& path) {
  static const std::regex re(R"(run_(\d+))");
  std::smatch m;
  if (std::regex_search(path, m, re)) {
    return std::stoi(m[1].str());
  }
  return -1;
}
}  // namespace

struct EcalEventBuilder final : Gaudi::Algorithm {
  using Gaudi::Algorithm::Algorithm;

  StatusCode initialize() override {
    if (m_inputFile.value().empty()) {
      return error() << "InputFile not set" << endmsg, StatusCode::FAILURE;
    }
    if (m_outputFile.value().empty()) {
      return error() << "OutputFile not set" << endmsg, StatusCode::FAILURE;
    }

    std::unique_ptr<TFile> fin(TFile::Open(m_inputFile.value().c_str(), "READ"));
    if (!fin || fin->IsZombie()) {
      return error() << "Cannot open input file: " << m_inputFile.value() << endmsg, StatusCode::FAILURE;
    }
    auto* tree = dynamic_cast<TTree*>(fin->Get(m_treeName.value().c_str()));
    if (!tree) {
      return error() << "Tree '" << m_treeName.value() << "' not found in " << m_inputFile.value() << endmsg,
             StatusCode::FAILURE;
    }

    // TreeBuffers is ~2.8MB (three [15][16][15][64] int arrays) -- too large
    // for a stack frame on the smaller worker-thread stacks Gaudi/k4run can
    // use (this already caused a real segfault in EcalRawDecoder before its
    // Acquisition buffer was moved to the heap). Heap-allocate explicitly.
    auto buf = std::make_unique<TreeBuffers>();
    bindReadBranches(*tree, *buf);

    k4siwecal::CalibrationTables calib;
    if (m_noCalibration.value()) {
      calib = k4siwecal::CalibrationTables::disabled();
    } else {
      if (m_pedestalFile.value().empty() || m_mipFile.value().empty()) {
        return error() << "PedestalFile and MipFile must both be set unless NoCalibration=true" << endmsg,
               StatusCode::FAILURE;
      }
      // Low-gain tables are optional but must be given as a PAIR: one without
      // the other would silently reconstruct saturated hits against a missing
      // pedestal or MIP scale.
      const bool hasPedLg = !m_pedestalFileLowGain.value().empty();
      const bool hasMipLg = !m_mipFileLowGain.value().empty();
      if (hasPedLg != hasMipLg) {
        return error() << "PedestalFileLowGain and MipFileLowGain must be set together (or both left empty "
                          "to disable the low-gain saturation recovery)"
                       << endmsg,
               StatusCode::FAILURE;
      }
      bool ok = false;
      std::string calibError;
      calib = k4siwecal::CalibrationTables::fromFiles(m_pedestalFile.value(), m_mipFile.value(),
                                                       m_pedestalFallback.value(), m_defaultMipFallback.value(), ok,
                                                       calibError, m_pedestalFileLowGain.value(),
                                                       m_mipFileLowGain.value());
      if (!ok) {
        return error() << calibError << endmsg, StatusCode::FAILURE;
      }
      if (calib.hasLowGain()) {
        info() << "Low-gain calibration loaded: hits with raw adc_high >= " << m_adcSaturationThreshold.value()
               << " take their energy from the low-gain branch" << endmsg;
      } else {
        warning() << "No low-gain calibration given: saturated high-gain hits will be under-read "
                     "(set PedestalFileLowGain/MipFileLowGain to enable saturation recovery)"
                  << endmsg;
      }
    }

    const k4siwecal::SlabGeometry geom = m_slabZFile.value().empty()
                                             ? k4siwecal::SlabGeometry::defaults()
                                             : k4siwecal::SlabGeometry::fromYamlFile(m_slabZFile.value());

    std::unique_ptr<k4siwecal::PadMap> padMap;
    if (!m_noMapping.value() && !m_padMapDefaultFile.value().empty()) {
      std::map<int, std::string> overrides;
      for (const auto& entry : m_padMapSlabOverrides.value()) {
        const auto colon = entry.find(':');
        if (colon == std::string::npos) {
          warning() << "Ignoring malformed PadMapSlabOverrides entry (expected 'slab:path'): " << entry << endmsg;
          continue;
        }
        overrides[std::stoi(entry.substr(0, colon))] = entry.substr(colon + 1);
      }
      padMap = std::make_unique<k4siwecal::PadMap>(k4siwecal::PadMap::fromFiles(m_padMapDefaultFile.value(), overrides));
    }

    k4siwecal::BuilderConfig cfg;
    cfg.skipBcidStart = m_skipBcidStart.value();
    cfg.dropBcids.clear();
    for (int v : m_dropBcids.value()) cfg.dropBcids.insert(v);
    cfg.mergeDelta = m_mergeDelta.value();
    cfg.minSlabsHit = m_minSlabsHit.value();
    cfg.bcidOverflow = m_bcidOverflow.value();
    cfg.badValue = m_badValue.value();
    cfg.adcUnderflowThreshold = m_adcUnderflowThreshold.value();
    cfg.adcSaturationThreshold = m_adcSaturationThreshold.value();
    cfg.maxHitsPerSca = m_maxHitsPerSca.value();
    cfg.maxHitsPerEvent = m_maxHitsPerEvent.value();

    const k4siwecal::EventBuilder builder(cfg, calib, geom, padMap.get());

    int runId = m_runId.value();
    if (runId < 0) runId = parseRunIdFromPath(m_inputFile.value());

    int thresholdDac = m_thresholdDacOverride.value();
    if (thresholdDac < 0 && tree->GetBranch("thresholdDac") != nullptr && tree->GetEntries() > 0) {
      int tdBuf = -1;
      tree->SetBranchAddress("thresholdDac", &tdBuf);
      tree->GetEntry(0);
      thresholdDac = tdBuf;
      tree->ResetBranchAddress(tree->GetBranch("thresholdDac"));
    }

    std::unique_ptr<TFile> fout(TFile::Open(m_outputFile.value().c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
      return error() << "Cannot create output file: " << m_outputFile.value() << endmsg, StatusCode::FAILURE;
    }
    fout->cd();
    EcalTreeWriter writer(cfg.maxHitsPerEvent, runId, thresholdDac);

    const Long64_t nEntries = tree->GetEntries();
    long long totalEvents = 0;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
      tree->GetEntry(entry);
      const k4siwecal::AcquisitionView acq(buf->nSlboards, buf->slboardId, buf->chipId, buf->bcid, buf->nhits,
                                            buf->hitbitHigh, buf->adcHigh, buf->adcLow);
      auto events = builder.build(acq);
      for (int eventIndex = 0; eventIndex < static_cast<int>(events.size()); ++eventIndex) {
        writer.write(events[eventIndex], static_cast<int>(entry), eventIndex);
      }
      totalEvents += static_cast<long long>(events.size());
    }

    fout->cd();
    fout->Write(nullptr, TObject::kOverwrite);
    fout->Close();
    info() << "EcalEventBuilder: wrote " << totalEvents << " event(s) from " << nEntries
           << " acquisition(s) to " << m_outputFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override { return StatusCode::SUCCESS; }

  // I/O
  Gaudi::Property<std::string> m_inputFile{this, "InputFile", "", "siwecaldecoded ROOT file"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "siwecaldecoded", "Input TTree name"};
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "", "Output ecal ROOT file"};
  Gaudi::Property<int> m_runId{this, "RunId", -1, "-1 = parse from InputFile basename (..._run_<N>...)"};
  Gaudi::Property<int> m_thresholdDacOverride{
      this, "ThresholdDacOverride", -1,
      "-1 = read from siwecaldecoded's thresholdDac branch (written by EcalRawDecoder from Run_Settings.txt); "
      "any other value forces that value / covers siwecaldecoded files without the branch"};

  // Calibration
  Gaudi::Property<std::string> m_pedestalFile{this, "PedestalFile", "", "Pedestal calibration text table"};
  Gaudi::Property<std::string> m_mipFile{this, "MipFile", "", "MIP calibration text table"};
  Gaudi::Property<std::string> m_pedestalFileLowGain{
      this, "PedestalFileLowGain", "",
      "Low-gain pedestal table. Optional, but must be set together with MipFileLowGain; both empty disables "
      "the low-gain saturation recovery and hit energy always comes from the high gain."};
  Gaudi::Property<std::string> m_mipFileLowGain{this, "MipFileLowGain", "", "Low-gain MIP calibration text table"};
  Gaudi::Property<bool> m_noCalibration{this, "NoCalibration", false, "Raw-ADC mode: pedestal=0, mip=1"};
  Gaudi::Property<double> m_pedestalFallback{this, "PedestalFallback", 250.0, ""};
  Gaudi::Property<double> m_defaultMipFallback{this, "DefaultMipFallback", 20.0, ""};

  // Geometry / mapping
  Gaudi::Property<std::string> m_padMapDefaultFile{this, "PadMapDefaultFile", "", "Default chip/channel->x,y map"};
  Gaudi::Property<std::vector<std::string>> m_padMapSlabOverrides{
      this, "PadMapSlabOverrides", {}, "Per-slab override maps as 'slab:path' entries, e.g. '12:mappings/...txt'"};
  Gaudi::Property<std::string> m_slabZFile{this, "SlabZFile", "",
                                           "slab_z_positions.yml; empty = compiled-in TB2026CERN defaults"};
  Gaudi::Property<bool> m_noMapping{this, "NoMapping", false, "Skip pad mapping: hit_x/hit_y stay NaN"};

  // BuilderConfig (matching siwecal_eventbuilder/config.py::BuilderConfig defaults)
  Gaudi::Property<int> m_skipBcidStart{this, "SkipBcidStart", 50, ""};
  Gaudi::Property<std::vector<int>> m_dropBcids{this, "DropBcids", {0, 901}, ""};
  Gaudi::Property<int> m_mergeDelta{this, "MergeDelta", 3, ""};
  Gaudi::Property<int> m_minSlabsHit{this, "MinSlabsHit", 10, ""};
  Gaudi::Property<int> m_bcidOverflow{this, "BcidOverflow", 4096, ""};
  Gaudi::Property<int> m_badValue{this, "BadValue", -999, ""};
  Gaudi::Property<int> m_adcUnderflowThreshold{this, "AdcUnderflowThreshold", 11, ""};
  Gaudi::Property<int> m_adcSaturationThreshold{
      this, "AdcSaturationThreshold", 1200,
      "Raw (not pedestal-subtracted) high-gain ADC at/above which the high-gain preamp is taken to be "
      "saturated and the energy is read from the low-gain branch instead. Only has an effect when the "
      "low-gain tables are loaded."};
  // config.py's BuilderConfig.max_hits_per_sca defaults to math.inf ("disabled"),
  // but Gaudi's confdb2 generator can't serialize an infinite double property
  // default (emits invalid Python "inf.0" and fails the build) -- use a huge
  // finite sentinel instead, which is behaviourally identical for this cut
  // (per-SCA hit counts realistically top out around 64).
  Gaudi::Property<double> m_maxHitsPerSca{this, "MaxHitsPerSca", 1.0e18, "math.inf equivalent (disabled)"};
  Gaudi::Property<int> m_maxHitsPerEvent{this, "MaxHitsPerEvent", 15360, "15 slabs * 16 chips * 64 channels"};

 private:
  // Bound branch buffers for the siwecaldecoded tree, matching the layout
  // written by EcalRawDecoder. Only the branches BcidClusterer/HitCollector
  // actually use are bound (badbcid/autogainbit_*/hitbit_low are read by the
  // Python AcquisitionReader but never referenced by the clustering/hit
  // collection algorithms themselves).
  struct TreeBuffers {
    int nSlboards = 0;
    int slboardId[kSlbDepth];
    int chipId[kSlbDepth][kSkirocsPerAsu];
    int bcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
    int nhits[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
    int hitbitHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
    int adcHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
    int adcLow[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  };

  static void bindReadBranches(TTree& tree, TreeBuffers& buf) {
    tree.SetBranchAddress("n_slboards", &buf.nSlboards);
    tree.SetBranchAddress("slboard_id", buf.slboardId);
    tree.SetBranchAddress("chipid", buf.chipId);
    tree.SetBranchAddress("bcid", buf.bcid);
    tree.SetBranchAddress("nhits", buf.nhits);
    tree.SetBranchAddress("hitbit_high", buf.hitbitHigh);
    tree.SetBranchAddress("adc_high", buf.adcHigh);
    tree.SetBranchAddress("adc_low", buf.adcLow);
  }

  // Port of siwecal_eventbuilder/root_io.py::EcalWriter -- branch-for-branch
  // identical schema (names, ROOT leaf-type strings, `[nhit_chan]` sizing),
  // since this feeds the EXISTING EcalToEDM4hep component unchanged.
  struct EcalTreeWriter {
    EcalTreeWriter(int maxHitsPerEvent, int runId, int thresholdDac)
        : m_maxHits(maxHitsPerEvent), m_run(runId), m_thresholdDac(thresholdDac) {
      m_hitSlab.resize(maxHitsPerEvent);
      m_hitChip.resize(maxHitsPerEvent);
      m_hitChan.resize(maxHitsPerEvent);
      m_hitSca.resize(maxHitsPerEvent);
      m_hitIsMasked.resize(maxHitsPerEvent);
      m_hitHg.resize(maxHitsPerEvent);
      m_hitLg.resize(maxHitsPerEvent);
      m_hitEnergy.resize(maxHitsPerEvent);
      m_hitX.resize(maxHitsPerEvent);
      m_hitY.resize(maxHitsPerEvent);
      m_hitZ.resize(maxHitsPerEvent);
      m_hitX0.resize(maxHitsPerEvent);

      m_tree = new TTree("ecal", "Reconstructed SiW-ECAL events");
      m_tree->Branch("run", &m_run, "run/I");
      m_tree->Branch("threshold_dac", &m_thresholdDac, "threshold_dac/I");
      m_tree->Branch("event", &m_event, "event/I");
      m_tree->Branch("spill", &m_spill, "spill/I");
      m_tree->Branch("bcid", &m_bcid, "bcid/I");
      m_tree->Branch("nhit_slab", &m_nSlab, "nhit_slab/I");
      m_tree->Branch("nhit_chip", &m_nChip, "nhit_chip/I");
      m_tree->Branch("nhit_chan", &m_nChan, "nhit_chan/I");
      m_tree->Branch("sum_hg", &m_sumHg, "sum_hg/F");
      m_tree->Branch("sum_energy", &m_sumEnergy, "sum_energy/F");
      m_tree->Branch("hit_slab", m_hitSlab.data(), "hit_slab[nhit_chan]/I");
      m_tree->Branch("hit_chip", m_hitChip.data(), "hit_chip[nhit_chan]/I");
      m_tree->Branch("hit_chan", m_hitChan.data(), "hit_chan[nhit_chan]/I");
      m_tree->Branch("hit_sca", m_hitSca.data(), "hit_sca[nhit_chan]/I");
      m_tree->Branch("hit_hg", m_hitHg.data(), "hit_hg[nhit_chan]/F");
      m_tree->Branch("hit_lg", m_hitLg.data(), "hit_lg[nhit_chan]/F");
      m_tree->Branch("hit_energy", m_hitEnergy.data(), "hit_energy[nhit_chan]/F");
      m_tree->Branch("hit_x", m_hitX.data(), "hit_x[nhit_chan]/F");
      m_tree->Branch("hit_y", m_hitY.data(), "hit_y[nhit_chan]/F");
      m_tree->Branch("hit_z", m_hitZ.data(), "hit_z[nhit_chan]/F");
      m_tree->Branch("hit_X0", m_hitX0.data(), "hit_X0[nhit_chan]/F");
      m_tree->Branch("hit_ismasked", m_hitIsMasked.data(), "hit_ismasked[nhit_chan]/I");
    }

    // Port of EcalWriter.write. Returns false (skip) if the event is empty or
    // exceeds the fixed per-hit buffer capacity, matching root_io.py:235.
    bool write(const k4siwecal::ReconstructedEvent& event, int spillIndex, int eventIndex) {
      if (event.hits.empty() || event.nChannels() > m_maxHits) return false;

      m_spill = spillIndex;
      m_event = spillIndex * 1000 + eventIndex;
      m_bcid = static_cast<int>(event.bcid);
      m_nChan = event.nChannels();
      m_nSlab = event.nSlabs();
      m_nChip = event.nChips();
      m_sumHg = static_cast<float>(event.sumAdcHigh());
      m_sumEnergy = static_cast<float>(event.sumEnergy());

      for (int i = 0; i < event.nChannels(); ++i) {
        const auto& hit = event.hits[i];
        m_hitSlab[i] = hit.slabPosition;
        m_hitChip[i] = hit.chipId;
        m_hitChan[i] = hit.channel;
        m_hitSca[i] = hit.sca;
        m_hitHg[i] = hit.adcHighPedsub;
        m_hitLg[i] = hit.adcLowPedsub;
        m_hitEnergy[i] = hit.energyMip;
        m_hitX[i] = hit.x;
        m_hitY[i] = hit.y;
        m_hitZ[i] = hit.z;
        m_hitX0[i] = hit.x0;
        m_hitIsMasked[i] = hit.isMasked ? 1 : 0;
      }
      m_tree->Fill();
      return true;
    }

    int m_maxHits;
    TTree* m_tree = nullptr;
    int m_run = -1, m_thresholdDac = -1, m_event = 0, m_spill = 0, m_bcid = 0;
    int m_nSlab = 0, m_nChip = 0, m_nChan = 0;
    float m_sumHg = 0.f, m_sumEnergy = 0.f;
    std::vector<int> m_hitSlab, m_hitChip, m_hitChan, m_hitSca, m_hitIsMasked;
    std::vector<float> m_hitHg, m_hitLg, m_hitEnergy, m_hitX, m_hitY, m_hitZ, m_hitX0;
  };
};

DECLARE_COMPONENT(EcalEventBuilder)
