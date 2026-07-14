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

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"

#include <limits>
#include <algorithm>
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
    // InputFiles (a list) or InputFile (one path) -- a run that was decoded in
    // parallel exists only as per-chunk files, and CONVERT's chunking is a
    // parallelisation detail, not a property of the data. Chaining them here is
    // what the calibrator's InputFiles already does, and it means a run never
    // has to be merged into one file first: TB2026CERN_run_000004's 1369 chunks
    // come to ~176 GB, which crosses ROOT's 100 GB TTree::fgMaxTreeSize and so
    // cannot even be held in a single tree. The OUTPUT is one `ecal` file either
    // way -- only the input side is chained.
    std::vector<std::string> inputs = m_inputFiles.value();
    if (inputs.empty() && !m_inputFile.value().empty()) inputs.push_back(m_inputFile.value());
    if (inputs.empty()) {
      return error() << "Set InputFiles (list) or InputFile (single path)" << endmsg, StatusCode::FAILURE;
    }
    if (m_outputFile.value().empty()) {
      return error() << "OutputFile not set" << endmsg, StatusCode::FAILURE;
    }

    // TChain IS a TTree, so everything downstream (bindReadBranches, GetEntry,
    // GetEntries) works unchanged; the chain just rolls over files for us.
    auto chain = std::make_unique<TChain>(m_treeName.value().c_str());
    for (const auto& path : inputs) {
      if (chain->AddFile(path.c_str()) == 0) {
        return error() << "Cannot add input file (missing, or holds no '" << m_treeName.value()
                       << "' tree): " << path << endmsg,
               StatusCode::FAILURE;
      }
    }
    if (chain->GetEntries() <= 0) {
      return error() << "No entries in '" << m_treeName.value() << "' across the " << inputs.size()
                     << " input file(s)" << endmsg,
             StatusCode::FAILURE;
    }
    if (inputs.size() > 1) {
      info() << "Chained " << inputs.size() << " input file(s), " << chain->GetEntries() << " acquisitions"
             << endmsg;
    }
    TTree* tree = chain.get();

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
               << " take their energy from the low-gain branch, anchored to the high gain through "
                  "adc_low - ped_lg = "
               << m_gainRatio.value() << "*(adc_high - ped_hg) + " << m_gainIntercept.value()
               << " (MIP_lg is NOT used; hit_energy_nocalib keeps the old MIP_lg-based value for comparison)"
               << endmsg;
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
    cfg.gainRatio = m_gainRatio.value();
    cfg.gainIntercept = m_gainIntercept.value();
    cfg.maxHitsPerSca = m_maxHitsPerSca.value();
    cfg.maxHitsPerEvent = m_maxHitsPerEvent.value();

    const k4siwecal::EventBuilder builder(cfg, calib, geom, padMap.get());

    int runId = m_runId.value();
    // Chunk files are named chunk_NNNN.root, so the run number is only in their
    // directory -- parse the full path, not the basename.
    if (runId < 0) runId = parseRunIdFromPath(inputs.front());

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

    // A run's raw chunks are decoded independently, so an acquisition that spans
    // a raw-chunk boundary comes out TWICE: once at the end of chunk N and once
    // at the start of chunk N+1, each copy holding only the part of the readout
    // that fell in that file. Measured on run_000060: 55 of 28,860 acquisitions
    // (0.19%), and every one of them a pair of CONSECUTIVE chain entries.
    //
    // The two copies are COMPLEMENTARY, not redundant: in all 55 pairs their
    // populated (slab, chip, SCA) cells are entirely disjoint, so splicing them
    // back together recovers the whole readout window and loses nothing. Dropping
    // either copy would throw away real hits (~0.29% of a run's) -- and dropping
    // "the second one" would be worse still, since the fragment is the FIRST of
    // the pair in 34 of those 55 (acq 1558 is 4 hits followed by 35).
    //
    // The pre-pass reads only acqNumber, so it costs one cheap scan of the chain.
    std::vector<char> mergeWithNext(static_cast<std::size_t>(nEntries), 0);
    std::vector<char> absorbed(static_cast<std::size_t>(nEntries), 0);
    long long spliced = 0;
    {
      tree->SetBranchStatus("*", false);
      tree->SetBranchStatus("acqNumber", true);
      int prevAcq = -1;
      for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        if (entry > 0 && buf->acqNumber == prevAcq && !absorbed[static_cast<std::size_t>(entry - 1)]) {
          mergeWithNext[static_cast<std::size_t>(entry - 1)] = 1;
          absorbed[static_cast<std::size_t>(entry)] = 1;
          ++spliced;
        }
        prevAcq = buf->acqNumber;
      }
      tree->SetBranchStatus("*", true);
    }
    if (spliced > 0) {
      info() << "EcalEventBuilder: spliced " << spliced
             << " acquisition(s) split across a raw-chunk boundary (both halves kept)" << endmsg;
    }

    // Only materialised when an acquisition actually needs splicing (a few dozen
    // per run): TreeBuffers is a couple of MB, far too big to copy per entry.
    auto splice = std::make_unique<TreeBuffers>();

    long long totalEvents = 0;
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
      if (absorbed[static_cast<std::size_t>(entry)]) {
        continue;  // already folded into the previous entry
      }
      tree->GetEntry(entry);

      const TreeBuffers* src = buf.get();
      if (mergeWithNext[static_cast<std::size_t>(entry)]) {
        *splice = *buf;
        tree->GetEntry(entry + 1);
        mergeAcquisition(*splice, *buf);
        src = splice.get();
      }

      const k4siwecal::AcquisitionView acq(src->nSlboards, src->slboardId, src->chipId, src->bcid, src->nhits,
                                            src->hitbitHigh, src->adcHigh, src->adcLow);
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
  Gaudi::Property<std::string> m_inputFile{this, "InputFile", "", "siwecaldecoded ROOT file (single)"};
  Gaudi::Property<std::vector<std::string>> m_inputFiles{
      this, "InputFiles", {},
      "siwecaldecoded ROOT files to chain, in order. Takes precedence over InputFile. Lets a run be event-built "
      "straight from its decoded chunks, with no merge step."};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "siwecaldecoded", "Input TTree name"};
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "", "Output ecal ROOT file"};
  Gaudi::Property<int> m_runId{this, "RunId", -1, "-1 = parse from the first input path (..._run_<N>...)"};
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
      this, "AdcSaturationThreshold", 1500,
      "Raw (not pedestal-subtracted) high-gain ADC at/above which the high-gain preamp is taken to be "
      "saturated and the energy is read from the low-gain branch instead. Only has an effect when the "
      "low-gain tables are loaded."};
  Gaudi::Property<double> m_gainRatio{
      this, "GainRatio", 0.0962,
      "k in adc_low - ped_lg = k*(adc_high - ped_hg) + c. Above AdcSaturationThreshold the hit energy is "
      "the low-gain reading anchored back onto the high gain through this line, so MIP_lg is never used "
      "(it is unmeasurable: a low-gain MIP is 2.08 ADC against ~1 ADC of noise). Measured 0.0974 on "
      "run_000012/th230 and 0.0967 on run_000072/th220 by tracing the band with medians."};
  Gaudi::Property<double> m_gainIntercept{
      this, "GainIntercept", 1.45,
      "c in the same line: the low-gain ADC left over at zero high-gain signal. Real (+1 to +2 ADC, always "
      "positive, seen by two independent robust methods) but only good to ~0.5 ADC -- the error is "
      "systematic, not statistical. It shifts the energy by c/k ~ 0.7 MIP, under 1% at the switch point."};
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
    int acqNumber = -1;
    int nSlboards = 0;
    int slboardId[kSlbDepth];
    int chipId[kSlbDepth][kSkirocsPerAsu];
    int bcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
    int nhits[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
    int hitbitHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
    int adcHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
    int adcLow[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  };

  /// Splice the other half of a boundary-split acquisition into `dst`.
  ///
  /// The two halves are disjoint (verified on every one of run_000060's 55 pairs:
  /// no populated (slab, chip, SCA) cell appears in both), so this is a union, not
  /// a sum -- a cell is taken from `src` only where `dst` has nothing. Written that
  /// way rather than as a blind add so that an overlap, if one ever occurs, keeps
  /// `dst`'s value instead of silently double-counting the hits.
  static void mergeAcquisition(TreeBuffers& dst, const TreeBuffers& src) {
    dst.nSlboards = std::max(dst.nSlboards, src.nSlboards);
    for (int slb = 0; slb < kSlbDepth; ++slb) {
      if (dst.slboardId[slb] < 0 && src.slboardId[slb] >= 0) {
        dst.slboardId[slb] = src.slboardId[slb];
      }
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        if (dst.chipId[slb][chip] < 0 && src.chipId[slb][chip] >= 0) {
          dst.chipId[slb][chip] = src.chipId[slb][chip];
        }
        for (int sca = 0; sca < kScasInSkiroc; ++sca) {
          if (dst.nhits[slb][chip][sca] > 0 || src.nhits[slb][chip][sca] <= 0) {
            continue;  // dst already holds this cell, or src has nothing to give
          }
          dst.nhits[slb][chip][sca] = src.nhits[slb][chip][sca];
          dst.bcid[slb][chip][sca] = src.bcid[slb][chip][sca];
          for (int ch = 0; ch < kChannelsInSkiroc; ++ch) {
            dst.hitbitHigh[slb][chip][sca][ch] = src.hitbitHigh[slb][chip][sca][ch];
            dst.adcHigh[slb][chip][sca][ch] = src.adcHigh[slb][chip][sca][ch];
            dst.adcLow[slb][chip][sca][ch] = src.adcLow[slb][chip][sca][ch];
          }
        }
      }
    }
  }

  static void bindReadBranches(TTree& tree, TreeBuffers& buf) {
    tree.SetBranchAddress("acqNumber", &buf.acqNumber);
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
      m_hitEnergyNoCalib.resize(maxHitsPerEvent);
      m_hitWEnergy.resize(maxHitsPerEvent);
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
      m_tree->Branch("sum_energy_nocalib", &m_sumEnergyNoCalib, "sum_energy_nocalib/F");
      m_tree->Branch("sum_w_energy", &m_sumWEnergy, "sum_w_energy/F");
      m_tree->Branch("hit_slab", m_hitSlab.data(), "hit_slab[nhit_chan]/I");
      m_tree->Branch("hit_chip", m_hitChip.data(), "hit_chip[nhit_chan]/I");
      m_tree->Branch("hit_chan", m_hitChan.data(), "hit_chan[nhit_chan]/I");
      m_tree->Branch("hit_sca", m_hitSca.data(), "hit_sca[nhit_chan]/I");
      m_tree->Branch("hit_hg", m_hitHg.data(), "hit_hg[nhit_chan]/F");
      m_tree->Branch("hit_lg", m_hitLg.data(), "hit_lg[nhit_chan]/F");
      m_tree->Branch("hit_energy", m_hitEnergy.data(), "hit_energy[nhit_chan]/F");
      // The old, MIP_lg-based energy, alongside the new one on the same events: the
      // two differ only for hits above AdcSaturationThreshold, so any difference in
      // a distribution built from them is entirely the saturated hits.
      m_tree->Branch("hit_energy_nocalib", m_hitEnergyNoCalib.data(), "hit_energy_nocalib[nhit_chan]/F");
      m_tree->Branch("hit_w_energy", m_hitWEnergy.data(), "hit_w_energy[nhit_chan]/F");
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
      m_sumEnergyNoCalib = static_cast<float>(event.sumEnergyNoCalib());
      m_sumWEnergy = static_cast<float>(event.sumWEnergy());

      for (int i = 0; i < event.nChannels(); ++i) {
        const auto& hit = event.hits[i];
        m_hitSlab[i] = hit.slabPosition;
        m_hitChip[i] = hit.chipId;
        m_hitChan[i] = hit.channel;
        m_hitSca[i] = hit.sca;
        m_hitHg[i] = hit.adcHighPedsub;
        m_hitLg[i] = hit.adcLowPedsub;
        m_hitEnergy[i] = hit.energyMip;
        m_hitEnergyNoCalib[i] = hit.energyMipNoCalib;
        m_hitWEnergy[i] = hit.wEnergy;
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
    float m_sumHg = 0.f, m_sumEnergy = 0.f, m_sumEnergyNoCalib = 0.f, m_sumWEnergy = 0.f;
    std::vector<int> m_hitSlab, m_hitChip, m_hitChan, m_hitSca, m_hitIsMasked;
    std::vector<float> m_hitHg, m_hitLg, m_hitEnergy, m_hitEnergyNoCalib, m_hitWEnergy, m_hitX, m_hitY,
        m_hitZ, m_hitX0;
  };
};

DECLARE_COMPONENT(EcalEventBuilder)
