/*
 * EcalRawDecoder: raw2root as a Gaudi component. Decodes the binary SLB DAQ
 * files (*_raw.bin_NNNN chunks) of one run into a `siwecaldecoded` ROOT tree,
 * the same non-podio legacy format already read by
 * siwecal_eventbuilder/root_io.py::AcquisitionReader, replacing the external
 * reference tool SiWECAL-TB-analysis/converter_SLB/SLBraw2ROOT.cc.
 *
 * Not a k4FWCore Producer/Transformer: the number of output tree entries
 * (decoded acquisition cycles) is only known after streaming and
 * cycle-buffering the whole input, and the output is a fixed-layout legacy
 * ROOT tree, not an EDM4hep/podio collection -- neither fits the "known
 * EvtMax, one podio collection per Gaudi event" k4FWCore contract. This is a
 * plain Gaudi::Algorithm that does all its work in initialize() (same
 * self-contained-I/O pattern EcalToEDM4hep already uses for its own
 * TFile/TTreeReader); execute() is a no-op.
 *
 * Chunk-file handling: by default (ResetStatePerInputFile=true) a fresh
 * cycle-buffer is used per input chunk file, matching the reference's
 * per-chunk-file SLBraw2ROOT instance (ConvertDirectorySL_Raw.cc:64) -- a
 * cycle split across a chunk boundary is not reassembled, same as the
 * reference. All chunks' decoded cycles are appended into one shared output
 * tree, replacing the reference's separate-file-per-chunk + hadd step.
 */
#include "k4SiWEcalReco/RunSettings.h"
#include "k4SiWEcalReco/SlbFrameDecoder.h"

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include <fstream>
#include <memory>
#include <string>
#include <vector>

struct EcalRawDecoder final : Gaudi::Algorithm {
  using Gaudi::Algorithm::Algorithm;

  StatusCode initialize() override {
    if (m_inputFiles.value().empty()) {
      return error() << "No InputFiles provided" << endmsg, StatusCode::FAILURE;
    }

    auto fout = std::unique_ptr<TFile>(TFile::Open(m_outputFile.value().c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
      return error() << "Cannot create output file: " << m_outputFile.value() << endmsg, StatusCode::FAILURE;
    }
    fout->cd();
    auto* tree = new TTree(m_treeName.value().c_str(), m_treeName.value().c_str());

    // Acquisition is ~5.5MB (fixed [15][16][15][64] int arrays) -- too large
    // for a stack frame on the smaller worker-thread stacks Gaudi/k4run can
    // use. Heap-allocate the tree-bound buffer explicitly.
    auto acq = std::make_unique<k4siwecal::Acquisition>();
    // Run-level metadata (constant across every entry): not part of
    // Acquisition (which models actual per-cycle decoded hardware state),
    // matching how threshold_dac/run are handled in
    // siwecal_eventbuilder/root_io.py::EcalWriter -- bound once, never
    // touched again, so every Fill() reads the same value.
    k4siwecal::RunSettings runSettings;
    if (!m_runSettingsFile.value().empty()) {
      runSettings = k4siwecal::parseRunSettings(m_runSettingsFile.value());
    }
    bindBranches(*tree, *acq, runSettings);

    long long totalCycles = 0;
    if (m_resetStatePerInputFile.value()) {
      for (const auto& path : m_inputFiles.value()) {
        std::vector<k4siwecal::Acquisition> flushed;
        if (!decodeFileWithRetry(path, flushed)) {
          return error() << "Could not decode " << path << " without repeated cycles after retries" << endmsg,
                 StatusCode::FAILURE;
        }
        for (const auto& a : flushed) {
          *acq = a;
          tree->Fill();
        }
        totalCycles += static_cast<long long>(flushed.size());
        info() << "EcalRawDecoder: decoded " << flushed.size() << " acquisition cycles from " << path << endmsg;
      }
    } else {
      k4siwecal::CycleAssembler assembler(m_maxReadoutCycleJump.value(), m_bcidThreshold.value());
      std::vector<k4siwecal::Acquisition> flushed;
      for (const auto& path : m_inputFiles.value()) {
        if (!streamFile(path, assembler, flushed)) {
          return error() << "Cannot open raw file: " << path << endmsg, StatusCode::FAILURE;
        }
      }
      assembler.drain(flushed);
      if (assembler.repeatedCycleCount() > 0) {
        warning() << "Repeated cycles detected across the combined input; MaxReadoutCycleJump="
                  << m_maxReadoutCycleJump.value()
                  << " may be too small (ResetStatePerInputFile=false disables the automatic retry)" << endmsg;
      }
      for (const auto& a : flushed) {
        *acq = a;
        tree->Fill();
      }
      totalCycles = static_cast<long long>(flushed.size());
    }

    fout->cd();
    fout->Write(nullptr, TObject::kOverwrite);
    fout->Close();
    info() << "EcalRawDecoder: wrote " << totalCycles << " total acquisition cycles to " << m_outputFile.value()
           << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override { return StatusCode::SUCCESS; }

  Gaudi::Property<std::vector<std::string>> m_inputFiles{
      this, "InputFiles", {}, "Ordered list of raw *_raw.bin_NNNN chunk files for one run"};
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "", "Output ROOT file (siwecaldecoded tree)"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "siwecaldecoded", "Output TTree name"};
  Gaudi::Property<int> m_maxReadoutCycleJump{
      this, "MaxReadoutCycleJump", 10,
      "Cycle-buffer flush depth; auto-retried x10 (up to a bounded number of attempts) on detected cycle "
      "repeats when ResetStatePerInputFile=true, matching the legacy converter's retry behaviour"};
  Gaudi::Property<int> m_bcidThreshold{this, "BcidThreshold", k4siwecal::kBcidThresDefault,
                                       "BCIDTHRES: retrigger/empty-event classification window"};
  Gaudi::Property<bool> m_zeroSuppression{this, "ZeroSuppression", false,
                                          "Store only hit channels instead of all 64 per SCA"};
  Gaudi::Property<bool> m_eudaqFormat{this, "EudaqFormat", false,
                                      "Use the 0xABCD EUDAQ frame sync instead of the plain 0xEEEEEEEE one"};
  Gaudi::Property<bool> m_computeBadBcid{this, "ComputeBadBcid", true, "Run the bad-BCID tagging state machine"};
  Gaudi::Property<bool> m_resetStatePerInputFile{
      this, "ResetStatePerInputFile", true,
      "Match legacy behaviour: fresh cycle-buffer state at each input chunk-file boundary"};
  Gaudi::Property<std::string> m_runSettingsFile{
      this, "RunSettingsFile", "",
      "Path to Run_Settings.txt (optional); if empty or unreadable, the new run-settings "
      "branches are written with their -1 sentinel default"};

 private:
  static void bindBranches(TTree& tree, k4siwecal::Acquisition& acq, k4siwecal::RunSettings& runSettings) {
    using k4siwecal::kChannelsInSkiroc;
    using k4siwecal::kScasInSkiroc;
    using k4siwecal::kSkirocsPerAsu;
    using k4siwecal::kSlbDepth;

    tree.Branch("acqNumber", &acq.acqNumber, "acqNumber/I");
    tree.Branch("n_slboards", &acq.nSlboards, "n_slboards/I");
    tree.Branch("slot", acq.slot, Form("slot[%d]/I", kSlbDepth));
    tree.Branch("slboard_id", acq.slboardId, Form("slboard_id[%d]/I", kSlbDepth));
    tree.Branch("chipid", acq.chipId, Form("chipid[%d][%d]/I", kSlbDepth, kSkirocsPerAsu));
    tree.Branch("nColumns", acq.numCol, Form("nColumns[%d][%d]/I", kSlbDepth, kSkirocsPerAsu));
    tree.Branch("startACQ", acq.startAcq, Form("startACQ[%d]/F", kSlbDepth));
    tree.Branch("rawTSD", acq.rawTsd, Form("rawTSD[%d]/I", kSlbDepth));
    tree.Branch("TSD", acq.tsd, Form("TSD[%d]/F", kSlbDepth));
    tree.Branch("rawAVDD0", acq.rawAvdd0, Form("rawAVDD0[%d]/I", kSlbDepth));
    tree.Branch("rawAVDD1", acq.rawAvdd1, Form("rawAVDD1[%d]/I", kSlbDepth));
    tree.Branch("AVDD0", acq.avdd0, Form("AVDD0[%d]/F", kSlbDepth));
    tree.Branch("AVDD1", acq.avdd1, Form("AVDD1[%d]/F", kSlbDepth));
    tree.Branch("bcid", acq.bcid, Form("bcid[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
    tree.Branch("corrected_bcid", acq.correctedBcid,
                Form("corrected_bcid[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
    tree.Branch("badbcid", acq.badbcid, Form("badbcid[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
    tree.Branch("nhits", acq.nhits, Form("nhits[%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc));
    tree.Branch("adc_low", acq.adcLow,
                Form("adc_low[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
    tree.Branch("adc_high", acq.adcHigh,
                Form("adc_high[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
    tree.Branch(
        "autogainbit_low", acq.autogainbitLow,
        Form("autogainbit_low[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
    tree.Branch(
        "autogainbit_high", acq.autogainbitHigh,
        Form("autogainbit_high[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
    tree.Branch("hitbit_low", acq.hitbitLow,
                Form("hitbit_low[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));
    tree.Branch("hitbit_high", acq.hitbitHigh,
                Form("hitbit_high[%d][%d][%d][%d]/I", kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc));

    // Run_Settings.txt-derived metadata (RunSettings.h), same value on every
    // entry -- see siwecal_eventbuilder/run_settings.py::read_threshold_dac
    // for the ThresholdDAC precedent this generalizes.
    tree.Branch("acqWindowMs", &runSettings.acqWindowMs, "acqWindowMs/F");
    tree.Branch("delayBetweenCycleMs", &runSettings.delayBetweenCycleMs, "delayBetweenCycleMs/F");
    tree.Branch("thresholdDac", &runSettings.thresholdDac, "thresholdDac/I");
    tree.Branch("holdDelay", &runSettings.holdDelay, "holdDelay/I");
    tree.Branch("fsPeakTime", &runSettings.fsPeakTime, "fsPeakTime/I");
    tree.Branch("gainSelectionThreshold", &runSettings.gainSelectionThreshold, "gainSelectionThreshold/I");
  }

  // Streams one raw file's frames into `assembler`, appending any flushed
  // acquisitions into `flushed`. Returns false only if the file cannot be
  // opened.
  bool streamFile(const std::string& path, k4siwecal::CycleAssembler& assembler,
                   std::vector<k4siwecal::Acquisition>& flushed) const {
    std::ifstream fin(path, std::ios::binary);
    if (!fin.is_open()) return false;
    k4siwecal::RawFrameReader reader(fin, m_eudaqFormat.value());
    std::vector<unsigned char> frameBytes;
    while (reader.nextFrame(frameBytes)) {
      assembler.addFrame(frameBytes, flushed);
    }
    return true;
  }

  // Decodes one chunk file with a fresh CycleAssembler, escalating
  // MaxReadoutCycleJump by x10 and re-reading from scratch (bounded to a
  // handful of attempts) if repeated cycles are detected, mirroring
  // ConvertDirectorySL_Raw.cc:66-71's `while(result==false)` retry loop.
  bool decodeFileWithRetry(const std::string& path, std::vector<k4siwecal::Acquisition>& flushed) const {
    int maxJump = m_maxReadoutCycleJump.value();
    constexpr int kMaxAttempts = 6;
    for (int attempt = 0; attempt < kMaxAttempts; ++attempt) {
      k4siwecal::CycleAssembler assembler(maxJump, m_bcidThreshold.value());
      flushed.clear();
      if (!streamFile(path, assembler, flushed)) {
        error() << "Cannot open raw file: " << path << endmsg;
        return false;
      }
      assembler.drain(flushed);
      if (assembler.repeatedCycleCount() == 0) {
        return true;
      }
      warning() << "Repeated cycles detected in " << path << " with MaxReadoutCycleJump=" << maxJump
                << "; retrying with a larger buffer window" << endmsg;
      maxJump *= 10;
    }
    return false;
  }
};

DECLARE_COMPONENT(EcalRawDecoder)
