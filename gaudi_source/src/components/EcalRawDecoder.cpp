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
#include "k4SiWEcalReco/AcquisitionTreeIO.h"
#include "k4SiWEcalReco/RunSettings.h"
#include "k4SiWEcalReco/SlbFrameDecoder.h"

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include "TFile.h"
#include "TString.h"
#include "TTree.h"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>

struct EcalRawDecoder final : Gaudi::Algorithm {
  using Gaudi::Algorithm::Algorithm;

  StatusCode initialize() override {
    if (m_inputFiles.value().empty()) {
      return error() << "No InputFiles provided" << endmsg, StatusCode::FAILURE;
    }

    // ZSTD-5 instead of ROOT's ZLIB-1 default: 18.6x smaller output on these
    // mostly-empty fixed-size arrays -- see k4siwecal::kDecodedFileCompression.
    auto fout = std::unique_ptr<TFile>(
        TFile::Open(m_outputFile.value().c_str(), "RECREATE", "", k4siwecal::kDecodedFileCompression));
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
    k4siwecal::bindAcquisitionBranches(*tree, *acq, runSettings);

    // Stream each flushed acquisition straight into the tree and drop it,
    // instead of collecting a whole file's worth in a vector first: an
    // Acquisition is ~5.5MB and a sparse (muon) chunk decodes into thousands
    // of them, so buffering them all peaked at >20GB per job. The sink keeps
    // memory flat (only the cycle-buffer window is ever live).
    long long totalCycles = 0;
    const auto sink = [&](const k4siwecal::Acquisition& a) {
      *acq = a;
      tree->Fill();
      ++totalCycles;
    };

    if (m_resetStatePerInputFile.value()) {
      for (const auto& path : m_inputFiles.value()) {
        const long long before = totalCycles;
        if (!decodeFileWithRetry(path, sink)) {
          return error() << "Could not decode " << path << " without repeated cycles after retries" << endmsg,
                 StatusCode::FAILURE;
        }
        info() << "EcalRawDecoder: decoded " << (totalCycles - before) << " acquisition cycles from " << path
               << endmsg;
      }
    } else {
      k4siwecal::CycleAssembler assembler(m_maxReadoutCycleJump.value(), m_bcidThreshold.value());
      for (const auto& path : m_inputFiles.value()) {
        if (!streamFile(path, assembler, sink)) {
          return error() << "Cannot open raw file: " << path << endmsg, StatusCode::FAILURE;
        }
      }
      drainAssembler(assembler, sink);
      if (assembler.repeatedCycleCount() > 0) {
        warning() << "Repeated cycles detected across the combined input; MaxReadoutCycleJump="
                  << m_maxReadoutCycleJump.value()
                  << " may be too small (ResetStatePerInputFile=false disables the automatic retry)" << endmsg;
      }
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

  // Streams one raw file's frames into `assembler`, invoking `sink` for each
  // acquisition the moment it flushes (the transient `flushed` vector holds at
  // most the one cycle a single addFrame can release, so memory stays flat).
  // Returns false only if the file cannot be opened.
  template <typename Sink>
  bool streamFile(const std::string& path, k4siwecal::CycleAssembler& assembler, const Sink& sink) const {
    std::ifstream fin(path, std::ios::binary);
    if (!fin.is_open()) return false;
    k4siwecal::RawFrameReader reader(fin, m_eudaqFormat.value());
    std::vector<unsigned char> frameBytes;
    std::vector<k4siwecal::Acquisition> flushed;
    while (reader.nextFrame(frameBytes)) {
      assembler.addFrame(frameBytes, flushed);
      for (const auto& a : flushed) sink(a);
      flushed.clear();
    }
    return true;
  }

  // Flushes the assembler's remaining buffered cycles at EOF through `sink`.
  // The transient vector is bounded by the cycle-buffer window (~MaxReadoutCycleJump
  // acquisitions), not the whole file.
  template <typename Sink>
  static void drainAssembler(k4siwecal::CycleAssembler& assembler, const Sink& sink) {
    std::vector<k4siwecal::Acquisition> flushed;
    assembler.drain(flushed);
    for (const auto& a : flushed) sink(a);
  }

  // Cheap key-only pre-scan: replays CycleAssembler's exact eviction rule
  // (evict the smallest-key cycle once more than `maxJump` distinct cycles are
  // buffered) reading ONLY each frame's cycle-id key (decodeCycleId, ~32
  // bytes/frame) -- never decoding the frame payload's Gray-coded ADC/BCID
  // data. Returns CycleAssembler::repeatedCycleCount() for this window without
  // paying the decode cost, so decodeFileWithRetry can pick the right window
  // BEFORE the single streaming decode pass fills the tree (a streaming pass
  // can't be un-filled, so the window must be settled up front). Open failures
  // return 0 and are surfaced later by the decode pass itself.
  std::size_t countRepeatedCyclesKeyOnly(const std::string& path, int maxJump) const {
    std::ifstream fin(path, std::ios::binary);
    if (!fin.is_open()) return 0;
    k4siwecal::RawFrameReader reader(fin, m_eudaqFormat.value());
    std::vector<unsigned char> frameBytes;
    std::set<int> buffered;      // cycle keys currently in the window (mirrors m_cycles)
    std::vector<int> seenKeys;   // every key as it first entered the window (mirrors m_seenKeys)
    while (reader.nextFrame(frameBytes)) {
      const int key = k4siwecal::decodeCycleId(frameBytes);
      if (buffered.insert(key).second) {
        seenKeys.push_back(key);
        if (static_cast<int>(buffered.size()) > maxJump) {
          buffered.erase(buffered.begin());  // smallest key: std::set is sorted, like m_cycles.begin()
        }
      }
    }
    std::vector<int> uniq = seenKeys;
    std::sort(uniq.begin(), uniq.end());
    uniq.erase(std::unique(uniq.begin(), uniq.end()), uniq.end());
    return seenKeys.size() - uniq.size();
  }

  // Decodes one chunk file, escalating MaxReadoutCycleJump by x10 (bounded to a
  // handful of attempts) if repeated cycles are detected, mirroring
  // ConvertDirectorySL_Raw.cc:66-71's `while(result==false)` retry loop. The
  // window is chosen by the key-only pre-scan first, then a single streaming
  // decode pass emits acquisitions through `sink` -- so a repeated-cycle window
  // never reaches the tree.
  template <typename Sink>
  bool decodeFileWithRetry(const std::string& path, const Sink& sink) const {
    int maxJump = m_maxReadoutCycleJump.value();
    constexpr int kMaxAttempts = 6;
    for (int attempt = 0; attempt < kMaxAttempts; ++attempt) {
      if (countRepeatedCyclesKeyOnly(path, maxJump) == 0) {
        k4siwecal::CycleAssembler assembler(maxJump, m_bcidThreshold.value());
        if (!streamFile(path, assembler, sink)) {
          error() << "Cannot open raw file: " << path << endmsg;
          return false;
        }
        drainAssembler(assembler, sink);
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
