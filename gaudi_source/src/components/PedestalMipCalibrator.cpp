/*
 * PedestalMipCalibrator: pedestal + MIP calibration as a Gaudi component,
 * generating the same calibration tables already consumed by
 * siwecal_eventbuilder/calibration.py::from_files (3-column-per-SCA
 * pedestal format; "layer chip channel mpv empv widthmpv chi2ndf nentries"
 * MIP format), replacing the external reference tool SiWECAL-TB-analysis.
 *
 * Ported method: SiWECAL-TB-analysis/SLBperformance/DecodedSLBAnalysis.cc::
 * NSlabsAnalysis (event selection + histogram filling: badbcid==0,
 * nhits<=MaxNhit, SimpleCoincidenceTagger, 400 bins 100.5-500.5 pedestal /
 * 900 bins 100.5-1000.5 MIP histograms per (slab,chip,channel,sca)) and
 * SLBperformance/TBchecks/analysis.h::pedanalysis / mipanalysis_summary
 * (pedestal_mode=1) for the fits -- verified to be the method that actually
 * produced this repo's deployed calibration/MuonCalib_it2/{pedestals,mips}
 * files (their header lines match verbatim). Do not confuse with the
 * simpler single-Gaussian pedestal fit or the noise-covariance
 * decomposition described in earlier, less precise exploration notes --
 * those are a different reference-tool code path not used to produce this
 * repo's actual calibration data.
 *
 * Mip mode needs its own per-SCA pedestal histograms too: the reference's
 * "pedestal_mode=1" on-the-fly subtraction uses a raw truncated mean of the
 * SAME per-SCA pedestal histogram filled in the same pass, NOT the fitted
 * Pedestal output file -- so Mip mode is self-contained and does not read
 * an external pedestal file.
 *
 * Plain Gaudi::Algorithm (not a k4FWCore Producer/Transformer): this is a
 * global reduction over an entire run (histograms/fits only final after
 * every entry is seen), and its output is ASCII calibration tables, not a
 * per-event EDM4hep collection. All work happens in initialize(), same
 * self-contained-I/O pattern as EcalToEDM4hep/EcalRawDecoder.
 */
#include "k4SiWEcalReco/PedestalMipCalib.h"
#include "k4SiWEcalReco/SlbFrameDecoder.h"  // kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TTree.h"

#include <cmath>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using k4siwecal::kChannelsInSkiroc;
using k4siwecal::kScasInSkiroc;
using k4siwecal::kSkirocsPerAsu;
using k4siwecal::kSlbDepth;

namespace {
inline std::size_t histIndex(int slab, int chip, int chn, int sca) {
  return ((static_cast<std::size_t>(slab) * kSkirocsPerAsu + chip) * kChannelsInSkiroc + chn) * kScasInSkiroc + sca;
}
constexpr std::size_t kNumHist =
    static_cast<std::size_t>(kSlbDepth) * kSkirocsPerAsu * kChannelsInSkiroc * kScasInSkiroc;

// Shared binning for the pedestal-like ("ped"/"ped_high"/"ped_low") and
// MIP-like ("mip"/"mip_high"/"mip_low") histogram grids. Used both when
// building them from a siwecaldecoded tree (Mode=Pedestal|Mip|Fill) and
// when reconstructing an empty placeholder for a zero-suppressed (missing)
// key while reading a merged histogram file (Mode=FitPedestal|FitMip) --
// binning must always match exactly or a hadd merge of same-named
// histograms across Fill jobs would silently corrupt/misalign bins.
constexpr int kPedNBins = 400;
constexpr double kPedLo = 100.5;
constexpr double kPedHi = 500.5;
constexpr int kMipNBins = 900;
constexpr double kMipLo = 100.5;
constexpr double kMipHi = 1000.5;
}  // namespace

struct PedestalMipCalibrator final : Gaudi::Algorithm {
  using Gaudi::Algorithm::Algorithm;

  StatusCode initialize() override {
    TH1::AddDirectory(kFALSE);  // histograms below are owned by us, not by whatever TFile is "current"

    const bool highGain = (m_gain.value() == "high");
    if (m_mode.value() == "Fill") {
      if (m_inputFiles.value().empty()) {
        return error() << "No InputFiles provided" << endmsg, StatusCode::FAILURE;
      }
      if (m_outputHistogramFile.value().empty()) {
        return error() << "Fill mode requires OutputHistogramFile" << endmsg, StatusCode::FAILURE;
      }
      return runFill();
    }
    if (m_mode.value() == "FitPedestal") {
      if (m_inputHistogramFile.value().empty()) {
        return error() << "FitPedestal mode requires InputHistogramFile" << endmsg, StatusCode::FAILURE;
      }
      return fitPedestalFromHistograms(highGain);
    }
    if (m_mode.value() == "FitMip") {
      if (m_inputHistogramFile.value().empty()) {
        return error() << "FitMip mode requires InputHistogramFile" << endmsg, StatusCode::FAILURE;
      }
      return fitMipFromHistograms(highGain);
    }

    if (m_inputFiles.value().empty()) {
      return error() << "No InputFiles provided" << endmsg, StatusCode::FAILURE;
    }
    if (m_mode.value() == "Pedestal") {
      return runPedestal(highGain);
    }
    if (m_mode.value() == "Mip") {
      return runMip(highGain);
    }
    return error() << "Unknown Mode '" << m_mode.value() << "' (expected Pedestal|Mip|Fill|FitPedestal|FitMip)"
                    << endmsg,
           StatusCode::FAILURE;
  }

  StatusCode execute(const EventContext&) const override { return StatusCode::SUCCESS; }

  Gaudi::Property<std::vector<std::string>> m_inputFiles{
      this, "InputFiles", {}, "One or more siwecaldecoded ROOT files (multiple = accumulated statistics)"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "siwecaldecoded", "Input TTree name"};
  Gaudi::Property<std::string> m_mode{this, "Mode", "Pedestal", "Pedestal|Mip"};
  Gaudi::Property<std::string> m_gain{this, "Gain", "high", "high|low"};
  Gaudi::Property<std::string> m_outputPedestalFile{this, "OutputPedestalFile", "", "Pedestal mode output path"};
  Gaudi::Property<std::string> m_outputMipFile{this, "OutputMipFile", "", "Mip mode output path"};
  Gaudi::Property<int> m_maxNhit{this, "MaxNhit", 1, "nhits<=this to accept the SCA slice"};
  Gaudi::Property<int> m_nSlabsHit{this, "NSlabsHit", 8,
                                   "Coincidence tagger requirement: total slabs (incl. this one) expected hit"};
  Gaudi::Property<double> m_minMipIntegral{this, "MinMipIntegral", 200.,
                                           "Combined per-channel MIP histogram integral required to ATTEMPT a fit. "
                                           "A channel below this (but with >0 entries) is not masked -- it is handed "
                                           "to the chip/global fallback pass. A channel is only left masked when it "
                                           "had exactly 0 entries."};
  Gaudi::Property<double> m_chipFallbackMinIntegral{
      this, "ChipFallbackMinIntegral", 2000.,
      "For a channel with too few entries of its own: combined-chip (all 64 channels' MIP histograms summed, "
      "port of the reference tool's hmip_chip cross-check) integral required to fit a chip-level fallback MPV. "
      "Matches SiWECAL-TB-analysis/SLBperformance/TBchecks/analysis.h's own hmip_chip threshold"};
  Gaudi::Property<double> m_pedestalMaxMeanAdc{
      this, "PedestalMaxMeanAdc", 300.0,
      "Pedestal-fit outlier ceiling: a fitted SCA mean that is NaN/Inf or exceeds this is treated as unusable "
      "and filled from the channel's other usable SCAs (matches calibration/clean_pedestals.py's MAX_MEAN)"};
  Gaudi::Property<double> m_mipFallbackMaxAdc{
      this, "MipFallbackMaxAdc", 300.0,
      "When the Langau fit fails (ndf<=0), fall back to the combined histogram's peak bin as the MPV if it is "
      "positive and at or below this ceiling; otherwise the channel is masked as before. Tune down for low "
      "gain (e.g. ~30) where a real MIP peak sits at a much smaller ADC value"};
  Gaudi::Property<double> m_maxMipChi2Ndf{
      this, "MaxMipChi2Ndf", 1.1,
      "Fit-quality gate: a Langau fit whose chi2/ndf exceeds this is rejected (treated as a bad fit and sent "
      "to the chip/global fallback), even if it converged and its MPV is in range. Applies to both the "
      "per-channel fit and the chip-level fallback fit."};
  Gaudi::Property<double> m_mipHistMaxTol{
      this, "MipHistMaxTol", 0.25,
      "Histogram-max consistency gate: if the fitted MPV disagrees with the combined histogram's peak bin by "
      "more than this fraction of the MPV (|hist_max - mpv| > tol*mpv), the fit is considered bad and the "
      "channel falls back to the histogram peak bin (status 2, empv=-2) instead of using the fitted MPV."};
  Gaudi::Property<double> m_maxMipChi2NdfFallback{
      this, "MaxMipChi2NdfFallback", 25.0,
      "LOOSE chi2/ndf ceiling for the chip- and slab-level fallback fits. chi2/ndf grows with statistics for "
      "an imperfect model (Langau vs the real MIP shape): a single channel sits ~0.3, a chip (~64x stats) "
      "~2-3, a slab (~1000x) ~8-17, all with a perfectly good MPV. So the tight per-channel MaxMipChi2Ndf "
      "would reject every chip/slab fit; this much looser ceiling only catches genuinely broken fallback fits."};
  Gaudi::Property<double> m_mipLowLim{
      this, "MipLowLim", 10.0,
      "'Looks like a real MIP' range check, high gain (pllohigh in the reference tool, originally a single "
      "lower bound of 45 -- widened to a [MipLowLim, MipHighLim] window after inspecting this run's actual "
      "peak histograms, where genuine MIP peaks sit around 20-45 ADC, well below the original 45 floor; the "
      "lower bound was further lowered to 10 to accept the low-side borderline channels that were being "
      "sent to fallback despite a converged fit)"};
  Gaudi::Property<double> m_mipHighLim{
      this, "MipHighLim", 50.0,
      "Upper bound of the 'looks like a real MIP' range check, high gain -- see MipLowLim"};
  Gaudi::Property<double> m_mipLowLimLowGain{
      this, "MipLowLimLowGain", 2.0,
      "'Looks like a real MIP' lower bound, low gain (pllolow in the reference tool was 5.; scaled down "
      "proportionally to MipLowLim/45 pending real low-gain data to check against)"};
  Gaudi::Property<double> m_mipHighLimLowGain{
      this, "MipHighLimLowGain", 20.0,
      "Upper bound of the 'looks like a real MIP' range check, low gain -- see MipLowLimLowGain"};
  Gaudi::Property<std::string> m_diagnosticsFile{
      this, "DiagnosticsFile", "",
      "Optional ROOT file with cross-check summary plots (2D layer/chip vs. sca/channel maps of the written "
      "value and a per-cell status code, plus 1D value distributions); empty (default) = no plots written"};
  Gaudi::Property<std::string> m_outputHistogramFile{
      this, "OutputHistogramFile", "",
      "Mode=Fill: output ROOT file with the ped_high/ped_low/mip_high/mip_low histogram grids (zero-suppressed: "
      "only non-empty histograms are written), no fit"};
  Gaudi::Property<std::string> m_inputHistogramFile{
      this, "InputHistogramFile", "",
      "Mode=FitPedestal|FitMip: input ROOT file with a (possibly hadd-merged) histogram grid written by "
      "Mode=Fill, read instead of building histograms from InputFiles trees"};

 private:
  // Bound branch buffers for one siwecaldecoded tree, matching the layout
  // written by EcalRawDecoder / SiWECAL-TB-analysis's SLBraw2ROOT.cc.
  struct TreeBuffers {
    int nSlboards = 0;
    int bcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
    int badbcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
    int nhits[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
    int hitbitLow[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
    int hitbitHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
    int adcLow[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
    int adcHigh[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  };

  static void bindReadBranches(TTree& tree, TreeBuffers& buf) {
    tree.SetBranchAddress("n_slboards", &buf.nSlboards);
    tree.SetBranchAddress("bcid", buf.bcid);
    tree.SetBranchAddress("badbcid", buf.badbcid);
    tree.SetBranchAddress("nhits", buf.nhits);
    tree.SetBranchAddress("hitbit_low", buf.hitbitLow);
    tree.SetBranchAddress("hitbit_high", buf.hitbitHigh);
    tree.SetBranchAddress("adc_low", buf.adcLow);
    tree.SetBranchAddress("adc_high", buf.adcHigh);
  }

  // Shared NSlabsAnalysis-style event loop over all InputFiles: for every
  // (slab,chip,sca) passing the badbcid/nhits/coincidence selection, calls
  // fn(buf, slab, chip, sca) with the TreeBuffers entry already loaded (both
  // gains' hitbit/adc arrays available) -- the caller decides what to do
  // per channel/gain. This is the same selection logic
  // forEachSelectedChannel used to run inline; factored out so Mode=Fill
  // can build all 4 histogram grids (both gains) in one tree pass instead
  // of one forEachSelectedChannel call per gain.
  template <typename Fn>
  StatusCode forEachSelectedSca(const Fn& fn) const {
    for (const auto& path : m_inputFiles.value()) {
      std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
      if (!file || file->IsZombie()) {
        return error() << "Cannot open input file: " << path << endmsg, StatusCode::FAILURE;
      }
      auto* tree = dynamic_cast<TTree*>(file->Get(m_treeName.value().c_str()));
      if (!tree) {
        return error() << "Tree '" << m_treeName.value() << "' not found in " << path << endmsg,
               StatusCode::FAILURE;
      }
      // TreeBuffers is ~3.7MB (fixed [15][16][15][64] int arrays) -- too
      // large for a stack frame on the smaller worker-thread stacks
      // Gaudi/k4run can use. Heap-allocate it explicitly.
      auto buf = std::make_unique<TreeBuffers>();
      bindReadBranches(*tree, *buf);
      const Long64_t nEntries = tree->GetEntries();
      for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        for (int slab = 0; slab < buf->nSlboards; ++slab) {
          for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
            for (int sca = 0; sca < kScasInSkiroc; ++sca) {
              if (buf->badbcid[slab][chip][sca] != 0 || buf->bcid[slab][chip][sca] < 0) continue;
              if (buf->nhits[slab][chip][sca] > m_maxNhit.value()) continue;

              const int bcidSeen = k4siwecal::simpleCoincidenceTagger(buf->bcid, buf->badbcid, buf->nSlboards, slab,
                                                                       buf->bcid[slab][chip][sca]);
              if (bcidSeen < (m_nSlabsHit.value() - 1)) continue;

              fn(*buf, slab, chip, sca);
            }
          }
        }
      }
    }
    return StatusCode::SUCCESS;
  }

  // Thin wrapper over forEachSelectedSca: same selection, but iterates the
  // 64 channels of one requested gain and hands the caller (slab, chip,
  // chn, sca, hitbit, adcValue) directly, matching the original
  // (pre-refactor) callback shape used by runPedestal/runMip. Behavior is
  // byte-for-byte identical to before this was split out.
  template <typename Fn>
  StatusCode forEachSelectedChannel(bool highGain, const Fn& fn) const {
    return forEachSelectedSca([&](TreeBuffers& buf, int slab, int chip, int sca) {
      for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
        const int hitbit = highGain ? buf.hitbitHigh[slab][chip][sca][chn] : buf.hitbitLow[slab][chip][sca][chn];
        const int adcValue = highGain ? buf.adcHigh[slab][chip][sca][chn] : buf.adcLow[slab][chip][sca][chn];
        fn(slab, chip, chn, sca, hitbit, adcValue);
      }
    });
  }

  static std::vector<std::unique_ptr<TH1F>> makeHistograms(const char* prefix, int nbins, double lo, double hi) {
    std::vector<std::unique_ptr<TH1F>> hist(kNumHist);
    for (int slab = 0; slab < kSlbDepth; ++slab) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          for (int sca = 0; sca < kScasInSkiroc; ++sca) {
            hist[histIndex(slab, chip, chn, sca)] =
                std::make_unique<TH1F>(Form("%s_s%d_c%d_ch%d_sca%d", prefix, slab, chip, chn, sca), "", nbins, lo, hi);
          }
        }
      }
    }
    return hist;
  }

  StatusCode runPedestal(bool highGain) {
    auto hist = makeHistograms("ped", kPedNBins, kPedLo, kPedHi);

    const auto sc = forEachSelectedChannel(highGain, [&](int slab, int chip, int chn, int sca, int hitbit, int adc) {
      if (hitbit == 0) {
        hist[histIndex(slab, chip, chn, sca)]->Fill(adc);
      }
    });
    if (!sc.isSuccess()) return sc;

    return writePedestalTable(hist);
  }

  // Fits + writes the pedestal table from an already-populated histogram
  // grid, regardless of whether it was just built from InputFiles trees
  // (runPedestal) or read back from a hadd-merged Mode=Fill output
  // (fitPedestalFromHistograms) -- identical code path either way, so the
  // fit/robustness/diagnostics logic can never drift between the two.
  StatusCode writePedestalTable(std::vector<std::unique_ptr<TH1F>>& hist) {
    std::ofstream fout(m_outputPedestalFile.value());
    if (!fout) {
      return error() << "Cannot open output file: " << m_outputPedestalFile.value() << endmsg, StatusCode::FAILURE;
    }
    fout << "#pedestal results (fit to a gaussian) remove channels/sca with two pedestals peaks from the analysis "
            ": PROTO15\n";
    fout << "#layer chip channel";
    for (int isca = 0; isca < kScasInSkiroc; ++isca) {
      fout << " ped_mean" << isca << " ped_error" << isca << " ped_width" << isca;
    }
    fout << "\n";

    // Optional cross-check plots (see DiagnosticsFile doc): 2D maps indexed
    // exactly like the reference tool's ped_all/width_all
    // (layer*20+chip, sca*100+channel), plus a per-cell status code map
    // (0=masked, 1=genuine fit, 2=filled from the channel's other usable
    // SCAs) so the new robustness policy's decisions are directly visible,
    // and 1D distributions of the written mean/width for a quick sanity
    // check of the whole run.
    std::unique_ptr<TFile> diagFile;
    TH2F* hMeanMap = nullptr;
    TH2F* hWidthMap = nullptr;
    TH2F* hStatusMap = nullptr;
    TH1F* hMeanDist = nullptr;
    TH1F* hWidthDist = nullptr;
    if (!m_diagnosticsFile.value().empty()) {
      diagFile.reset(TFile::Open(m_diagnosticsFile.value().c_str(), "RECREATE"));
      if (!diagFile || diagFile->IsZombie()) {
        return error() << "Cannot create diagnostics file: " << m_diagnosticsFile.value() << endmsg,
               StatusCode::FAILURE;
      }
      diagFile->cd();
      hMeanMap = new TH2F("pedestal_mean_map", "Pedestal mean; layer*20+chip; sca*100+channel", 300, -0.5, 299.5,
                           1500, -0.5, 1499.5);
      hWidthMap = new TH2F("pedestal_width_map", "Pedestal width; layer*20+chip; sca*100+channel", 300, -0.5, 299.5,
                            1500, -0.5, 1499.5);
      hStatusMap = new TH2F("pedestal_status_map",
                             "Pedestal SCA status (0=masked, 1=fit, 2=filled-from-average); "
                             "layer*20+chip; sca*100+channel",
                             300, -0.5, 299.5, 1500, -0.5, 1499.5);
      hMeanDist = new TH1F("pedestal_mean_dist", "Pedestal mean distribution (unmasked SCAs); ADC", 200, 0.0, 400.0);
      hWidthDist = new TH1F("pedestal_width_dist", "Pedestal width distribution (unmasked SCAs); ADC", 100, 0.0, 20.0);
    }

    // Robustness policy, ported from calibration/clean_pedestals.py's
    // "previous correction" (already validated on real data) directly into
    // the generation step, rather than as a separate post-processing pass:
    //   - A fitted SCA is "usable" only if the fit converged (kGood) AND its
    //     mean is finite, strictly positive, and at or below
    //     PedestalMaxMeanAdc. NaN/Inf and out-of-range means (e.g. a stray
    //     high-charge tail biasing the Gaussian fit) are NOT usable.
    //   - If a channel has at least one usable SCA: usable SCAs keep their
    //     own fitted (mean, error, width); every unusable SCA is filled with
    //     the plain average of the usable SCAs' mean/width, flagged with the
    //     existing error=-10 sentinel (calibration.py::_read_pedestal_file
    //     only inspects the mean to decide masking, so this sentinel is for
    //     human/debugging visibility, not the masking decision itself).
    //   - If a channel has NO usable SCA at all: the whole channel is
    //     written as 15x "0 -10 0" -- an all-zero-means row is exactly the
    //     masking convention calibration.py::_read_pedestal_file already
    //     checks for ("canal a 0 cuenta como masked").
    const double maxMean = m_pedestalMaxMeanAdc.value();
    for (int layer = 0; layer < kSlbDepth; ++layer) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          double mean[kScasInSkiroc] = {0};
          double error[kScasInSkiroc] = {0};
          double width[kScasInSkiroc] = {0};
          bool usable[kScasInSkiroc] = {false};

          for (int sca = 0; sca < kScasInSkiroc; ++sca) {
            auto& h = *hist[histIndex(layer, chip, chn, sca)];
            const auto r = k4siwecal::fitPedestalSca(h);
            if (r.status == k4siwecal::PedestalFitResult::Status::kGood && std::isfinite(r.mean) && r.mean > 0.0 &&
                r.mean <= maxMean) {
              mean[sca] = r.mean;
              // Minos error estimation can fail (non-finite r.error) even when
              // the mean itself fitted fine; the mean/width are still usable,
              // so keep them and just zero out the unusable error rather than
              // writing a NaN that would propagate as text into the table.
              error[sca] = std::isfinite(r.error) ? r.error : 0.0;
              width[sca] = r.width;
              usable[sca] = true;
            }
          }

          int nUsable = 0;
          double sumMean = 0.0;
          double sumWidth = 0.0;
          for (int sca = 0; sca < kScasInSkiroc; ++sca) {
            if (usable[sca]) {
              ++nUsable;
              sumMean += mean[sca];
              sumWidth += width[sca];
            }
          }

          const double diagX = layer * 20 + chip;
          fout << layer << " " << chip << " " << chn << " ";
          if (nUsable == 0) {
            for (int sca = 0; sca < kScasInSkiroc; ++sca) {
              fout << "0 -10 0 ";
              if (hMeanMap != nullptr) {
                const double diagY = sca * 100 + chn;
                hMeanMap->Fill(diagX, diagY, 0.0);
                hWidthMap->Fill(diagX, diagY, 0.0);
                hStatusMap->Fill(diagX, diagY, 0.0);
              }
            }
          } else {
            const double avgMean = sumMean / nUsable;
            const double avgWidth = sumWidth / nUsable;
            for (int sca = 0; sca < kScasInSkiroc; ++sca) {
              const double m = usable[sca] ? mean[sca] : avgMean;
              const double w = usable[sca] ? width[sca] : avgWidth;
              fout << m << " " << (usable[sca] ? error[sca] : -10.0) << " " << w << " ";
              if (hMeanMap != nullptr) {
                const double diagY = sca * 100 + chn;
                hMeanMap->Fill(diagX, diagY, m);
                hWidthMap->Fill(diagX, diagY, w);
                hStatusMap->Fill(diagX, diagY, usable[sca] ? 1.0 : 2.0);
                hMeanDist->Fill(m);
                hWidthDist->Fill(w);
              }
            }
          }
          fout << "\n";
        }
      }
    }
    if (diagFile) {
      diagFile->cd();
      // TH1::AddDirectory(kFALSE) at the top of initialize() (needed so the
      // per-SCA fit histograms above aren't auto-attached to whatever TFile
      // is "current") also applies to these diagnostic histograms, so cd()
      // doesn't auto-register them -- Write() each one explicitly.
      hMeanMap->Write();
      hWidthMap->Write();
      hStatusMap->Write();
      hMeanDist->Write();
      hWidthDist->Write();
      diagFile->Close();
      info() << "PedestalMipCalibrator: wrote pedestal diagnostics to " << m_diagnosticsFile.value() << endmsg;
    }
    info() << "PedestalMipCalibrator: wrote pedestal table to " << m_outputPedestalFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode runMip(bool highGain) {
    auto pedHist = makeHistograms("mipped", kPedNBins, kPedLo, kPedHi);
    auto mipHist = makeHistograms("mip", kMipNBins, kMipLo, kMipHi);

    const auto sc = forEachSelectedChannel(highGain, [&](int slab, int chip, int chn, int sca, int hitbit, int adc) {
      const auto idx = histIndex(slab, chip, chn, sca);
      if (hitbit == 0) {
        pedHist[idx]->Fill(adc);
      } else if (hitbit == 1) {
        mipHist[idx]->Fill(adc);
      }
    });
    if (!sc.isSuccess()) return sc;

    return writeMipTable(pedHist, mipHist, highGain);
  }

  // One (layer,chip,channel) MIP result, before the chip/global fallback
  // pass below is applied. status: 0=masked/lowstats/bad-fit (mpv==0, needs
  // fallback), 1=genuine fit, 2=histogram-peak fallback. (status 3 -- "fit
  // succeeded but the MPV is outside the plausible MIP range" -- is no longer
  // produced: such a fit is now considered bad and routed to the chip/global
  // fallback like any other low-stats channel, rather than kept as-is.)
  struct MipEntry {
    double mpv = 0., empv = -10., width = 0., chi2ndf = 0., integral = 0.;
    int status = 0;
  };

  // Fits + writes the MIP table from already-populated pedestal-like and
  // MIP-like histogram grids -- same rationale as writePedestalTable:
  // identical code path whether the grids were just built (runMip) or read
  // back from a merged Mode=Fill file (fitMipFromHistograms).
  StatusCode writeMipTable(std::vector<std::unique_ptr<TH1F>>& pedHist, std::vector<std::unique_ptr<TH1F>>& mipHist,
                            bool highGain) {
    std::ofstream fout(m_outputMipFile.value());
    if (!fout) {
      return error() << "Cannot open output file: " << m_outputMipFile.value() << endmsg, StatusCode::FAILURE;
    }
    fout << "#mip results PROTO15-TB2022-03\n";
    fout << "#layer chip channel mpv empv widthmpv chi2ndf nentries\n";
    // empv sentinel scheme: -10 = masked (mpv=0) -- the channel had exactly 0
    // entries (it never fired this run) or no fallback was available anywhere;
    // -5 = Langau fit failed (ndf<=0) and the histogram-peak fallback below
    // was also unreasonable (masked); -6 = filled from a Langau fit to the
    // whole SLAB's combined histogram (all 16 chips, tried after the chip fit
    // and before the global average); -4 = no fit possible for this channel,
    // its chip, or its slab, filled from the GLOBAL average of every
    // genuinely-fit channel in the whole calibration; -3 = no usable
    // statistics for this channel, filled from a Langau fit to the WHOLE
    // CHIP's combined MIP histogram (all 64 channels x 15 SCA -- MIP response
    // is expected to be reasonably uniform across channels of the same
    // ASIC/chip, so this is a much better substitute than leaving the channel
    // at 0); -2 = fit failed but the histogram-peak fallback was accepted (mpv
    // is a rough peak estimate, not a real fit -- width/chi2ndf are left at 0);
    // >=0 = a genuine fit; for these the written mpv is the Langau fit CURVE's
    // MPV parameter (not the histogram peak bin). A fit that converges but is
    // bad -- MPV outside the [MipLowLim, MipHighLim] window, or chi2/ndf above
    // MaxMipChi2Ndf -- is NOT written as its own value: it is treated as a bad
    // fit and routed to the -3/-4 fallback instead.
    //
    // IMPORTANT: the -3/-4 chip/global fallbacks are applied to any channel
    // that has SOME raw statistics (Integral > 0) but either not enough to fit
    // (Integral <= MinMipIntegral), or whose own fit was bad (out of range or
    // chi2/ndf too high). Only a channel with Integral==0 (it never fired at all in this
    // run) is left masked (mpv=0, empv=-10) permanently -- it never reaches
    // pass 2/3, however dry its chip or the global average is. There is
    // nothing to calibrate
    // for a channel that produced zero MIP-tagged hits, and inventing a
    // value for it would silently misrepresent a dead/unread channel as a
    // working one downstream (CalibrationTables::readPedestalFile-style
    // consumers only see mpv, not the raw entry count).

    // Range check, not just a floor (see MipLowLim/MipHighLim doc): a real
    // MIP peak for this run/threshold sits well below the reference tool's
    // original 45 ADC floor, so a bare lower bound alone would accept
    // nothing; an upper bound also guards against a peak bin picked from
    // noise/tail far out in the histogram.
    const double lowlim = highGain ? m_mipLowLim.value() : m_mipLowLimLowGain.value();
    const double highlim = highGain ? m_mipHighLim.value() : m_mipHighLimLowGain.value();

    // Pass 1: compute the raw per-channel result (fit, histogram-peak
    // fallback, or masked) exactly as before, but store it instead of
    // writing it out immediately -- the chip/global fallback pass below
    // needs every channel's result available before it can average them.
    std::vector<MipEntry> entries(static_cast<size_t>(kSlbDepth) * kSkirocsPerAsu * kChannelsInSkiroc);
    auto chanIndex = [](int layer, int chip, int chn) {
      return (static_cast<size_t>(layer) * kSkirocsPerAsu + chip) * kChannelsInSkiroc + chn;
    };

    for (int layer = 0; layer < kSlbDepth; ++layer) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          // Combine the 15 per-SCA MIP histograms into one, subtracting
          // each SCA's own on-the-fly (truncated-mean) pedestal -- port of
          // mipanalysis_summary's pedestal_mode=1 combination loop. Bin
          // range/loop bounds (0..899, i.e. the underflow bin plus all but
          // the last valid bin) match the reference exactly, quirks
          // included.
          TH1F hmips("hmips_tmp", "", 900, -100.5, 799.5);
          for (int sca = 0; sca < kScasInSkiroc; ++sca) {
            const auto idx = histIndex(layer, chip, chn, sca);
            TH1F& mh = *mipHist[idx];
            TH1F& ph = *pedHist[idx];
            double pedMean = 0.;
            if (ph.GetEntries() > 0) {
              ph.GetXaxis()->SetRangeUser(ph.GetMean() - 20, ph.GetMean() + 20);
              pedMean = ph.GetMean();
            }
            if (pedMean <= 0) continue;
            for (int k = 0; k < 900; ++k) {
              const double y = mh.GetBinContent(k);
              if (y > 0) {
                hmips.Fill(static_cast<int>(mh.GetXaxis()->GetBinCenter(k) - pedMean), y);
              }
            }
          }

          MipEntry entry;
          entry.integral = hmips.Integral();
          if (entry.integral > m_minMipIntegral.value()) {
            const auto r = k4siwecal::fitLangau(hmips, highGain);
            if (r.ndf > 0) {
              // Report the Langau fit's own MPV (the fitted curve's most-
              // probable value, r.mpv = fit parameter "MP") for a genuine fit.
              // The fit must actually be GOOD to be accepted: enough entries in
              // the peak region, MPV inside the plausible MIP window, AND
              // chi2/ndf under the ceiling. Any of these failing -> treat as a
              // BAD fit and hand the channel to the chip/global fallback pass
              // (`entry` is left at its default status 0, integral set, so
              // pass 2/3 fill it; only a 0-entry channel is ever masked).
              const double histMax = hmips.GetXaxis()->GetBinCenter(hmips.GetMaximumBin());
              double mip1 = 0.;
              for (int k = 0; k < 900; ++k) {
                if (hmips.GetBinCenter(k) < r.mpv * 4) {
                  mip1 += hmips.GetBinContent(k);
                }
              }
              if (mip1 > 4. && r.mpv > lowlim && r.mpv < highlim && r.chi2ndf <= m_maxMipChi2Ndf.value()) {
                if (std::fabs(histMax - r.mpv) <= m_mipHistMaxTol.value() * r.mpv) {
                  // Fit is good AND consistent with the visible peak: keep it.
                  entry = {r.mpv, r.empv, r.width, r.chi2ndf, entry.integral, 1};
                } else if (histMax > 0.0 && histMax <= m_mipFallbackMaxAdc.value()) {
                  // Fit passed the quality gates but its MPV disagrees with the
                  // histogram's peak bin by more than MipHistMaxTol -- distrust
                  // the fit and fall back to the histogram peak (status 2).
                  entry = {histMax, -2., 0., 0., entry.integral, 2};
                }
                // else: fit inconsistent AND histMax itself unreasonable ->
                // leave status 0 for the chip/global fallback.
              }
            } else {
              // Fit failed (ndf<=0): instead of masking immediately, check
              // the combined histogram's peak bin as a rough MPV estimate
              // (the same idea as the pre-Gaudi Python
              // Calibration._compute_mips histogram-peak method). Accepted
              // only if positive and within MipFallbackMaxAdc; otherwise
              // leave it masked for the chip/global fallback pass below.
              const int peakBin = hmips.GetMaximumBin();
              const double peakValue = hmips.GetXaxis()->GetBinCenter(peakBin);
              if (peakValue > 0.0 && peakValue <= m_mipFallbackMaxAdc.value()) {
                entry = {peakValue, -2., 0., 0., entry.integral, 2};
              } else {
                entry.status = 0;  // masked/lowstats -- resolved in pass 2/3
              }
            }
          }
          entries[chanIndex(layer, chip, chn)] = entry;
        }
      }
    }

    // Pass 2: chip-level fallback. A masked channel (status==0) borrows the
    // MPV of a Langau fit to the WHOLE CHIP's combined MIP histogram (all
    // 64 channels x 15 SCA summed into one, each SCA still pedestal-
    // subtracted the same way as the per-channel hmips in pass 1) -- port
    // of the reference tool's hmip_chip cross-check
    // (SiWECAL-TB-analysis/SLBperformance/TBchecks/analysis.h, ~line 368-
    // 480), which fits that same combined histogram once it clears its own
    // (much higher) integral threshold. A real fit on ~64x the statistics
    // of a single channel is a far more robust fallback than averaging the
    // handful (often zero or one) of already-fit donor channels in the
    // chip, which is all an earlier version of this pass did.
    int nChipFallback = 0;
    for (int layer = 0; layer < kSlbDepth; ++layer) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        bool needsFallback = false;
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          if (entries[chanIndex(layer, chip, chn)].status == 0) {
            needsFallback = true;
            break;
          }
        }
        if (!needsFallback) continue;  // every channel in this chip already has its own value

        TH1F hmipsChip("hmips_chip_tmp", "", 900, -100.5, 799.5);
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          for (int sca = 0; sca < kScasInSkiroc; ++sca) {
            const auto idx = histIndex(layer, chip, chn, sca);
            TH1F& mh = *mipHist[idx];
            TH1F& ph = *pedHist[idx];
            double pedMean = 0.;
            if (ph.GetEntries() > 0) {
              ph.GetXaxis()->SetRangeUser(ph.GetMean() - 20, ph.GetMean() + 20);
              pedMean = ph.GetMean();
            }
            if (pedMean <= 0) continue;
            for (int k = 0; k < 900; ++k) {
              const double y = mh.GetBinContent(k);
              if (y > 0) {
                hmipsChip.Fill(static_cast<int>(mh.GetXaxis()->GetBinCenter(k) - pedMean), y);
              }
            }
          }
        }

        if (hmipsChip.Integral() <= m_chipFallbackMinIntegral.value()) continue;  // still dry -- pass 2.5/3 handle it
        const auto rChip = k4siwecal::fitLangau(hmipsChip, highGain);
        // The chip-level fit must converge, have its MPV in range, and pass the
        // LOOSE fallback chi2 ceiling (a chip histogram's chi2/ndf is naturally
        // ~2-3 for a good fit -- see MaxMipChi2NdfFallback). If not, pass 2.5
        // (slab) / pass 3 (global) handle these channels. Report the fit MPV.
        if (rChip.ndf <= 0 || rChip.chi2ndf > m_maxMipChi2NdfFallback.value() ||
            rChip.mpv <= lowlim || rChip.mpv >= highlim) {
          continue;
        }
        const double peakMpvChip = rChip.mpv;
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          MipEntry& e = entries[chanIndex(layer, chip, chn)];
          if (e.status != 0) continue;
          // A channel with literally zero raw MIP entries (as opposed to
          // "some entries but below MinMipIntegral") never triggered at all
          // in this run -- there is nothing to calibrate, and no amount of
          // fallback should invent a value for it. Leave it masked (mpv=0)
          // permanently; do not count it as a fallback.
          if (e.integral == 0.) continue;
          e.mpv = peakMpvChip;
          e.empv = -3.;
          e.status = 4;
          ++nChipFallback;
        }
      }
    }

    // Pass 2.5: slab (layer) level fallback. For any channel the chip fallback
    // could not fill (chip histogram too dry, or its fit rejected), fit the
    // WHOLE SLAB's combined MIP histogram (all 16 chips x 64 channels x 15 SCA,
    // each SCA pedestal-subtracted the same way) BEFORE giving up to the flat
    // detector-wide global average -- a per-slab MPV is far more local/correct
    // for a low-MPV slab than the global mean. Uses the loose fallback chi2
    // ceiling (a slab has ~1000x a channel's statistics, chi2/ndf ~8-17 even
    // for a good fit).
    int nSlabFallback = 0;
    for (int layer = 0; layer < kSlbDepth; ++layer) {
      bool needsFallback = false;
      for (int chip = 0; chip < kSkirocsPerAsu && !needsFallback; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          const MipEntry& e = entries[chanIndex(layer, chip, chn)];
          if (e.status == 0 && e.integral > 0.) {
            needsFallback = true;
            break;
          }
        }
      }
      if (!needsFallback) continue;

      TH1F hmipsSlab("hmips_slab_tmp", "", 900, -100.5, 799.5);
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          for (int sca = 0; sca < kScasInSkiroc; ++sca) {
            const auto idx = histIndex(layer, chip, chn, sca);
            TH1F& mh = *mipHist[idx];
            TH1F& ph = *pedHist[idx];
            double pedMean = 0.;
            if (ph.GetEntries() > 0) {
              ph.GetXaxis()->SetRangeUser(ph.GetMean() - 20, ph.GetMean() + 20);
              pedMean = ph.GetMean();
            }
            if (pedMean <= 0) continue;
            for (int k = 0; k < 900; ++k) {
              const double y = mh.GetBinContent(k);
              if (y > 0) {
                hmipsSlab.Fill(static_cast<int>(mh.GetXaxis()->GetBinCenter(k) - pedMean), y);
              }
            }
          }
        }
      }

      if (hmipsSlab.Integral() <= m_chipFallbackMinIntegral.value()) continue;
      const auto rSlab = k4siwecal::fitLangau(hmipsSlab, highGain);
      if (rSlab.ndf <= 0 || rSlab.chi2ndf > m_maxMipChi2NdfFallback.value() ||
          rSlab.mpv <= lowlim || rSlab.mpv >= highlim) {
        continue;  // slab fit unusable -- pass 3 (global) handles these
      }
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          MipEntry& e = entries[chanIndex(layer, chip, chn)];
          if (e.status != 0 || e.integral == 0.) continue;
          e.mpv = rSlab.mpv;
          e.empv = -6.;
          e.status = 6;
          ++nSlabFallback;
        }
      }
    }

    // No global fallback. A channel that neither it, nor its chip, nor its slab
    // could fit is MASKED -- it does not inherit a detector-wide average MPV.
    //
    // That average was a fiction. It pooled every fitted channel in the detector,
    // including slab 12, the chip-on-board: by design its gain is ~1.65x every
    // other slab's (MPV ~35 vs ~21.5), so the "average MIP" it produced belongs
    // to no real channel anywhere. The one channel that ever reached it in th220
    // (slab 14, chip 0, channel 27) was handed 22.62 while its own slab's mean is
    // 25.11 and its immediate neighbours fit at 23-26 -- a ~10% error on every hit
    // that channel would ever record. A wrong MIP is worse than no MIP: a masked
    // channel is excluded downstream and its absence is visible, whereas a
    // fabricated one silently miscalibrates every event it touches.
    //
    // If a channel fired but nothing at any granularity could be fitted from it,
    // that is a fact worth surfacing, not papering over.
    int nStillMasked = 0;
    int nZeroEntryMasked = 0;
    for (const auto& e : entries) {
      if (e.status != 0) continue;
      if (e.integral == 0.) {
        ++nZeroEntryMasked;  // never fired in this run
      } else {
        ++nStillMasked;      // fired, but unfittable at channel, chip AND slab level
      }
    }

    // Optional cross-check plots (see DiagnosticsFile doc): 2D maps indexed
    // like the reference tool's mpv_all (layer*20+chip, channel), plus a
    // per-channel status code map (0=masked/nothing available anywhere,
    // 1=genuine fit, 2=histogram-peak fallback, 3=fit-not-miplike,
    // 4=chip-fit fallback, 6=slab-fit fallback; there is no global fallback) so the fallback
    // policy's decisions are directly visible, and 1D distributions of the
    // written MPV/chi2ndf.
    std::unique_ptr<TFile> diagFile;
    TH2F* hMpvMap = nullptr;
    TH2F* hStatusMap = nullptr;
    TH1F* hMpvDist = nullptr;
    TH1F* hChi2NdfDist = nullptr;
    if (!m_diagnosticsFile.value().empty()) {
      diagFile.reset(TFile::Open(m_diagnosticsFile.value().c_str(), "RECREATE"));
      if (!diagFile || diagFile->IsZombie()) {
        return error() << "Cannot create diagnostics file: " << m_diagnosticsFile.value() << endmsg,
               StatusCode::FAILURE;
      }
      diagFile->cd();
      hMpvMap = new TH2F("mip_mpv_map", "MIP MPV; layer*20+chip; channel", 300, -0.5, 299.5, 64, -0.5, 63.5);
      hStatusMap = new TH2F("mip_status_map",
                             "MIP status (0=masked, 1=fit, 2=histfallback, "
                             "4=chip-fallback, 5=global-fallback, 6=slab-fallback); layer*20+chip; channel",
                             300, -0.5, 299.5, 64, -0.5, 63.5);
      hMpvDist = new TH1F("mip_mpv_dist", "MIP MPV distribution (fit or fallback); ADC", 200, 0.0, 400.0);
      hChi2NdfDist = new TH1F("mip_chi2ndf_dist", "MIP fit chi2/ndf distribution (genuine fits only)", 100, 0.0, 20.0);
    }

    for (int layer = 0; layer < kSlbDepth; ++layer) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          const MipEntry& e = entries[chanIndex(layer, chip, chn)];
          fout << layer << " " << chip << " " << chn << " " << e.mpv << " " << e.empv << " " << e.width << " "
               << e.chi2ndf << " " << e.integral << "\n";
          if (hMpvMap != nullptr) {
            const double diagX = layer * 20 + chip;
            const double diagY = chn;
            hMpvMap->Fill(diagX, diagY, e.mpv);
            hStatusMap->Fill(diagX, diagY, static_cast<double>(e.status));
            if (e.status != 0) hMpvDist->Fill(e.mpv);
            if (e.status == 1) hChi2NdfDist->Fill(e.chi2ndf);
          }
        }
      }
    }
    if (diagFile) {
      diagFile->cd();
      // See the matching comment in runPedestal(): AddDirectory(kFALSE)
      // means cd() alone doesn't register these -- Write() explicitly.
      hMpvMap->Write();
      hStatusMap->Write();
      hMpvDist->Write();
      hChi2NdfDist->Write();
      diagFile->Close();
      info() << "PedestalMipCalibrator: wrote MIP diagnostics to " << m_diagnosticsFile.value() << endmsg;
    }
    info() << "PedestalMipCalibrator: MIP fallback summary: " << nChipFallback << " channel(s) filled from their "
           << "chip fit, " << nSlabFallback << " channel(s) filled from their slab fit, " << nZeroEntryMasked
           << " channel(s) masked (zero raw MIP entries -- never fired in this run)" << endmsg;
    if (nStillMasked > 0) {
      // Not a detail: these channels DID fire, and nothing at any granularity
      // could be fitted from them. They used to be quietly handed the
      // detector-wide average MPV. Now they are masked, and said out loud.
      warning() << "PedestalMipCalibrator: " << nStillMasked
                << " channel(s) had MIP statistics but could not be fitted at channel, chip OR slab level -- "
                << "MASKED (no MIP). Investigate before trusting this calibration." << endmsg;
    }
    info() << "PedestalMipCalibrator: wrote MIP table to " << m_outputMipFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  // Mode=Fill: one forEachSelectedSca pass builds all 4 histogram grids
  // (both gains x ped/mip) at once -- avoids re-reading the same input
  // tree once per gain/mode the way the Condor pipeline's Fill stage would
  // otherwise need to. "ped_high"/"ped_low" here are the same thing
  // runPedestal calls "ped" and runMip calls "mipped" (identical binning,
  // identical hitbit==0 fill condition) -- built once, reused by both
  // fitPedestalFromHistograms and fitMipFromHistograms's on-the-fly
  // pedestal subtraction.
  StatusCode runFill() {
    auto pedHigh = makeHistograms("ped_high", kPedNBins, kPedLo, kPedHi);
    auto pedLow = makeHistograms("ped_low", kPedNBins, kPedLo, kPedHi);
    auto mipHigh = makeHistograms("mip_high", kMipNBins, kMipLo, kMipHi);
    auto mipLow = makeHistograms("mip_low", kMipNBins, kMipLo, kMipHi);

    const auto sc = forEachSelectedSca([&](TreeBuffers& buf, int slab, int chip, int sca) {
      for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
        const auto idx = histIndex(slab, chip, chn, sca);
        const int hitbitHigh = buf.hitbitHigh[slab][chip][sca][chn];
        const int hitbitLow = buf.hitbitLow[slab][chip][sca][chn];
        if (hitbitHigh == 0) {
          pedHigh[idx]->Fill(buf.adcHigh[slab][chip][sca][chn]);
        } else if (hitbitHigh == 1) {
          mipHigh[idx]->Fill(buf.adcHigh[slab][chip][sca][chn]);
        }
        if (hitbitLow == 0) {
          pedLow[idx]->Fill(buf.adcLow[slab][chip][sca][chn]);
        } else if (hitbitLow == 1) {
          mipLow[idx]->Fill(buf.adcLow[slab][chip][sca][chn]);
        }
      }
    });
    if (!sc.isSuccess()) return sc;

    std::unique_ptr<TFile> fout(TFile::Open(m_outputHistogramFile.value().c_str(), "RECREATE"));
    if (!fout || fout->IsZombie()) {
      return error() << "Cannot create output histogram file: " << m_outputHistogramFile.value() << endmsg,
             StatusCode::FAILURE;
    }
    fout->cd();
    // Zero-suppression is mandatory, not an optimization: at kNumHist=230400
    // cells per grid x 4 grids, an unsuppressed file is ~300MB even
    // completely empty (TH1F/TAxis object overhead dominates over bin
    // content at this object count) -- a single large run's worth of
    // per-chunk Fill files would otherwise blow past the available EOS
    // scratch. A missing key is treated as an empty histogram on the read
    // side (readHistogramGrid), so this is safe.
    std::size_t written = 0;
    std::size_t total = 0;
    for (auto* grid : {&pedHigh, &pedLow, &mipHigh, &mipLow}) {
      for (auto& h : *grid) {
        ++total;
        if (h->GetEntries() > 0) {
          h->Write();
          ++written;
        }
      }
    }
    fout->Close();
    info() << "PedestalMipCalibrator: wrote " << written << "/" << total << " non-empty histograms to "
           << m_outputHistogramFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  // Reconstructs one histogram grid from a (possibly hadd-merged) Mode=Fill
  // output file. A missing key means every contributing Fill job saw zero
  // entries for that cell (zero-suppressed on write) -- reconstructed as an
  // empty histogram, not an error, so the fit code downstream sees exactly
  // what it would have if the grid had never been zero-suppressed.
  std::vector<std::unique_ptr<TH1F>> readHistogramGrid(TFile& file, const char* prefix, int nbins, double lo,
                                                        double hi) const {
    std::vector<std::unique_ptr<TH1F>> hist(kNumHist);
    for (int slab = 0; slab < kSlbDepth; ++slab) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        for (int chn = 0; chn < kChannelsInSkiroc; ++chn) {
          for (int sca = 0; sca < kScasInSkiroc; ++sca) {
            const auto idx = histIndex(slab, chip, chn, sca);
            auto* h = file.Get<TH1F>(Form("%s_s%d_c%d_ch%d_sca%d", prefix, slab, chip, chn, sca));
            if (h != nullptr) {
              // AddDirectory(kFALSE) (top of initialize()) means this isn't
              // owned by `file` -- safe to keep past file.Close().
              hist[idx] = std::unique_ptr<TH1F>(h);
            } else {
              hist[idx] = std::make_unique<TH1F>("", "", nbins, lo, hi);
            }
          }
        }
      }
    }
    return hist;
  }

  // Mode=FitPedestal: reads InputHistogramFile (the hadd-merged output of
  // one or more Mode=Fill runs) instead of building histograms from
  // InputFiles trees, then reuses writePedestalTable verbatim -- same fit,
  // robustness policy and diagnostics as the sequential Mode=Pedestal path.
  StatusCode fitPedestalFromHistograms(bool highGain) {
    std::unique_ptr<TFile> fin(TFile::Open(m_inputHistogramFile.value().c_str(), "READ"));
    if (!fin || fin->IsZombie()) {
      return error() << "Cannot open input histogram file: " << m_inputHistogramFile.value() << endmsg,
             StatusCode::FAILURE;
    }
    auto hist = readHistogramGrid(*fin, highGain ? "ped_high" : "ped_low", kPedNBins, kPedLo, kPedHi);
    fin->Close();
    return writePedestalTable(hist);
  }

  // Mode=FitMip: same idea as fitPedestalFromHistograms, but also reads the
  // ped_high/ped_low grid (Mode=Fill's replacement for runMip's "mipped"
  // grid -- identical binning/fill condition, see runFill) for the
  // on-the-fly pedestal subtraction writeMipTable already does.
  StatusCode fitMipFromHistograms(bool highGain) {
    std::unique_ptr<TFile> fin(TFile::Open(m_inputHistogramFile.value().c_str(), "READ"));
    if (!fin || fin->IsZombie()) {
      return error() << "Cannot open input histogram file: " << m_inputHistogramFile.value() << endmsg,
             StatusCode::FAILURE;
    }
    auto pedHist = readHistogramGrid(*fin, highGain ? "ped_high" : "ped_low", kPedNBins, kPedLo, kPedHi);
    auto mipHist = readHistogramGrid(*fin, highGain ? "mip_high" : "mip_low", kMipNBins, kMipLo, kMipHi);
    fin->Close();
    return writeMipTable(pedHist, mipHist, highGain);
  }
};

DECLARE_COMPONENT(PedestalMipCalibrator)
