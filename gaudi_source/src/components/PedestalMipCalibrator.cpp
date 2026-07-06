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
  Gaudi::Property<double> m_minMipIntegral{this, "MinMipIntegral", 500.,
                                           "Combined per-channel MIP histogram integral required to attempt a fit"};
  Gaudi::Property<double> m_pedestalMaxMeanAdc{
      this, "PedestalMaxMeanAdc", 300.0,
      "Pedestal-fit outlier ceiling: a fitted SCA mean that is NaN/Inf or exceeds this is treated as unusable "
      "and filled from the channel's other usable SCAs (matches calibration/clean_pedestals.py's MAX_MEAN)"};
  Gaudi::Property<double> m_mipFallbackMaxAdc{
      this, "MipFallbackMaxAdc", 300.0,
      "When the Langau fit fails (ndf<=0), fall back to the combined histogram's peak bin as the MPV if it is "
      "positive and at or below this ceiling; otherwise the channel is masked as before. Tune down for low "
      "gain (e.g. ~30) where a real MIP peak sits at a much smaller ADC value"};
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
              error[sca] = r.error;
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
    // empv sentinel scheme: -10 = too few entries (Integral<=MinMipIntegral,
    // masked); -5 = Langau fit failed (ndf<=0) and the histogram-peak
    // fallback below was also unreasonable (masked); -2 = fit failed but the
    // histogram-peak fallback was accepted (mpv is a rough peak estimate,
    // not a real fit -- width/chi2ndf are left at 0); -1 = fit succeeded but
    // didn't pass the "looks like a real MIP" check; >=0 = a genuine fit.

    const double lowlim = highGain ? 45. : 5.;  // pllohigh / pllolow (TBchecks/analysis.h)

    // Optional cross-check plots (see DiagnosticsFile doc): 2D maps indexed
    // like the reference tool's mpv_all (layer*20+chip, channel), plus a
    // per-channel status code map (0=masked/low-stats, 1=genuine fit,
    // 2=histogram-peak fallback used, 3=fit succeeded but doesn't look like
    // a real MIP) so the new robustness policy's decisions are directly
    // visible, and 1D distributions of the written MPV/chi2ndf.
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
                             "MIP status (0=masked/lowstats, 1=fit, 2=histfallback, 3=fit-not-miplike); "
                             "layer*20+chip; channel",
                             300, -0.5, 299.5, 64, -0.5, 63.5);
      hMpvDist = new TH1F("mip_mpv_dist", "MIP MPV distribution (fit or fallback); ADC", 200, 0.0, 400.0);
      hChi2NdfDist = new TH1F("mip_chi2ndf_dist", "MIP fit chi2/ndf distribution (genuine fits only)", 100, 0.0, 20.0);
    }

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

          const double diagX = layer * 20 + chip;
          const double diagY = chn;
          if (hmips.Integral() > m_minMipIntegral.value()) {
            const auto r = k4siwecal::fitLangau(hmips, highGain);
            if (r.ndf > 0) {
              double mip1 = 0.;
              for (int k = 0; k < 900; ++k) {
                const double center = hmips.GetBinCenter(k);
                if (center < r.mpv * 4) {
                  mip1 += hmips.GetBinContent(k);
                }
              }
              if (mip1 > 4. && r.mpv > lowlim) {
                fout << layer << " " << chip << " " << chn << " " << r.mpv << " " << r.empv << " " << r.width << " "
                     << r.chi2ndf << " " << hmips.Integral() << "\n";
                if (hMpvMap != nullptr) {
                  hMpvMap->Fill(diagX, diagY, r.mpv);
                  hStatusMap->Fill(diagX, diagY, 1.0);
                  hMpvDist->Fill(r.mpv);
                  hChi2NdfDist->Fill(r.chi2ndf);
                }
              } else {
                fout << layer << " " << chip << " " << chn << " " << r.mpv << " " << -1 << " " << r.width << " "
                     << r.chi2ndf << " " << hmips.Integral() << "\n";
                if (hMpvMap != nullptr) {
                  hMpvMap->Fill(diagX, diagY, r.mpv);
                  hStatusMap->Fill(diagX, diagY, 3.0);
                  hMpvDist->Fill(r.mpv);
                }
              }
            } else {
              // Fit failed (ndf<=0): instead of masking immediately, check
              // the combined histogram's peak bin as a rough MPV estimate
              // (the same idea as the pre-Gaudi Python
              // Calibration._compute_mips histogram-peak method). Accepted
              // only if positive and within MipFallbackMaxAdc; otherwise
              // fall back to masking as before.
              const int peakBin = hmips.GetMaximumBin();
              const double peakValue = hmips.GetXaxis()->GetBinCenter(peakBin);
              if (peakValue > 0.0 && peakValue <= m_mipFallbackMaxAdc.value()) {
                fout << layer << " " << chip << " " << chn << " " << peakValue << " " << -2 << " " << 0 << " " << 0
                     << " " << hmips.Integral() << "\n";
                if (hMpvMap != nullptr) {
                  hMpvMap->Fill(diagX, diagY, peakValue);
                  hStatusMap->Fill(diagX, diagY, 2.0);
                  hMpvDist->Fill(peakValue);
                }
              } else {
                fout << layer << " " << chip << " " << chn << " " << 0 << " " << -5 << " " << 0 << " " << 0 << " "
                     << 0 << "\n";
                if (hMpvMap != nullptr) {
                  hMpvMap->Fill(diagX, diagY, 0.0);
                  hStatusMap->Fill(diagX, diagY, 0.0);
                }
              }
            }
          } else {
            fout << layer << " " << chip << " " << chn << " " << 0 << " " << -10 << " " << 0 << " " << 0 << " " << 0
                 << "\n";
            if (hMpvMap != nullptr) {
              hMpvMap->Fill(diagX, diagY, 0.0);
              hStatusMap->Fill(diagX, diagY, 0.0);
            }
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
