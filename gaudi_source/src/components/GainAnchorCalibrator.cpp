/*
 * GainAnchorCalibrator: the LG->HG anchor calibration as a Gaudi component.
 *
 * Above the high-gain switch (pedestal-subtracted adc_high >= 1500) the event
 * builder takes a hit's energy from the low gain, through the anchor line
 *
 *     adc_low - ped_lg = k * (adc_high - ped_hg) + c       (inverted per hit)
 *
 * (EcalEventBuilder GainRatio / GainIntercept). This algorithm MEASURES that
 * line from decoded data and writes the anchor file the reconstruction reads
 * (<calib-dir>/anchor/th<N>/gain_anchor_th<N>.txt, see
 * gaudi_jobs/condor/generate_reco_dag.py::_read_anchor). Until it existed the
 * anchor was fitted by hand with diagnostics/plot_gain_anchor_fit.py and its
 * numbers copied into that file; the method here is that script's, unchanged
 * (GainAnchorFit.h), so the deployed files are reproduced.
 *
 * Input: siwecaldecoded chunks of an ELECTRON run (muons never reach the
 * fiducial band: a MIP is ~19 high-gain ADC), and the high- and low-gain
 * pedestal tables the reconstruction uses. Every fired channel (hitbit_high)
 * of a live SCA with a usable pedestal in both gains fills the 2D histogram;
 * the line goes through the slice medians in the fiducial band.
 *
 * Optional occupancy window (MinHitBits / MaxHitBits, hit bits in the
 * chip-SCA): the chip baseline drops with the charge it collects and the
 * anchor slope with it (0.100-0.103 in quiet chips, 0.088 at >= 30 fired
 * channels, th210), so an anchor meant for the saturated hits -- which sit at
 * 12-24 fired channels -- can be traced on that population only. The default
 * window (0-64) takes every hit, as the deployed anchors did.
 *
 * Plain Gaudi::Algorithm, not a k4FWCore Consumer: the input is the
 * siwecaldecoded tree, not an EDM4hep collection in the event store, and the
 * result is a reduction over the whole run, final only after every entry is
 * seen. All work happens in initialize(), the pattern of PedestalMipCalibrator
 * and the other stages upstream of EcalToEDM4hep; execute() is empty.
 */
#include "k4SiWEcalReco/GainAnchorFit.h"

#include "Gaudi/Algorithm.h"
#include "Gaudi/Property.h"

#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TNamed.h"
#include "TParameter.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

using k4siwecal::kAnchorChannels;
using k4siwecal::kAnchorChips;
using k4siwecal::kAnchorScas;
using k4siwecal::kAnchorSlabs;

struct GainAnchorCalibrator final : Gaudi::Algorithm {
  using Gaudi::Algorithm::Algorithm;

  StatusCode initialize() override {
    TH1::AddDirectory(kFALSE);
    if (m_inputFiles.value().empty()) return fail("No InputFiles provided");
    if (m_pedestalFile.value().empty() || m_pedestalFileLg.value().empty())
      return fail("PedestalFile and PedestalFileLowGain are both required");
    if (m_outputFile.value().empty()) return fail("OutputFile is required");
    if (m_fidLow.value() >= m_fidHigh.value()) return fail("FiducialLow must be below FiducialHigh");
    if (m_minHitBits.value() > m_maxHitBits.value()) return fail("MinHitBits must not exceed MaxHitBits");

    std::string err;
    k4siwecal::PedestalGrid pedHg, pedLg;
    if (!pedHg.read(m_pedestalFile.value(), err) || !pedLg.read(m_pedestalFileLg.value(), err)) return fail(err);

    k4siwecal::GainAnchorHistogram hist;
    long long nAcq = 0, nHits = 0;
    const StatusCode sc = fill(pedHg, pedLg, hist, nAcq, nHits);
    if (!sc.isSuccess()) return sc;

    const auto fit = k4siwecal::fitGainAnchor(hist, m_fidLow.value(), m_fidHigh.value(), m_minSliceEntries.value());
    if (!fit.ok) return fail("anchor fit failed: " + fit.error + " (" + std::to_string(nHits) + " hits)");

    // nHits counts every selected fired channel; the histogram keeps those inside its range (adc_high -50..2400,
    // adc_low -25..320), which is the number the reference script reports and the one written as n_hits.
    const long long nOutside = nHits - hist.entries();
    nHits = hist.entries();
    info() << "GainAnchorCalibrator: " << nAcq << " acquisitions, " << nHits << " fired channels in the "
           << "histogram range (" << nOutside << " outside); "
           << fit.sliceX.size() << " slices in [" << m_fidLow.value() << ", " << m_fidHigh.value() << "]" << endmsg;
    info() << "GainAnchorCalibrator: adc_low - ped_lg = " << fit.k << " * (adc_high - ped_hg) + " << fit.c
           << "   (1/k = " << 1.0 / fit.k << ")" << endmsg;

    const StatusCode wsc = writeAnchorFile(fit, nAcq, nHits, nOutside);
    if (!wsc.isSuccess()) return wsc;
    if (!m_diagnosticsFile.value().empty()) return writeDiagnostics(hist, fit);
    return StatusCode::SUCCESS;
  }

  StatusCode execute(const EventContext&) const override { return StatusCode::SUCCESS; }

 private:
  StatusCode fail(const std::string& msg) const {
    error() << "GainAnchorCalibrator: " << msg << endmsg;
    return StatusCode::FAILURE;
  }

  // Branch buffers: ~7.4 MB, heap-allocated (worker-thread stacks are small).
  struct Buffers {
    int slboardId[kAnchorSlabs];
    int nColumns[kAnchorSlabs][kAnchorChips];
    int adcHigh[kAnchorSlabs][kAnchorChips][kAnchorScas][kAnchorChannels];
    int adcLow[kAnchorSlabs][kAnchorChips][kAnchorScas][kAnchorChannels];
    int hitbitHigh[kAnchorSlabs][kAnchorChips][kAnchorScas][kAnchorChannels];
  };

  StatusCode fill(const k4siwecal::PedestalGrid& pedHg, const k4siwecal::PedestalGrid& pedLg,
                  k4siwecal::GainAnchorHistogram& hist, long long& nAcq, long long& nHits) const {
    TChain chain(m_treeName.value().c_str());
    for (const auto& f : m_inputFiles.value()) {
      if (chain.Add(f.c_str(), 0) == 0) return fail("cannot add input file " + f);
    }
    const Long64_t nEntries = chain.GetEntries();
    if (nEntries <= 0) return fail("no entries in the input chunks");
    auto buf = std::make_unique<Buffers>();
    chain.SetBranchStatus("*", 0);
    for (const char* b : {"slboard_id", "nColumns", "adc_high", "adc_low", "hitbit_high"}) chain.SetBranchStatus(b, 1);
    chain.SetBranchAddress("slboard_id", buf->slboardId);
    chain.SetBranchAddress("nColumns", buf->nColumns);
    chain.SetBranchAddress("adc_high", buf->adcHigh);
    chain.SetBranchAddress("adc_low", buf->adcLow);
    chain.SetBranchAddress("hitbit_high", buf->hitbitHigh);

    const Long64_t stop = m_maxAcq.value() > 0 ? std::min<Long64_t>(nEntries, m_maxAcq.value()) : nEntries;
    const int minBits = m_minHitBits.value(), maxBits = m_maxHitBits.value();
    const bool occupancyCut = (minBits > 0 || maxBits < kAnchorChannels);
    for (Long64_t e = 0; e < stop; ++e) {
      if (chain.GetEntry(e) <= 0) return fail("read error at entry " + std::to_string(e));
      ++nAcq;
      for (int slot = 0; slot < kAnchorSlabs; ++slot) {
        const int slab = buf->slboardId[slot];  // the pedestal tables are indexed by slab
        if (slab < 0 || slab >= kAnchorSlabs) continue;
        for (int chip = 0; chip < kAnchorChips; ++chip) {
          const int ncol = buf->nColumns[slot][chip];
          for (int sca = 0; sca < kAnchorScas && sca < ncol; ++sca) {
            const int* hb = buf->hitbitHigh[slot][chip][sca];
            if (occupancyCut) {
              int bits = 0;
              for (int ch = 0; ch < kAnchorChannels; ++ch) bits += (hb[ch] == 1);
              if (bits < minBits || bits > maxBits) continue;
            }
            for (int ch = 0; ch < kAnchorChannels; ++ch) {
              if (hb[ch] != 1) continue;
              const int ah = buf->adcHigh[slot][chip][sca][ch];
              const int al = buf->adcLow[slot][chip][sca][ch];
              if (ah <= 0 || al <= 0) continue;
              const double ph = pedHg.get(slab, chip, sca, ch);
              const double pl = pedLg.get(slab, chip, sca, ch);
              if (ph <= 0.0 || pl <= 0.0) continue;
              hist.fill(ah - ph, al - pl);
              ++nHits;
            }
          }
        }
      }
    }
    return StatusCode::SUCCESS;
  }

  StatusCode writeAnchorFile(const k4siwecal::GainAnchorFitResult& fit, long long nAcq, long long nHits,
                             long long nOutside) const {
    std::ofstream out(m_outputFile.value());
    if (!out) return fail("cannot write " + m_outputFile.value());
    char now[32];
    const std::time_t t = std::time(nullptr);
    std::strftime(now, sizeof now, "%Y-%m-%d %H:%M", std::localtime(&t));
    const std::string th = m_thresholdLabel.value().empty() ? "?" : m_thresholdLabel.value();
    out << "# LG -> HG anchoring line for th" << th << " (SiW-ECAL TB2026CERN), measured by GainAnchorCalibrator\n"
        << "# (gaudi_source/src/components/GainAnchorCalibrator.cpp) on " << now << ".\n"
        << "# Maps a low-gain reading back onto the high gain in the fiducial region:\n"
        << "#     adc_low - ped_lg = k * (adc_high - ped_hg) + c\n"
        << "# The event builder inverts it for saturated hits and divides by MIP_hg:\n"
        << "#     adc_high_equiv = (adc_low - ped_lg - c) / k ;  energy = adc_high_equiv / MIP_hg\n"
        << "# A hit switches to the low gain when (adc_high - ped_hg) >= switch (pedestal-subtracted ADC).\n"
        << "# Fitted by the MEDIAN of adc_low per adc_high slice in the fiducial band\n"
        << "# (adc_high - ped in [" << m_fidLow.value() << ", " << m_fidHigh.value() << "]), " << fit.sliceX.size()
        << " slices of 12.5 ADC with >= " << m_minSliceEntries.value() << " entries.\n"
        << "# Fired channels (hitbit_high) of live SCAs, " << nHits << " of them in " << nAcq << " acquisitions;"
        << " hit bits per chip-SCA in [" << m_minHitBits.value() << ", " << m_maxHitBits.value() << "].\n"
        << "# Pedestals: " << m_pedestalFile.value() << "\n"
        << "#            " << m_pedestalFileLg.value() << "\n";
    char kbuf[32], cbuf[32];
    std::snprintf(kbuf, sizeof kbuf, "%.5f", fit.k);
    std::snprintf(cbuf, sizeof cbuf, "%.3f", fit.c);
    out << "k " << kbuf << "\n"
        << "c " << cbuf << "\n"
        << "switch " << m_switchAdc.value() << "\n"
        << "source_run " << (m_sourceLabel.value().empty() ? "unknown" : m_sourceLabel.value()) << "\n"
        << "n_hits " << nHits << "\n"
        << "n_hits_outside_range " << nOutside << "\n"
        << "n_acquisitions " << nAcq << "\n"
        << "n_slices " << fit.sliceX.size() << "\n"
        << "min_hit_bits " << m_minHitBits.value() << "\n"
        << "max_hit_bits " << m_maxHitBits.value() << "\n";
    info() << "GainAnchorCalibrator: wrote " << m_outputFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  StatusCode writeDiagnostics(const k4siwecal::GainAnchorHistogram& h, const k4siwecal::GainAnchorFitResult& fit) const {
    std::unique_ptr<TFile> f(TFile::Open(m_diagnosticsFile.value().c_str(), "RECREATE"));
    if (!f || f->IsZombie()) return fail("cannot create diagnostics file " + m_diagnosticsFile.value());
    using H = k4siwecal::GainAnchorHistogram;
    TH2F h2("anchor_hist", "fired channels;adc_high - ped_hg [ADC];adc_low - ped_lg [ADC]", H::kNX, H::kXLo,
            H::kXLo + H::kNX * H::kXStep, H::kNY, H::kYLo, H::kYLo + H::kNY * H::kYStep);
    for (int ix = 0; ix < H::kNX; ++ix)
      for (int iy = 0; iy < H::kNY; ++iy)
        if (const auto n = h.count(ix, iy)) h2.SetBinContent(ix + 1, iy + 1, static_cast<double>(n));
    std::vector<double> ax, am;
    for (int ix = 0; ix < H::kNX; ++ix) {
      const double m = h.sliceMedian(ix, m_minSliceEntries.value());
      if (!std::isnan(m)) {
        ax.push_back(H::xCentre(ix));
        am.push_back(m);
      }
    }
    TGraph gAll(static_cast<int>(ax.size()), ax.data(), am.data());
    gAll.SetName("slice_medians");
    gAll.SetTitle("median adc_low per adc_high slice;adc_high - ped_hg [ADC];adc_low - ped_lg [ADC]");
    TGraph gFit(static_cast<int>(fit.sliceX.size()), fit.sliceX.data(), fit.sliceMedian.data());
    gFit.SetName("fit_points");
    gFit.SetTitle("slices in the fiducial band (fitted)");
    TParameter<double> pk("k", fit.k), pc("c", fit.c), plo("fiducial_low", m_fidLow.value()),
        phi("fiducial_high", m_fidHigh.value()), psw("switch", m_switchAdc.value());
    TNamed src("source_run", m_sourceLabel.value().c_str());
    f->cd();
    h2.Write();
    gAll.Write();
    gFit.Write();
    pk.Write();
    pc.Write();
    plo.Write();
    phi.Write();
    psw.Write();
    src.Write();
    f->Close();
    info() << "GainAnchorCalibrator: wrote diagnostics to " << m_diagnosticsFile.value() << endmsg;
    return StatusCode::SUCCESS;
  }

  Gaudi::Property<std::vector<std::string>> m_inputFiles{this, "InputFiles", {},
                                                         "siwecaldecoded chunks of an electron run"};
  Gaudi::Property<std::string> m_treeName{this, "TreeName", "siwecaldecoded", "Input TTree name"};
  Gaudi::Property<std::string> m_pedestalFile{this, "PedestalFile", "",
                                              "High-gain pedestal table (PedestalMipCalibrator format), the one "
                                              "the reconstruction uses"};
  Gaudi::Property<std::string> m_pedestalFileLg{this, "PedestalFileLowGain", "", "Low-gain pedestal table"};
  Gaudi::Property<std::string> m_outputFile{this, "OutputFile", "",
                                            "Anchor file to write (gain_anchor_th<N>.txt: 'k', 'c', 'switch', "
                                            "'source_run' lines, read by generate_reco_dag.py)"};
  Gaudi::Property<std::string> m_diagnosticsFile{this, "DiagnosticsFile", "",
                                                 "Optional ROOT file: the 2D histogram, the slice medians, the "
                                                 "fitted points and k, c"};
  Gaudi::Property<long long> m_maxAcq{this, "MaxAcquisitions", 8000,
                                      "Acquisitions read, in InputFiles order (0 = all). 8000 is what the "
                                      "deployed anchors used"};
  Gaudi::Property<double> m_fidLow{this, "FiducialLow", 150.0,
                                   "Lower edge of the fiducial band, pedestal-subtracted adc_high: above the "
                                   "low-gain noise"};
  Gaudi::Property<double> m_fidHigh{this, "FiducialHigh", 1200.0,
                                    "Upper edge: below the switch and the onset of high-gain non-linearity "
                                    "(1 % at raw adc_high ~1270)"};
  Gaudi::Property<long long> m_minSliceEntries{this, "MinEntriesPerSlice", 50,
                                               "adc_high slices with fewer fired channels are not used"};
  Gaudi::Property<int> m_minHitBits{this, "MinHitBits", 0,
                                    "Occupancy window: hit bits in the chip-SCA, lower edge (0 = no cut)"};
  Gaudi::Property<int> m_maxHitBits{this, "MaxHitBits", 64, "Occupancy window, upper edge (64 = no cut)"};
  Gaudi::Property<double> m_switchAdc{this, "SwitchAdc", 1500.0,
                                      "Written to the file as 'switch': EcalEventBuilder.AdcSaturationThreshold"};
  Gaudi::Property<std::string> m_sourceLabel{this, "SourceLabel", "", "Run name written as 'source_run'"};
  Gaudi::Property<std::string> m_thresholdLabel{this, "ThresholdLabel", "", "th<N> label, for the file header"};
};

DECLARE_COMPONENT(GainAnchorCalibrator)
