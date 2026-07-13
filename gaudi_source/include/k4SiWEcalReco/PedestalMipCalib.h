/*
 * PedestalMipCalib: pedestal and MIP fit math ported from the external
 * reference tool SiWECAL-TB-analysis, ROOT-dependent (TH1F/TF1/TSpectrum)
 * but Gaudi-independent, kept separate from PedestalMipCalibrator.cpp so it
 * stays testable with synthetic histograms (same separation principle as
 * EcalShowerVars.h).
 *
 * Verified against the actually-deployed calibration files in this repo
 * (calibration/MuonCalib_it2/{pedestals,mips}/thNNN/Pedestal_or_MIP_*.txt):
 * their header lines
 * ("#pedestal results (fit to a gaussian) remove channels/sca with two
 * pedestals peaks from the analysis : PROTO15" and "#mip results
 * PROTO15-TB2022-03") match SiWECAL-TB-analysis/SLBperformance/TBchecks/
 * analysis.h's `pedanalysis` and `mipanalysis_summary` (pedestal_mode=1)
 * exactly -- that is the method ported here, NOT the simpler single-Gaussian
 * pedestal fit or the noise-covariance decomposition described in earlier
 * exploration notes (those come from a different reference tool variant
 * that does not match this repo's deployed files).
 *
 * Event selection and histogram binning (badbcid==0, nhits<=maxnhit,
 * SimpleCoincidenceTagger, 400 bins 100.5-500.5 for pedestal, 900 bins
 * 100.5-1000.5 for MIP) come from SiWECAL-TB-analysis/SLBperformance/
 * DecodedSLBAnalysis.cc::NSlabsAnalysis, which fills the histograms these
 * fits are applied to.
 */
#pragma once

#include "k4SiWEcalReco/SlbFrameDecoder.h"  // kSlbDepth, kSkirocsPerAsu, kScasInSkiroc

#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSpectrum.h"

#include <algorithm>
#include <cmath>

namespace k4siwecal {

// ---------------------------------------------------------------------------
// Coincidence tagger (SiWECAL-TB-analysis/include/utils.h:8-29, verbatim
// port). Single-pass: only needs the current tree entry's already-loaded
// bcid/badbcid arrays, no cross-event buffering. Counts how many OTHER
// slboards (besides `ilayer`) have at least one (chip, sca) with
// badbcid==0, bcid>=0 and |bcid-bcidRef|<2.
inline int simpleCoincidenceTagger(const int bcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc],
                                    const int badbcid[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc], int nSlboards,
                                    int ilayer, int bcidRef) {
  int bcidSeen = 0;
  for (int islb = 0; islb < nSlboards; ++islb) {
    if (islb == ilayer) continue;
    int seenThisSlb = 0;
    for (int isca = 0; isca < kScasInSkiroc; ++isca) {
      for (int ichip = 0; ichip < kSkirocsPerAsu; ++ichip) {
        if (badbcid[islb][ichip][isca] != 0 || bcid[islb][ichip][isca] < 0) continue;
        if (std::abs(bcid[islb][ichip][isca] - bcidRef) < 2) ++seenThisSlb;
      }
    }
    if (seenThisSlb > 0) ++bcidSeen;
  }
  return bcidSeen;
}

// ---------------------------------------------------------------------------
// Port of TBchecks/analysis.h::pedanalysis's per-(layer,chip,channel,sca)
// fit. Operates on one already-filled pedestal histogram (400 bins,
// 100.5-500.5, see DecodedSLBAnalysis.cc::NSlabsAnalysis).
struct PedestalFitResult {
  double mean = 0.;
  double error = 0.;
  double width = 0.;
  double widthError = 0.;
  // kGood: npeaks==1, dual-window Gaussian fit applied (this SCA's row uses
  //   the fit values directly).
  // kMultiPeak: histogram had >50 entries but TSpectrum found npeaks!=1;
  //   the caller falls back to the weighted average of the kGood SCAs for
  //   this channel, written with error sentinel -5 (utils.h/pedanalysis
  //   convention).
  // kLowStats: histogram had <=50 entries; same fallback, error sentinel -10.
  enum class Status { kGood, kMultiPeak, kLowStats } status = Status::kLowStats;
};

inline PedestalFitResult fitPedestalSca(TH1F& hist) {
  PedestalFitResult r;
  if (hist.GetEntries() <= 50) {
    r.status = PedestalFitResult::Status::kLowStats;
    return r;
  }

  TSpectrum spectrum;
  const int npeaks = spectrum.Search(&hist, 2, "", 0.8);
  if (npeaks != 1) {
    r.status = PedestalFitResult::Status::kMultiPeak;
    return r;
  }

  const double* peakX = spectrum.GetPositionX();
  const double* peakY = spectrum.GetPositionY();
  double meanPeakHigher = 0.;
  double meanHighHigher = 0.;
  for (int ipeak = 0; ipeak < npeaks; ++ipeak) {
    if (peakY[ipeak] > meanHighHigher && peakX[ipeak] > 180) {
      meanHighHigher = peakY[ipeak];
      meanPeakHigher = peakX[ipeak];
    }
  }

  TF1 f0("f0", "gaus", meanPeakHigher - 8, meanPeakHigher + 8);
  TF1 f1("f1", "gaus", meanPeakHigher - 16, meanPeakHigher + 16);
  hist.Fit(&f0, "RQE");
  hist.Fit(&f1, "RQE");

  const TF1& chosen = (f0.GetParameter(2) < f1.GetParameter(2)) ? f0 : f1;
  r.mean = chosen.GetParameter(1);
  r.width = chosen.GetParameter(2);
  r.error = chosen.GetParError(1);
  r.widthError = chosen.GetParError(2);
  r.status = PedestalFitResult::Status::kGood;
  return r;
}

// ---------------------------------------------------------------------------
// Landau (scale/MP/area) convolved with a Gaussian (sigma), CERNLIB-style
// numeric convolution. Verbatim port of utils.h::langaufun.
inline Double_t langaufun(Double_t* x, Double_t* par) {
  const Double_t invsq2pi = 0.3989422804014;
  const Double_t mpshift = -0.22278298;
  const Double_t np = 100.0;
  const Double_t sc = 5.0;

  const Double_t mpc = par[1] - mpshift * par[0];
  const Double_t xlow = x[0] - sc * par[3];
  const Double_t xupp = x[0] + sc * par[3];
  const Double_t step = (xupp - xlow) / np;

  Double_t sum = 0.0;
  for (Double_t i = 1.0; i <= np / 2; i += 1.0) {
    Double_t xx = xlow + (i - .5) * step;
    Double_t fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);

    xx = xupp - (i - .5) * step;
    fland = TMath::Landau(xx, mpc, par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0], xx, par[3]);
  }
  return par[2] * step * sum * invsq2pi / par[3];
}

// Port of TBchecks/analysis.h::resultfit (initial Landau-only fit for
// starting values, then the full Landau*Gauss fit via utils.h::langaufit).
// `highGain` selects the gain-dependent parameter limits/fit range exactly
// as resultfit() does for gain=="high" vs. other.
struct MipFitResult {
  double mpv = 0.;
  double empv = 0.;
  double width = 0.;
  double chi2ndf = 0.;
  int ndf = 0;
  bool ok = false;
};

inline MipFitResult fitLangau(TH1F& hist, bool highGain) {
  MipFitResult result;

  double pllo[4];
  double plhi[4];
  if (highGain) {
    pllo[0] = 0.1;
    pllo[1] = 14.;
    pllo[2] = 1.0;
    pllo[3] = 0.5;
    plhi[0] = 20.0;
    plhi[1] = 220.0;
    plhi[2] = 1.0e8;
    plhi[3] = 10.0;
  } else {
    pllo[0] = 0.1;
    pllo[1] = 2.;
    pllo[2] = 1.0;
    pllo[3] = 0.05;
    plhi[0] = 5.;
    plhi[1] = 15.0;
    plhi[2] = 1.0e8;
    plhi[3] = 2.0;
  }

  double fr[2];
  if (highGain) {
    // Fit over the full positive ADC range so the Langau describes the whole
    // MIP distribution (peak + Landau tail), not just a narrow window around
    // the mean. The earlier narrow range ([mean-0.8*RMS, mean*1.3]) biased the
    // fitted MPV low (it sat on the rising edge, left of the visible peak) and
    // made chi2/ndf meaningless -- computed over too few bins, so a fit that
    // missed the peak still scored a good chi2. A full-range fit fixes both.
    fr[0] = 0.;
    fr[1] = 100.;
    hist.GetXaxis()->SetRangeUser(fr[0], fr[1]);
  } else {
    fr[0] = std::max(hist.GetMean() - 2. * hist.GetRMS(), 2.);
    hist.GetXaxis()->SetRangeUser(fr[0], 25.);
    fr[1] = hist.GetMean() * 1.1 + 2. * hist.GetRMS();
  }

  TF1 fitlandau("fitlandau_tmp", "landau", fr[0], fr[1]);
  hist.Fit(&fitlandau, "QRE");
  hist.GetXaxis()->SetRangeUser(fr[0], fr[1]);

  double sv[4];
  sv[0] = std::min(std::max(fitlandau.GetParameter(2), pllo[0]), plhi[0]);
  sv[1] = std::min(std::max(fitlandau.GetParameter(1), pllo[1]), plhi[1]);
  sv[2] = hist.Integral("width") * 1.2;
  sv[3] = std::min(std::max(hist.GetRMS() * 0.1, pllo[3]), plhi[3]);

  TF1 langau("langau_tmp", langaufun, fr[0], fr[1], 4);
  langau.SetParameters(sv);
  langau.SetParNames("Width", "MP", "Area", "GSigma");
  for (int i = 0; i < 4; ++i) {
    langau.SetParLimits(i, pllo[i], plhi[i]);
  }
  hist.Fit(&langau, "RBQM");

  result.width = langau.GetParameter(0);
  result.mpv = langau.GetParameter(1);
  result.empv = langau.GetParError(1);
  result.ndf = langau.GetNDF();
  result.chi2ndf = result.ndf > 0 ? langau.GetChisquare() / result.ndf : 0.;
  result.ok = true;
  return result;
}

}  // namespace k4siwecal
