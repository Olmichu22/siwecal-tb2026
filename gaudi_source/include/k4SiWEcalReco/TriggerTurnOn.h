/*
 * TriggerTurnOn: the trigger turn-on of a channel, chip or slab, measured from
 * the calibration's own histograms and fitted the way the simulation uses it.
 *
 * The calibration Fill pass (PedestalMipCalibrator Mode=Fill) splits every
 * selected readout by its fast-shaper bit: hitbit_high == 0 into ped_high,
 * hitbit_high == 1 into mip_high, per (slab, chip, channel, SCA) in raw ADC.
 * After subtracting each SCA's pedestal, mip / (ped + mip) bin by bin IS the
 * probability that the discriminator fired at that slow-shaper amplitude --
 * the turn-on RealDigitizer (siwecal_k4sim) has to reproduce:
 *
 *     P(x) = plateau * 0.5 * (1 + erf((x - mu) / (sqrt(2) * sigma)))
 *
 * mu = discriminator level, sigma = spread of the fast channel, plateau =
 * efficiency above threshold, all in ADC above the pedestal. The fit is a
 * binomial maximum likelihood over every bin with entries
 * (TBinomialEfficiencyFitter): a chi2 fit with a minimum-entries cut drops the
 * sparse plateau bins of a low-statistics cell and can settle on a low plateau
 * with a narrow turn-on.
 *
 * Gaudi-independent (ROOT only), kept apart from PedestalMipCalibrator.cpp
 * like PedestalMipCalib.h.
 */
#pragma once

#include "TBinomialEfficiencyFitter.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TH1D.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace k4siwecal {

// Amplitude axis of the turn-on histograms: ADC above the pedestal.
constexpr int kTurnOnNBins = 200;
constexpr double kTurnOnLo = -20.0;
constexpr double kTurnOnHi = 180.0;
// Fit range: below the turn-on (where P = 0 pins the baseline) to well onto the plateau.
constexpr double kTurnOnFitLo = -5.0;
constexpr double kTurnOnFitHi = 150.0;

inline double turnOnModel(double x, double mu, double sigma, double plateau) {
  return plateau * 0.5 * (1.0 + TMath::Erf((x - mu) / (std::sqrt(2.0) * sigma)));
}

// The fitted function keeps P inside [1e-7, 1 - 1e-7]: a stray hit-bit readout
// far below the threshold (P = 0 exactly) would otherwise make the binomial
// log-likelihood infinite and abort the minimisation.
inline Double_t turnOnTF1(Double_t* x, Double_t* p) {
  const double v = turnOnModel(x[0], p[0], p[1], p[2]);
  return std::min(1.0 - 1e-7, std::max(1e-7, v));
}

struct TurnOnFit {
  // status: 0 = own fit; 1 = too few entries; 2 = fit failed. The caller
  // decides the fallback (channel -> chip -> slab -> all slabs).
  int status = 1;
  double mu = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double plateau = std::numeric_limits<double>::quiet_NaN();
  double nPass = 0.0;   // readouts with the hit bit, fit range
  double nTotal = 0.0;  // all readouts, fit range
};

// Raw starting values from the efficiency itself: plateau = efficiency well
// above the half-height point, mu = first bin over half the plateau.
inline void turnOnStart(const TH1D& pass, const TH1D& total, double& mu0, double& plateau0) {
  mu0 = std::numeric_limits<double>::quiet_NaN();
  plateau0 = std::numeric_limits<double>::quiet_NaN();
  double p = 0.0, t = 0.0;
  for (int b = 1; b <= total.GetNbinsX(); ++b) {
    const double x = total.GetBinCenter(b);
    if (x > 40.0 && x < kTurnOnFitHi) {
      p += pass.GetBinContent(b);
      t += total.GetBinContent(b);
    }
  }
  plateau0 = t > 0.0 ? p / t : 0.95;
  if (!(plateau0 > 0.05)) plateau0 = 0.95;
  for (int b = 1; b <= total.GetNbinsX(); ++b) {
    const double n = total.GetBinContent(b);
    if (n < 5.0 || total.GetBinCenter(b) < 0.0) continue;
    if (pass.GetBinContent(b) / n > 0.5 * plateau0) {
      mu0 = total.GetBinCenter(b);
      break;
    }
  }
  if (!std::isfinite(mu0)) mu0 = 20.0;
}

// Fit one turn-on. With fixSigma/fixPlateau finite, only mu is free (a
// channel inherits its chip's width and plateau; three free parameters on a
// single channel's statistics are noise).
inline TurnOnFit fitTurnOn(TH1D& pass, TH1D& total, double minPass, double minTotal,
                           double fixSigma = std::numeric_limits<double>::quiet_NaN(),
                           double fixPlateau = std::numeric_limits<double>::quiet_NaN()) {
  TurnOnFit r;
  const int b0 = total.FindBin(kTurnOnFitLo), b1 = total.FindBin(kTurnOnFitHi);
  r.nPass = pass.Integral(b0, b1);
  r.nTotal = total.Integral(b0, b1);
  if (r.nPass < minPass || r.nTotal < minTotal) {
    r.status = 1;
    return r;
  }
  double mu0, plateau0;
  turnOnStart(pass, total, mu0, plateau0);

  TF1 f("turnon", turnOnTF1, kTurnOnFitLo, kTurnOnFitHi, 3);
  f.SetParameters(mu0, 3.0, plateau0);
  f.SetParLimits(0, 0.0, kTurnOnFitHi);
  f.SetParLimits(1, 0.3, 50.0);
  f.SetParLimits(2, 0.0, 1.0);
  if (std::isfinite(fixSigma)) f.FixParameter(1, fixSigma);
  if (std::isfinite(fixPlateau)) f.FixParameter(2, fixPlateau);

  TBinomialEfficiencyFitter fitter(&pass, &total);
  const TFitResultPtr res = fitter.Fit(&f, "RSQ");
  const int fitStatus = res;  // TFitResultPtr converts to the minimizer status
  if (fitStatus != 0 || !std::isfinite(f.GetParameter(0))) {
    r.status = 2;
    return r;
  }
  r.status = 0;
  r.mu = f.GetParameter(0);
  r.sigma = f.GetParameter(1);
  r.plateau = f.GetParameter(2);
  return r;
}

}  // namespace k4siwecal
