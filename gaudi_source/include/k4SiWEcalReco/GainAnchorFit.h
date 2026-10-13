/*
 * GainAnchorFit: pure C++ (no ROOT/Gaudi dependency) core of the LG->HG anchor
 * calibration, used by GainAnchorCalibrator. Kept ROOT/Gaudi-independent so it
 * can be unit-tested with synthetic data (tests/test_gain_anchor_fit.cpp), same
 * separation principle as SlbFrameDecoder.h and EcalShowerVars.h.
 *
 * The anchor maps a low-gain reading onto the high gain:
 *
 *     adc_low - ped_lg = k * (adc_high - ped_hg) + c
 *
 * EcalEventBuilder inverts it for every hit past the high-gain switch
 * (BuilderConfig::gainRatio / gainIntercept). It is traced, not least-squares
 * fitted to the hits: the fired channels fill a 2D histogram of
 * (adc_high - ped_hg, adc_low - ped_lg); in each adc_high slice inside the
 * fiducial band the MEDIAN adc_low is taken, and a straight line goes through
 * those medians. Least squares on the hits themselves is pulled off the band by
 * its tails (c = +3.13 where the medians say +1.60, see EventBuilder.h).
 *
 * This is a byte-for-byte port of diagnostics/plot_gain_anchor_fit.py, the
 * script the deployed anchor files were produced with: same bins (12.5 ADC in
 * adc_high from -50 to 2400, 0.5 ADC in adc_low from -25 to 320), same median
 * (the centre of the first adc_low bin where the cumulative count reaches half
 * the slice), same 50-entry floor per slice, same fiducial band (150-1200,
 * inclusive on the slice centres), same pedestal-table reading (a per-SCA mean
 * is kept only if 0 < mean < 500).
 */
#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace k4siwecal {

// Detector dimensions of the decoded tree and the pedestal tables.
inline constexpr int kAnchorSlabs = 15;
inline constexpr int kAnchorChips = 16;
inline constexpr int kAnchorScas = 15;
inline constexpr int kAnchorChannels = 64;

// Per-(slab, chip, sca, channel) pedestal means read from a
// PedestalMipCalibrator table ("layer chip channel" then mean/error/width per
// SCA). 0 = no usable pedestal, which excludes the cell.
class PedestalGrid {
 public:
  PedestalGrid() : m_mean(static_cast<std::size_t>(kAnchorSlabs) * kAnchorChips * kAnchorScas * kAnchorChannels, 0.0) {}

  bool read(const std::string& path, std::string& error) {
    std::ifstream in(path);
    if (!in) {
      error = "cannot open pedestal table " + path;
      return false;
    }
    std::string line;
    std::size_t rows = 0;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      std::vector<double> f;
      double v;
      while (ss >> v) f.push_back(v);
      if (f.size() < 4) continue;
      const int slab = static_cast<int>(f[0]), chip = static_cast<int>(f[1]), chn = static_cast<int>(f[2]);
      if (slab < 0 || slab >= kAnchorSlabs || chip < 0 || chip >= kAnchorChips || chn < 0 || chn >= kAnchorChannels)
        continue;
      for (int sca = 0; sca < kAnchorScas; ++sca) {
        const std::size_t col = 3 + 3 * static_cast<std::size_t>(sca);
        if (col >= f.size()) break;
        const double m = f[col];
        at(slab, chip, sca, chn) = (std::isfinite(m) && m > 0.0 && m < 500.0) ? m : 0.0;
      }
      ++rows;
    }
    if (rows == 0) {
      error = "no channel rows in pedestal table " + path;
      return false;
    }
    return true;
  }

  double get(int slab, int chip, int sca, int chn) const { return m_mean[index(slab, chip, sca, chn)]; }
  double& at(int slab, int chip, int sca, int chn) { return m_mean[index(slab, chip, sca, chn)]; }

 private:
  static std::size_t index(int slab, int chip, int sca, int chn) {
    return ((static_cast<std::size_t>(slab) * kAnchorChips + chip) * kAnchorScas + sca) * kAnchorChannels + chn;
  }
  std::vector<double> m_mean;
};

// The 2D histogram the anchor is traced on.
class GainAnchorHistogram {
 public:
  static constexpr double kXLo = -50.0;
  static constexpr double kXStep = 12.5;
  static constexpr int kNX = 196;  // -50 .. 2400
  static constexpr double kYLo = -25.0;
  static constexpr double kYStep = 0.5;
  static constexpr int kNY = 690;  // -25 .. 320

  GainAnchorHistogram() : m_counts(static_cast<std::size_t>(kNX) * kNY, 0) {}

  // x = adc_high - ped_hg, y = adc_low - ped_lg. Out-of-range entries are dropped.
  void fill(double x, double y) {
    const int ix = static_cast<int>(std::floor((x - kXLo) / kXStep));
    const int iy = static_cast<int>(std::floor((y - kYLo) / kYStep));
    if (ix < 0 || ix >= kNX || iy < 0 || iy >= kNY) return;
    ++m_counts[static_cast<std::size_t>(ix) * kNY + iy];
    ++m_entries;
  }

  std::int64_t count(int ix, int iy) const { return m_counts[static_cast<std::size_t>(ix) * kNY + iy]; }
  std::int64_t entries() const { return m_entries; }
  static double xCentre(int ix) { return kXLo + (ix + 0.5) * kXStep; }
  static double yCentre(int iy) { return kYLo + (iy + 0.5) * kYStep; }

  std::int64_t sliceEntries(int ix) const {
    std::int64_t t = 0;
    for (int iy = 0; iy < kNY; ++iy) t += count(ix, iy);
    return t;
  }

  // Median adc_low of one adc_high slice: the centre of the first bin whose
  // cumulative count reaches half the slice. NaN below minEntries.
  double sliceMedian(int ix, std::int64_t minEntries) const {
    const std::int64_t tot = sliceEntries(ix);
    if (tot < minEntries || tot == 0) return std::numeric_limits<double>::quiet_NaN();
    const double half = 0.5 * static_cast<double>(tot);
    std::int64_t cum = 0;
    for (int iy = 0; iy < kNY; ++iy) {
      cum += count(ix, iy);
      if (static_cast<double>(cum) >= half) return yCentre(iy);
    }
    return yCentre(kNY - 1);
  }

 private:
  std::vector<std::int64_t> m_counts;
  std::int64_t m_entries = 0;
};

struct GainAnchorFitResult {
  bool ok = false;
  std::string error;
  double k = std::numeric_limits<double>::quiet_NaN();
  double c = std::numeric_limits<double>::quiet_NaN();
  // The slice centres and medians the line goes through (fiducial band only).
  std::vector<double> sliceX, sliceMedian;
};

// Straight line through the slice medians whose centre lies in [fidLo, fidHi].
inline GainAnchorFitResult fitGainAnchor(const GainAnchorHistogram& h, double fidLo, double fidHi,
                                         std::int64_t minEntries) {
  GainAnchorFitResult r;
  for (int ix = 0; ix < GainAnchorHistogram::kNX; ++ix) {
    const double xc = GainAnchorHistogram::xCentre(ix);
    if (xc < fidLo || xc > fidHi) continue;
    const double m = h.sliceMedian(ix, minEntries);
    if (std::isnan(m)) continue;
    r.sliceX.push_back(xc);
    r.sliceMedian.push_back(m);
  }
  const std::size_t n = r.sliceX.size();
  if (n < 2) {
    r.error = "fewer than 2 adc_high slices with >= " + std::to_string(minEntries) + " entries in the fiducial band";
    return r;
  }
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (std::size_t i = 0; i < n; ++i) {
    sx += r.sliceX[i];
    sy += r.sliceMedian[i];
    sxx += r.sliceX[i] * r.sliceX[i];
    sxy += r.sliceX[i] * r.sliceMedian[i];
  }
  const double dn = static_cast<double>(n);
  const double den = dn * sxx - sx * sx;
  if (den == 0.0) {
    r.error = "degenerate slice positions";
    return r;
  }
  r.k = (dn * sxy - sx * sy) / den;
  r.c = (sy - r.k * sx) / dn;
  r.ok = std::isfinite(r.k) && std::isfinite(r.c) && r.k > 0.0;
  if (!r.ok) r.error = "non-physical line (k <= 0 or not finite)";
  return r;
}

}  // namespace k4siwecal
