/*
 * CalibrationTables: pure C++ (no ROOT/Gaudi dependency) port of the
 * read-side of siwecal_eventbuilder/calibration.py::Calibration --
 * `from_files` plus the `pedestal()`/`mip()`/`is_masked()` lookups applied
 * per hit by hit_collector.py. The GENERATION side (`_compute_pedestals`/
 * `_compute_mips`, i.e. producing these text tables from data) is already
 * ported separately in PedestalMipCalib.h / PedestalMipCalibrator.cpp; this
 * header only consumes their output format.
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace k4siwecal {

class CalibrationTables {
 public:
  // Default-constructed state is equivalent to disabled() (m_enabled=false),
  // so a `CalibrationTables calib;` local can be assigned via disabled()/
  // fromFiles() afterwards without needing a two-phase factory dance.
  CalibrationTables() = default;

  // Raw-ADC mode: pedestal()=0, mip()=1, isMasked()=false always. Port of
  // Calibration.disabled().
  static CalibrationTables disabled() {
    CalibrationTables t;
    t.m_enabled = false;
    return t;
  }

  // Port of Calibration.from_files: loads both text tables, unions their
  // masked-channel sets (a channel masked in EITHER file is masked in the
  // result). Sets `ok=false` and `error` on a malformed pedestal file --
  // calibration.py::_read_pedestal_file raises ValueError in that case
  // rather than silently substituting a value, and this port preserves that
  // hard-failure behaviour (the caller should fail the Gaudi component).
  //
  // The low-gain paths are OPTIONAL: pass them empty and hasLowGain() stays
  // false, pedestalLg()/mipLg() are never consulted, and the hit energy is
  // computed from the high gain exactly as before low-gain support existed.
  // A threshold whose low-gain tables have not been produced therefore keeps
  // working instead of silently reconstructing against absent calibration.
  static CalibrationTables fromFiles(const std::string& pedestalPath, const std::string& mipPath,
                                      double pedestalFallback, double defaultMipFallback, bool& ok,
                                      std::string& error, const std::string& pedestalPathLg = "",
                                      const std::string& mipPathLg = "") {
    CalibrationTables t;
    t.m_enabled = true;
    t.m_pedestalFallback = pedestalFallback;

    std::set<std::tuple<int, int, int>> pedMasked;
    if (!t.readPedestalFile(pedestalPath, t.m_pedestalMap, pedMasked, error)) {
      ok = false;
      return t;
    }

    std::set<std::tuple<int, int, int>> mipMasked;
    t.readMipFile(mipPath, defaultMipFallback, t.m_mipMap, t.m_defaultMip, mipMasked);

    t.m_masked = pedMasked;
    t.m_masked.insert(mipMasked.begin(), mipMasked.end());

    if (!pedestalPathLg.empty() && !mipPathLg.empty()) {
      std::set<std::tuple<int, int, int>> pedMaskedLg;
      if (!t.readPedestalFile(pedestalPathLg, t.m_pedestalMapLg, pedMaskedLg, error)) {
        ok = false;
        return t;
      }
      std::set<std::tuple<int, int, int>> mipMaskedLg;
      t.readMipFile(mipPathLg, defaultMipFallback, t.m_mipMapLg, t.m_defaultMipLg, mipMaskedLg);
      t.m_maskedLg = pedMaskedLg;
      t.m_maskedLg.insert(mipMaskedLg.begin(), mipMaskedLg.end());
      t.m_hasLowGain = true;
    }

    ok = true;
    return t;
  }

  bool enabled() const { return m_enabled; }

  // True once low-gain tables were loaded; false keeps the pre-low-gain
  // behaviour (energy always from the high gain).
  bool hasLowGain() const { return m_hasLowGain; }

  // Pedestal for a channel/SCA (0 in raw mode, fallback if unknown).
  double pedestal(int slabId, int chipId, int channel, int sca) const {
    if (!m_enabled) return 0.0;
    auto it = m_pedestalMap.find({slabId, chipId, channel, sca});
    return it != m_pedestalMap.end() ? it->second : m_pedestalFallback;
  }

  double pedestalLg(int slabId, int chipId, int channel, int sca) const {
    if (!m_enabled || !m_hasLowGain) return 0.0;
    auto it = m_pedestalMapLg.find({slabId, chipId, channel, sca});
    return it != m_pedestalMapLg.end() ? it->second : m_pedestalFallback;
  }

  // MIP scale for a channel (1 in raw mode, default-median if uncalibrated).
  double mip(int slabId, int chipId, int channel) const {
    if (!m_enabled) return 1.0;
    auto it = m_mipMap.find({slabId, chipId, channel});
    return it != m_mipMap.end() ? it->second : m_defaultMip;
  }

  double mipLg(int slabId, int chipId, int channel) const {
    if (!m_enabled || !m_hasLowGain) return 1.0;
    auto it = m_mipMapLg.find({slabId, chipId, channel});
    return it != m_mipMapLg.end() ? it->second : m_defaultMipLg;
  }

  // True for channels masked in EITHER the pedestal or MIP file of the HIGH
  // gain. Deliberately NOT unioned with the low-gain masks (which is what
  // calibration.py does): a channel whose low-gain MIP fit failed is still a
  // perfectly good high-gain channel, and unioning would zero out its energy
  // for every unsaturated hit -- i.e. the overwhelming majority of them.
  // Saturated hits on a low-gain-masked channel are handled in buildHit().
  bool isMasked(int slabId, int chipId, int channel) const {
    if (!m_enabled) return false;
    return m_masked.count({slabId, chipId, channel}) > 0;
  }

  // True when the LOW gain has no usable calibration for this channel, so a
  // saturated hit there cannot be recovered.
  bool isMaskedLg(int slabId, int chipId, int channel) const {
    if (!m_enabled || !m_hasLowGain) return false;
    return m_maskedLg.count({slabId, chipId, channel}) > 0;
  }

 private:
  // Port of _read_pedestal_file: "layer chip channel mean0 err0 wid0 mean1
  // err1 wid1 ... (15 SCAs)". Only the mean of each triplet is used (Python's
  // `items[3::3][:15]`: start at index 3, every 3rd item, first 15 values). A
  // pad's 15 means must be either all-zero (masked) or all finite+nonzero
  // (calibrated) -- anything else is a hard error, matching Python's
  // ValueError (no silent substitution for a malformed pad).
  bool readPedestalFile(const std::string& path, std::map<std::tuple<int, int, int, int>, double>& out,
                        std::set<std::tuple<int, int, int>>& masked, std::string& error) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
      error = "Cannot open pedestal file: " + path;
      return false;
    }
    std::string line;
    while (std::getline(fin, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream iss(line);
      std::vector<std::string> items;
      std::string tok;
      while (iss >> tok) items.push_back(tok);
      if (items.size() < 3) continue;
      const int layer = std::stoi(items[0]);
      const int chip = std::stoi(items[1]);
      const int channel = std::stoi(items[2]);

      std::vector<double> means;
      for (std::size_t i = 3; i < items.size() && means.size() < 15; i += 3) {
        means.push_back(std::stod(items[i]));
      }
      const bool allZero = std::all_of(means.begin(), means.end(), [](double m) { return m == 0.0; });
      if (allZero) {
        masked.insert({layer, chip, channel});
        continue;
      }
      const bool allValid =
          std::all_of(means.begin(), means.end(), [](double m) { return std::isfinite(m) && m != 0.0; });
      if (allValid) {
        for (std::size_t sca = 0; sca < means.size(); ++sca) {
          out[{layer, chip, channel, static_cast<int>(sca)}] = means[sca];
        }
        continue;
      }
      error = "Invalid pedestal calibration for pad (layer=" + std::to_string(layer) +
              ", chip=" + std::to_string(chip) + ", channel=" + std::to_string(channel) + ") in " + path +
              ": SCA means must be all zero (masked) or all finite and non-zero (calibrated).";
      return false;
    }
    return true;
  }

  // Port of _read_mip_file: "layer chip channel mpv ...". mpv>0 -> calibrated;
  // else masked. Default MIP = median of loaded mpv values (np.median
  // semantics: average the two middle values for an even count), or
  // defaultMipFallback if no channel was calibrated at all.
  void readMipFile(const std::string& path, double defaultMipFallback,
                    std::map<std::tuple<int, int, int>, double>& out, double& defaultMipOut,
                    std::set<std::tuple<int, int, int>>& masked) {
    std::ifstream fin(path);
    if (!fin.is_open()) {
      defaultMipOut = defaultMipFallback;
      return;
    }
    std::string line;
    std::vector<double> mpvValues;
    while (std::getline(fin, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream iss(line);
      int layer, chip, channel;
      double mpv;
      if (!(iss >> layer >> chip >> channel >> mpv)) continue;
      if (mpv > 0) {
        out[{layer, chip, channel}] = mpv;
        mpvValues.push_back(mpv);
      } else {
        masked.insert({layer, chip, channel});
      }
    }
    defaultMipOut = mpvValues.empty() ? defaultMipFallback : median(mpvValues);
  }

  static double median(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    const std::size_t n = values.size();
    if (n % 2 == 1) return values[n / 2];
    return 0.5 * (values[n / 2 - 1] + values[n / 2]);
  }

  bool m_enabled = false;
  bool m_hasLowGain = false;
  double m_pedestalFallback = 250.0;
  double m_defaultMip = 1.0;
  double m_defaultMipLg = 1.0;
  std::map<std::tuple<int, int, int, int>, double> m_pedestalMap, m_pedestalMapLg;
  std::map<std::tuple<int, int, int>, double> m_mipMap, m_mipMapLg;
  std::set<std::tuple<int, int, int>> m_masked, m_maskedLg;
};

}  // namespace k4siwecal
