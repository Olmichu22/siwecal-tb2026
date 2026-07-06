/*
 * PadMapGeometry: pure C++ (no ROOT/Gaudi dependency) port of
 * siwecal_eventbuilder/pad_map.py (transverse chip/channel -> x/y lookup)
 * and siwecal_eventbuilder/geometry.py's slab_z/slab_x0 (longitudinal
 * position, cumulative W radiation lengths).
 */
#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace k4siwecal {

// Port of pad_map.py::load_mapping_file: "chip x0 y0 channel x y" (one row
// per (chip,channel); x0/y0 are parsed by the reference format but unused --
// only chip/channel/x/y are kept). Skips blank lines, "#" comments, and the
// header row (detected the same way as Python: the first token must parse as
// an integer, optionally negative).
inline std::map<std::pair<int, int>, std::pair<float, float>> loadPadMappingFile(const std::string& path) {
  std::map<std::pair<int, int>, std::pair<float, float>> mapping;
  std::ifstream fin(path);
  if (!fin.is_open()) return mapping;
  std::string line;
  while (std::getline(fin, line)) {
    std::size_t start = line.find_first_not_of(" \t\r\n");
    if (start == std::string::npos || line[start] == '#') continue;
    std::istringstream iss(line);
    std::vector<std::string> items;
    std::string tok;
    while (iss >> tok) items.push_back(tok);
    if (items.size() < 6) continue;
    // items[0].lstrip("-").isdigit() guard: skip the header row ("chip ...").
    std::string first = items[0];
    std::size_t digitsStart = (!first.empty() && first[0] == '-') ? 1 : 0;
    if (digitsStart >= first.size()) continue;
    bool allDigits = true;
    for (std::size_t i = digitsStart; i < first.size(); ++i) {
      if (!std::isdigit(static_cast<unsigned char>(first[i]))) {
        allDigits = false;
        break;
      }
    }
    if (!allDigits) continue;

    const int chip = std::stoi(items[0]);
    const int channel = std::stoi(items[3]);
    const float x = std::stof(items[4]);
    const float y = std::stof(items[5]);
    mapping[{chip, channel}] = {x, y};
  }
  return mapping;
}

// Port of pad_map.py::PadMap. `position()` is looked up with the *geometric*
// chip index (0..15), NOT the hardware chip_id written to hit_chip -- see
// EventBuilder.h's HitCollector for where the two are kept deliberately
// distinct.
class PadMap {
 public:
  static PadMap fromFiles(const std::string& defaultFile, const std::map<int, std::string>& slabOverrides) {
    PadMap map;
    map.m_default = loadPadMappingFile(defaultFile);
    for (const auto& [slab, path] : slabOverrides) {
      map.m_slabMaps[slab] = loadPadMappingFile(path);
    }
    return map;
  }

  std::pair<float, float> position(int slab, int chip, int channel) const {
    const auto slabIt = m_slabMaps.find(slab);
    const auto& mapping = (slabIt != m_slabMaps.end()) ? slabIt->second : m_default;
    const auto it = mapping.find({chip, channel});
    if (it == mapping.end()) {
      return {std::numeric_limits<float>::quiet_NaN(), std::numeric_limits<float>::quiet_NaN()};
    }
    return it->second;
  }

 private:
  std::map<std::pair<int, int>, std::pair<float, float>> m_default;
  std::map<int, std::map<std::pair<int, int>, std::pair<float, float>>> m_slabMaps;
};

// Defaults mirror geometry.py::DEFAULT_SLAB_Z_MM / DEFAULT_SLAB_W_THICKNESS_MM
// (TB2026CERN 15-slab configuration). The live mappings/slab_z_positions.yml
// (note: NEGATIVE z values, 0..-165mm) is the actual source of truth used by
// the current pipeline and overrides these defaults when provided -- do not
// "fix" the sign to match these positive defaults.
inline constexpr double kWX0Mm = 3.5;  // radiation length of tungsten [mm]
inline constexpr std::array<double, 15> kDefaultSlabZMm = {0.0,  11.0,  22.0,  33.0,  44.0,  55.0,  66.0, 77.0,
                                                            88.0, 99.0,  110.0, 132.0, 143.0, 154.0, 165.0};
inline constexpr std::array<double, 15> kDefaultSlabWThicknessMm = {2.8, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2,
                                                                     4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6};

class SlabGeometry {
 public:
  static SlabGeometry defaults() {
    SlabGeometry g;
    g.m_slabZMm.assign(kDefaultSlabZMm.begin(), kDefaultSlabZMm.end());
    g.m_slabWThicknessMm.assign(kDefaultSlabWThicknessMm.begin(), kDefaultSlabWThicknessMm.end());
    return g;
  }

  // Port of geometry.py::load_slab_z_mm + load_slab_w_thickness_mm: a
  // hand-rolled line parser for slab_z_positions.yml's trivial fixed
  // structure (two top-level keys, each a flat YAML list of floats with
  // trailing "# comment"s) -- avoids adding a yaml-cpp dependency for a
  // two-key, no-nesting file. Falls back to `defaults()` field-by-field if a
  // key is absent from the file (matches load_slab_z_mm/
  // load_slab_w_thickness_mm each defaulting to an empty tuple, which
  // DetectorGeometry then... note: unlike Python (which would leave an empty
  // tuple and return NaN via out-of-range slab_z()/slab_x0()), this port
  // falls back to the compiled-in defaults per-field so a partially-written
  // override file doesn't silently NaN every hit -- flagged for verification
  // against the live file, which is expected to always provide both keys.
  static SlabGeometry fromYamlFile(const std::string& path) {
    SlabGeometry g = defaults();
    std::ifstream fin(path);
    if (!fin.is_open()) return g;

    std::vector<double>* current = nullptr;
    std::vector<double> parsedZ, parsedW;
    std::string line;
    while (std::getline(fin, line)) {
      const auto hashPos = line.find('#');
      std::string content = (hashPos == std::string::npos) ? line : line.substr(0, hashPos);
      const std::size_t start = content.find_first_not_of(" \t\r\n");
      if (start == std::string::npos) continue;
      const std::size_t end = content.find_last_not_of(" \t\r\n");
      const std::string trimmed = content.substr(start, end - start + 1);

      if (trimmed.rfind("slab_z_mm:", 0) == 0) {
        current = &parsedZ;
        continue;
      }
      if (trimmed.rfind("w_thickness_mm:", 0) == 0) {
        current = &parsedW;
        continue;
      }
      if (!trimmed.empty() && trimmed[0] == '-' && current != nullptr) {
        try {
          current->push_back(std::stod(trimmed.substr(1)));
        } catch (const std::exception&) {
          // malformed line: skip, keep whatever was parsed so far
        }
      }
    }
    if (!parsedZ.empty()) g.m_slabZMm = parsedZ;
    if (!parsedW.empty()) g.m_slabWThicknessMm = parsedW;
    return g;
  }

  double slabZ(int slab) const {
    if (slab >= 0 && static_cast<std::size_t>(slab) < m_slabZMm.size()) return m_slabZMm[slab];
    return std::numeric_limits<double>::quiet_NaN();
  }

  // Cumulative radiation lengths of W traversed up to and including `slab`.
  double slabX0(int slab) const {
    if (slab < 0 || static_cast<std::size_t>(slab) >= m_slabWThicknessMm.size()) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    double sum = 0.0;
    for (int i = 0; i <= slab; ++i) sum += m_slabWThicknessMm[i];
    return sum / kWX0Mm;
  }

 private:
  std::vector<double> m_slabZMm;
  std::vector<double> m_slabWThicknessMm;
};

}  // namespace k4siwecal
