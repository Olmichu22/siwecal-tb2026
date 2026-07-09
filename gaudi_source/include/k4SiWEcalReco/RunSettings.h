/*
 * RunSettings: pure C++ (no ROOT/Gaudi dependency) parser for the DAQ
 * `Run_Settings.txt` file written alongside each raw run, generalizing
 * siwecal_eventbuilder/run_settings.py::read_threshold_dac to the other
 * per-chip fields and the two once-per-file fields requested for the
 * siwecaldecoded tree.
 *
 * Per-chip fields (ThresholdDAC, HoldDelay, FSPeakTime,
 * GainSelectionThreshold) appear once per "## ChipIndex: ..." line -- one
 * per chip per slab, so potentially hundreds of occurrences per file. They
 * are aggregated by MODE (most frequent value), exactly matching
 * read_threshold_dac's `Counter(values).most_common(1)[0][0]`: this lets the
 * bulk of data-taking lines outvote occasional per-chip fine-tuning entries,
 * rather than just taking the first occurrence. Ties resolve to whichever
 * value was seen first (matching Python's stable-sort tie-breaking in
 * Counter.most_common).
 *
 * ACQWindow and DelayBetweenCycle appear once per file (in the
 * "== TriggerType: ... ==" line), so the first occurrence is used directly.
 */
#pragma once

#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace k4siwecal {

struct RunSettings {
  float acqWindowMs = -1.f;
  float delayBetweenCycleMs = -1.f;
  int thresholdDac = -1;
  int holdDelay = -1;
  int fsPeakTime = -1;
  int gainSelectionThreshold = -1;
};

namespace detail {

inline std::vector<int> findAllInts(const std::string& text, const std::regex& re) {
  std::vector<int> result;
  for (auto it = std::sregex_iterator(text.begin(), text.end(), re); it != std::sregex_iterator(); ++it) {
    result.push_back(std::stoi((*it)[1].str()));
  }
  return result;
}

inline bool findFirstFloat(const std::string& text, const std::regex& re, float& out) {
  std::smatch m;
  if (std::regex_search(text, m, re)) {
    out = std::stof(m[1].str());
    return true;
  }
  return false;
}

// Port of Python's Counter(values).most_common(1)[0][0]: most frequent value,
// ties broken by first-seen order (Python's sort is stable and Counter
// preserves insertion order, so most_common's reverse-count sort keeps
// equal-count elements in their original relative order).
inline int modeOfInts(const std::vector<int>& values) {
  if (values.empty()) return -1;
  std::vector<int> firstSeenOrder;
  std::map<int, int> counts;
  for (int v : values) {
    if (counts.find(v) == counts.end()) firstSeenOrder.push_back(v);
    ++counts[v];
  }
  int bestValue = firstSeenOrder.front();
  int bestCount = counts[bestValue];
  for (int v : firstSeenOrder) {
    if (counts[v] > bestCount) {
      bestCount = counts[v];
      bestValue = v;
    }
  }
  return bestValue;
}

}  // namespace detail

// Returns default-initialized (-1 sentinel) fields if the file cannot be
// opened, matching read_threshold_dac's "-1 on failure" convention.
inline RunSettings parseRunSettings(const std::string& path) {
  RunSettings settings;
  std::ifstream fin(path);
  if (!fin.is_open()) return settings;
  std::stringstream buffer;
  buffer << fin.rdbuf();
  const std::string text = buffer.str();

  static const std::regex kAcqWindowRe(R"(ACQWindow:\s*([0-9.]+))");
  static const std::regex kDelayBetweenCycleRe(R"(DelayBetweenCycle:\s*([0-9.]+))");
  static const std::regex kThresholdDacRe(R"(ThresholdDAC:\s*(\d+))");
  static const std::regex kHoldDelayRe(R"(HoldDelay:\s*(\d+))");
  static const std::regex kFsPeakTimeRe(R"(FSPeakTime:\s*(\d+))");
  static const std::regex kGainSelectionThresholdRe(R"(GainSelectionThreshold:\s*(\d+))");

  detail::findFirstFloat(text, kAcqWindowRe, settings.acqWindowMs);
  detail::findFirstFloat(text, kDelayBetweenCycleRe, settings.delayBetweenCycleMs);
  settings.thresholdDac = detail::modeOfInts(detail::findAllInts(text, kThresholdDacRe));
  settings.holdDelay = detail::modeOfInts(detail::findAllInts(text, kHoldDelayRe));
  settings.fsPeakTime = detail::modeOfInts(detail::findAllInts(text, kFsPeakTimeRe));
  settings.gainSelectionThreshold = detail::modeOfInts(detail::findAllInts(text, kGainSelectionThresholdRe));
  return settings;
}

}  // namespace k4siwecal
