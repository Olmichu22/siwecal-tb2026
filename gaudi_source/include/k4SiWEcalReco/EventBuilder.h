/*
 * EventBuilder: pure C++ (no ROOT/Gaudi dependency) port of
 * siwecal_eventbuilder's BCID clustering + hit collection + event assembly
 * (bcid_clustering.py::BcidClusterer, hit_collector.py::HitCollector,
 * event_builder.py::EventBuilder, config.py::BuilderConfig, models.py).
 *
 * This is the same "logic in a ROOT/Gaudi-light header, orchestration in the
 * .cpp" split already used for SlbFrameDecoder.h/EcalRawDecoder.cpp and
 * PedestalMipCalib.h/PedestalMipCalibrator.cpp -- kept testable with
 * synthetic AcquisitionView data.
 *
 * Two deliberate design notes carried over verbatim from the Python source:
 *
 * 1. The BCID overflow-adjustment map (buildOverflowMap, ported from
 *    bcid_clustering.py::_build_overflow_map) does NOT use the raw2root
 *    converter's `corrected_bcid` branch. That branch's overflow correction
 *    (SlbFrameDecoder.h::recordFrame) resets its `previousBcid`/`loopBcid`
 *    state PER FRAME (i.e. per slab/chip), not per whole acquisition cycle,
 *    which injects spurious +4096 offsets that would fragment real events.
 *    This class recomputes its own overflow tracking from the raw `bcid`
 *    branch instead (bcid_clustering.py:14-20).
 *
 * 2. In HitCollector::buildHit, the pad-map lookup uses the *geometric* chip
 *    loop index (0..15, matching how the map file is keyed -- see
 *    pad_map.py:18-20), while the calibration lookup and the Hit's chipId
 *    field (which becomes the output `hit_chip` branch) use the *hardware*
 *    chip_id read from the tree's `chipid[slab][chip]` branch
 *    (hit_collector.py:87-90, models.py:27). These are numerically equal in
 *    every currently-known SLB configuration (SlbFrameDecoder.h::recordFrame
 *    self-assigns `chipId[slab][chip] = chip`), but the two must be kept
 *    distinct rather than collapsed, matching the Python source exactly.
 */
#pragma once

#include "k4SiWEcalReco/CalibrationTables.h"
#include "k4SiWEcalReco/PadMapGeometry.h"
#include "k4SiWEcalReco/SlbFrameDecoder.h"  // kSlbDepth, kSkirocsPerAsu, kScasInSkiroc, kChannelsInSkiroc

#include <array>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>

namespace k4siwecal {

// Port of config.py::BuilderConfig (only the fields event building actually
// uses -- badbcid_max_good/default_workers are calibration-only/pipeline-only
// knobs not read by BcidClusterer/HitCollector, so they're omitted here).
struct BuilderConfig {
  int skipBcidStart = 50;
  std::set<long long> dropBcids = {0, 901};
  int mergeDelta = 3;
  int minSlabsHit = 10;
  long long bcidOverflow = 4096;
  long long badValue = -999;
  int adcUnderflowThreshold = 11;
  // Raw (NOT pedestal-subtracted) high-gain ADC at/above which the high-gain
  // preamp is taken to be saturated, so the hit energy is computed from the
  // low-gain branch instead.
  //
  // 1500, measured -- not the 1900 this used to be, and not the 1200 of
  // hit_collector.py. Using the LOW gain as the ruler (it is linear over this
  // whole range; the high gain is the one that saturates), slice in adc_low, take
  // the median adc_high per slice, and compare against the straight line the two
  // trace out well below any turn-over. The high gain falls short of that line by:
  //
  //     raw adc_high    1200     1500     1700     1900
  //     non-linearity   0.4%     1.2%     1.9%     3.9%      (run_000012 and
  //                                                            run_000072 agree)
  //
  // At the old 1900 the high gain was ALREADY 4% low, and every hit between 1500
  // and 1900 was read with that one-sided bias -- a systematic under-read, which
  // does not average out over the hits of an event the way noise does.
  //
  // 1900 was chosen to keep hits AWAY from the low gain, because the low gain was
  // calibrated through MIP_lg and MIP_lg is broken (see gainRatio below). That
  // reason is gone: the low gain is now anchored to the high gain and MIP_lg is
  // never used, so handing it a hit is cheap. What it costs is noise -- the
  // low-gain noise (~1 ADC) is amplified by 1/k = 10.4 -- and that trade crosses
  // over right around here:
  //
  //     switch          1200     1500     1900
  //     bias (syst.)    0.35%    1.15%    3.85%
  //     LG noise (rnd)  0.87%    0.69%    0.55%
  //
  // Measured effect of the move on run_000012: mean +0.16%, sigma/mu unchanged
  // (0.3499 -> 0.3500). Small -- but note that is itself the point. The anchored
  // energy barely depends on where the switch sits (+0.16% between 1500 and 1900)
  // whereas the old MIP_lg energy moved -2.68% over the same change. A correct
  // intercalibration is insensitive to the switch point; that insensitivity is the
  // evidence that it IS correct.
  int adcSaturationThreshold = 1500;

  // GAIN INTERCALIBRATION. This is how the hit energy is calibrated above the
  // saturation threshold; energyMipNoCalib keeps the old way, for comparison only.
  //
  // The old way calibrated the low-gain branch on its OWN MIP scale, and that scale
  // is not measurable: a MIP in the low gain is 2.08 ADC against ~1 ADC of noise
  // (S/N ~ 2), so the Landau fit lands on 3.2-4.5 ADC instead -- inflated 1.5x to
  // 2x -- and every run inflates it. Energy read through MIP_lg is wrong by that
  // factor.
  //
  // So do not calibrate the low gain against itself. Anchor it to the high gain,
  // which IS well calibrated, using the line the two branches trace out in the
  // FIDUCIAL REGION, where both are valid: adc_high - pedestal in [150, 1200].
  // That upper bound is not arbitrary -- it sits below the switch (no point fitting
  // where the high gain is no longer used) AND below the onset of high-gain
  // non-linearity, so the bend the switch exists to avoid cannot contaminate the
  // line meant to replace it.
  //
  //     adc_low - ped_lg  =  k * (adc_high - ped_hg)  +  c
  //
  // Inverting it turns a low-gain reading into the high-gain ADC the unsaturated
  // preamp WOULD have given, and the energy then comes from MIP_hg alone:
  //
  //     adc_high_equiv = (adc_low - ped_lg - c) / k
  //     energy         = adc_high_equiv / MIP_hg
  //
  // MIP_lg never appears. This is standard dual-gain intercalibration, not a fudge:
  // the low gain's job is to extend the dynamic range, not to be a second
  // independent calibration.
  //
  // The line is traced by the MEDIAN of adc_low in each slice of adc_high, never by
  // a least-squares fit to the hits: least squares on this band is dragged off the
  // data by its tails and reports c = +3.13 where the medians say +1.60.
  //
  // k = 0.0962 (1/k = 10.4) is solid: run_000012 (th230) gives 0.0961 and
  // run_000072 (th220) 0.0963 -- two runs, two thresholds, agreeing to 0.2%.
  //
  // c = +1.45 ADC is real (always positive, seen by two independent robust methods:
  // this line, and the low-gain baseline measured directly on hits that fired with
  // no high-gain signal) but only good to ~0.5 ADC. The error is SYSTEMATIC -- the
  // value moves with the fit region -- not the +/-0.03 a least-squares covariance
  // would claim. It matters little: c shifts the energy by c/k ~ 15 high-gain ADC
  // ~ 0.7 MIP, against a switch point of 1500 ADC ~ 58 MIP, so ~1%.
  double gainRatio = 0.0962;     // k: adc_low per unit adc_high
  double gainIntercept = 1.45;   // c: low-gain ADC at zero high-gain signal
  double maxHitsPerSca = 1.0e18;  // math.inf equivalent; see EcalEventBuilder.cpp for why not literal infinity
  int maxHitsPerEvent = 15360;  // 15 slabs * 16 chips * 64 channels
};

using ChipKey = std::pair<int, int>;  // (slab, chip)

// Port of models.py::Hit.
struct Hit {
  int slabPosition = 0;  // -> hit_slab
  int chipId = 0;        // hardware chip ID -> hit_chip
  int channel = 0;       // -> hit_chan
  int sca = 0;            // -> hit_sca
  float adcHighPedsub = 0.f;  // -> hit_hg
  float adcLowPedsub = 0.f;   // -> hit_lg
  // The hit energy, in MIP units. Above the saturation threshold it comes from the
  // low gain anchored to the high gain (see BuilderConfig::gainRatio). Everything
  // downstream -- sumEnergy, wEnergy, sumWEnergy, the shower variables -- is built
  // from THIS one.
  float energyMip = 0.f;      // -> hit_energy
  // The same hit, calibrated the old way: above the threshold, the low gain on its
  // own MIP_lg scale. Kept ONLY so the two can be compared on the same events; it
  // is not used to build anything. MIP_lg is inflated 1.5x-2x (see gainRatio), so
  // this systematically under-reads the energy of every saturated hit.
  float energyMipNoCalib = 0.f;  // -> hit_energy_nocalib
  float wEnergy = 0.f;        // -> hit_w_energy (sampling-corrected: energyMip * W[slab]/X0)
  float x = std::numeric_limits<float>::quiet_NaN();   // -> hit_x
  float y = std::numeric_limits<float>::quiet_NaN();   // -> hit_y
  float z = std::numeric_limits<float>::quiet_NaN();   // -> hit_z
  float x0 = std::numeric_limits<float>::quiet_NaN();  // -> hit_X0
  bool isMasked = false;                                // -> hit_ismasked
};

// Port of models.py::BcidWindow.
struct BcidWindow {
  long long bcidLabel = 0;
  long long startRaw = 0;
  long long stopRaw = 0;
  std::map<ChipKey, std::vector<int>> chipToScas;  // active SCAs per (slab,chip)
};

// Port of models.py::ReconstructedEvent.
struct ReconstructedEvent {
  long long bcid = 0;
  std::vector<Hit> hits;

  int nChannels() const { return static_cast<int>(hits.size()); }
  int nSlabs() const {
    std::set<int> slabs;
    for (const auto& h : hits) slabs.insert(h.slabPosition);
    return static_cast<int>(slabs.size());
  }
  int nChips() const {
    std::set<std::pair<int, int>> chips;
    for (const auto& h : hits) chips.insert({h.slabPosition, h.chipId});
    return static_cast<int>(chips.size());
  }
  double sumAdcHigh() const {
    double s = 0.0;
    for (const auto& h : hits) s += h.adcHighPedsub;
    return s;
  }
  double sumEnergy() const {
    double s = 0.0;
    for (const auto& h : hits) s += h.energyMip;
    return s;
  }
  // The event energy the OLD way, for the side-by-side comparison. Not used to
  // build anything.
  double sumEnergyNoCalib() const {
    double s = 0.0;
    for (const auto& h : hits) s += h.energyMipNoCalib;
    return s;
  }
  double sumWEnergy() const {
    double s = 0.0;
    for (const auto& h : hits) s += h.wEnergy;
    return s;
  }
};

// Non-owning read-only view over one already-loaded siwecaldecoded tree
// entry's bound branch buffers, matching root_io.py::Acquisition's accessor
// shape (but binding fixed C arrays directly, like
// PedestalMipCalibrator.cpp::TreeBuffers, rather than Python's per-leaf
// TLeaf::GetValue -- simpler and faster in C++). Only the branches actually
// used by BcidClusterer/HitCollector are bound (matches
// AcquisitionReader._LEAF_NAMES minus `badbcid`, which is read in Python but
// never referenced by the clustering/collection algorithms themselves).
class AcquisitionView {
 public:
  AcquisitionView(int nSlboards, const int (&slboardId)[kSlbDepth], const int (&chipId)[kSlbDepth][kSkirocsPerAsu],
                  const int (&bcid)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc],
                  const int (&nhits)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc],
                  const int (&hitbitHigh)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc],
                  const int (&adcHigh)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc],
                  const int (&adcLow)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc])
      : m_nSlboards(nSlboards),
        m_slboardId(slboardId),
        m_chipId(chipId),
        m_bcid(bcid),
        m_nhits(nhits),
        m_hitbitHigh(hitbitHigh),
        m_adcHigh(adcHigh),
        m_adcLow(adcLow) {}

  int nSlabPositions() const { return m_nSlboards; }
  int slboardId(int slab) const { return m_slboardId[slab]; }
  int chipId(int slab, int chip) const { return m_chipId[slab][chip]; }
  int nHits(int slab, int chip, int sca) const { return m_nhits[slab][chip][sca]; }
  long long rawBcid(int slab, int chip, int sca) const { return m_bcid[slab][chip][sca]; }
  int hitbitHigh(int slab, int chip, int sca, int channel) const { return m_hitbitHigh[slab][chip][sca][channel]; }
  int adcHigh(int slab, int chip, int sca, int channel) const { return m_adcHigh[slab][chip][sca][channel]; }
  int adcLow(int slab, int chip, int sca, int channel) const { return m_adcLow[slab][chip][sca][channel]; }

  // Port of root_io.py::Acquisition.raw_bcid_matrix: flattens to
  // (n_chip_rows, n_scas), row = slab*n_chips_per_slab + chip.
  using BcidMatrix = std::array<std::array<long long, kScasInSkiroc>, kSlbDepth * kSkirocsPerAsu>;
  BcidMatrix rawBcidMatrix() const {
    BcidMatrix m{};
    for (int slab = 0; slab < kSlbDepth; ++slab) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        const int row = slab * kSkirocsPerAsu + chip;
        for (int sca = 0; sca < kScasInSkiroc; ++sca) {
          m[row][sca] = m_bcid[slab][chip][sca];
        }
      }
    }
    return m;
  }

 private:
  int m_nSlboards;
  const int (&m_slboardId)[kSlbDepth];
  const int (&m_chipId)[kSlbDepth][kSkirocsPerAsu];
  const int (&m_bcid)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
  const int (&m_nhits)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc];
  const int (&m_hitbitHigh)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  const int (&m_adcHigh)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
  const int (&m_adcLow)[kSlbDepth][kSkirocsPerAsu][kScasInSkiroc][kChannelsInSkiroc];
};

// Port of bcid_clustering.py::BcidClusterer.
class BcidClusterer {
 public:
  explicit BcidClusterer(const BuilderConfig& cfg) : m_cfg(cfg) {}

  // Port of find_windows (bcid_clustering.py:43-71).
  std::vector<BcidWindow> findWindows(const AcquisitionView& acq) const {
    const auto matrix = acq.rawBcidMatrix();
    const auto overflowMap = buildOverflowMap(matrix);
    const auto chipToScasBcid = collectValidBcids(acq, matrix);

    std::vector<BcidWindow> windows;
    if (chipToScasBcid.empty()) return windows;

    for (const auto& [startRaw, stopRaw] : mergeIntoWindows(chipToScasBcid)) {
      auto active = chipsInWindow(chipToScasBcid, startRaw, stopRaw);
      std::set<int> slabs;
      for (const auto& [key, scas] : active) slabs.insert(key.first);
      if (static_cast<int>(slabs.size()) < m_cfg.minSlabsHit) continue;

      BcidWindow w;
      w.bcidLabel = overflowCorrectedLabel(startRaw, stopRaw, overflowMap);
      w.startRaw = startRaw;
      w.stopRaw = stopRaw;
      w.chipToScas = std::move(active);
      windows.push_back(std::move(w));
    }
    return windows;
  }

 private:
  using BcidMatrix = AcquisitionView::BcidMatrix;

  // Port of _build_overflow_map (bcid_clustering.py:75-108). See the file
  // header's note 1 for why this deliberately avoids `corrected_bcid`.
  std::map<long long, long long> buildOverflowMap(const BcidMatrix& matrix) const {
    constexpr int kRows = kSlbDepth * kSkirocsPerAsu;
    constexpr int kSteps = kScasInSkiroc - 1;  // 14 sca-to-sca transitions

    // n_cycles[row][k]: cumulative count of BCID decreases up to and
    // including transition k (sca k -> sca k+1).
    std::array<std::array<long long, kSteps>, kRows> nCycles{};
    std::array<std::array<bool, kSteps>, kRows> startsNewCycle{};

    for (int r = 0; r < kRows; ++r) {
      long long cum = 0;
      for (int k = 0; k < kSteps; ++k) {
        const bool decreased = (matrix[r][k + 1] - matrix[r][k]) < 0;
        cum += decreased ? 1 : 0;
        nCycles[r][k] = cum;
        // The very first step (k==0) is never counted as a cycle start
        // (matches starts_new_cycle[:, 1:] leaving column 0 at its default
        // False -- Python's starts_new_cycle array has the same width as
        // n_cycles, and only indices 1.. are ever set).
        startsNewCycle[r][k] = (k == 0) ? false : (nCycles[r][k] != nCycles[r][k - 1]);
        if (matrix[r][k + 1] == m_cfg.badValue) startsNewCycle[r][k] = false;
      }
    }

    int totalStarts = 0;
    for (const auto& row : startsNewCycle)
      for (bool b : row) totalStarts += b ? 1 : 0;

    std::map<long long, long long> overflowMap;
    if (totalStarts > 1) {
      std::set<long long> adjustedValues;
      for (int r = 0; r < kRows; ++r) {
        for (int k = 0; k < kSteps; ++k) {
          if (startsNewCycle[r][k]) {
            adjustedValues.insert(matrix[r][k + 1] + m_cfg.bcidOverflow * nCycles[r][k]);
          }
        }
      }
      for (long long adjusted : adjustedValues) {
        const long long key = ((adjusted % m_cfg.bcidOverflow) + m_cfg.bcidOverflow) % m_cfg.bcidOverflow;
        overflowMap[key] = adjusted;
      }
    }
    return overflowMap;
  }

  // Port of _collect_valid_bcids (bcid_clustering.py:110-138). Iterates all
  // SCAs regardless of nColumns, matching the Python's own comment.
  std::map<ChipKey, std::map<int, long long>> collectValidBcids(const AcquisitionView& acq,
                                                                  const BcidMatrix& matrix) const {
    std::map<ChipKey, std::map<int, long long>> result;
    for (int slab = 0; slab < acq.nSlabPositions(); ++slab) {
      for (int chip = 0; chip < kSkirocsPerAsu; ++chip) {
        if (acq.chipId(slab, chip) < 0) continue;
        for (int sca = 0; sca < kScasInSkiroc; ++sca) {
          const int nHits = acq.nHits(slab, chip, sca);
          if (nHits == 0 || static_cast<double>(nHits) > m_cfg.maxHitsPerSca) continue;
          const long long rawBcid = matrix[slab * kSkirocsPerAsu + chip][sca];
          if (rawBcid < 0 || rawBcid < m_cfg.skipBcidStart || m_cfg.dropBcids.count(rawBcid)) continue;
          result[{slab, chip}][sca] = rawBcid;
        }
      }
    }
    return result;
  }

  // Port of _merge_into_windows (bcid_clustering.py:140-164).
  std::vector<std::pair<long long, long long>> mergeIntoWindows(
      const std::map<ChipKey, std::map<int, long long>>& chipToScasBcid) const {
    std::set<long long> uniqueBcids;
    for (const auto& [key, scas] : chipToScasBcid)
      for (const auto& [sca, bcid] : scas) uniqueBcids.insert(bcid);

    std::vector<std::pair<long long, long long>> windows;
    if (uniqueBcids.empty()) return windows;

    auto it = uniqueBcids.begin();
    long long start = *it;
    long long stop = *it;
    for (++it; it != uniqueBcids.end(); ++it) {
      const long long bcid = *it;
      if (bcid - stop < m_cfg.mergeDelta) {
        stop = bcid;
      } else {
        windows.emplace_back(start, stop);
        start = stop = bcid;
      }
    }
    windows.emplace_back(start, stop);
    return windows;
  }

  // Port of _chips_in_window (bcid_clustering.py:166-175).
  static std::map<ChipKey, std::vector<int>> chipsInWindow(
      const std::map<ChipKey, std::map<int, long long>>& chipToScasBcid, long long startRaw, long long stopRaw) {
    std::map<ChipKey, std::vector<int>> active;
    for (const auto& [key, scas] : chipToScasBcid) {
      std::vector<int> matched;
      for (const auto& [sca, bcid] : scas) {
        if (bcid >= startRaw && bcid <= stopRaw) matched.push_back(sca);
      }
      if (!matched.empty()) active[key] = std::move(matched);
    }
    return active;
  }

  // Port of _overflow_corrected_label (bcid_clustering.py:177-188). Integer
  // division below is floor division for non-negative operands (the only
  // case that arises here: BCID values and bcidOverflow are non-negative),
  // matching Python's `//`.
  long long overflowCorrectedLabel(long long startRaw, long long stopRaw,
                                    const std::map<long long, long long>& overflowMap) const {
    long long unwrappedMax = std::numeric_limits<long long>::min();
    for (long long bcid = startRaw; bcid <= stopRaw; ++bcid) {
      const auto it = overflowMap.find(bcid);
      const long long val = (it != overflowMap.end()) ? it->second : bcid;
      if (val > unwrappedMax) unwrappedMax = val;
    }
    const long long nCycles = unwrappedMax / m_cfg.bcidOverflow;
    return startRaw + nCycles * m_cfg.bcidOverflow;
  }

  const BuilderConfig& m_cfg;
};

// Port of hit_collector.py::HitCollector.
class HitCollector {
 public:
  HitCollector(const BuilderConfig& cfg, const CalibrationTables& calib, const SlabGeometry& geom,
               const PadMap* padMap)
      : m_cfg(cfg), m_calib(calib), m_geom(geom), m_padMap(padMap) {}

  // Port of collect (hit_collector.py:33-45).
  std::vector<Hit> collect(const AcquisitionView& acq, const BcidWindow& window) const {
    std::vector<Hit> hits;
    for (const auto& [key, scas] : window.chipToScas) {
      const int slab = key.first;
      const int chip = key.second;  // geometric loop index -- used for pad-map lookup below
      const int chipId = acq.chipId(slab, chip);  // hardware ID -- used for calibration + output branch
      if (chipId < 0) continue;
      const int slabId = acq.slboardId(slab);
      for (const auto& [channel, sca] : bestScaPerChannel(acq, slab, chip, scas)) {
        auto hit = buildHit(acq, slab, chip, sca, channel, slabId, chipId);
        if (hit.has_value()) hits.push_back(*hit);
      }
    }
    return hits;
  }

 private:
  // Port of _best_sca_per_channel (hit_collector.py:49-65): for each
  // triggered channel, keep the SCA with the largest hitbit_high (strictly
  // greater, so the first-seen SCA wins ties).
  std::map<int, int> bestScaPerChannel(const AcquisitionView& acq, int slab, int chip,
                                        const std::vector<int>& scas) const {
    std::map<int, int> bestSca;
    std::map<int, int> bestResponse;
    for (int sca : scas) {
      for (int channel = 0; channel < kChannelsInSkiroc; ++channel) {
        const int response = acq.hitbitHigh(slab, chip, sca, channel);
        if (response <= 0) continue;
        const auto it = bestResponse.find(channel);
        if (it == bestResponse.end() || response > it->second) {
          bestResponse[channel] = response;
          bestSca[channel] = sca;
        }
      }
    }
    return bestSca;
  }

  // Port of _build_hit (hit_collector.py:67-106).
  std::optional<Hit> buildHit(const AcquisitionView& acq, int slab, int chip, int sca, int channel, int slabId,
                               int chipId) const {
    const int adcHigh = acq.adcHigh(slab, chip, sca, channel);
    if (adcHigh <= m_cfg.adcUnderflowThreshold) return std::nullopt;
    const int adcLow = acq.adcLow(slab, chip, sca, channel);

    const double pedestalHg = m_calib.pedestal(slabId, chipId, channel, sca);
    const double mipHg = m_calib.mip(slabId, chipId, channel);
    // Low-gain pedestal, NOT the high-gain one: the two branches have their own
    // baselines, so subtracting the high-gain pedestal from adc_low (as this
    // did before low-gain support) yields a meaningless hit_lg. Without low-gain
    // tables there is nothing better to subtract, so that old behaviour is kept
    // rather than silently turning hit_lg into an unsubtracted raw ADC.
    const double pedestalLg =
        m_calib.hasLowGain() ? m_calib.pedestalLg(slabId, chipId, channel, sca) : pedestalHg;
    const double mipLg = m_calib.mipLg(slabId, chipId, channel);

    const double adcHighPedsub = adcHigh - pedestalHg;
    const double adcLowPedsub = adcLow - pedestalLg;

    // The high-gain preamp saturates on large deposits, flattening adc_high and so
    // under-reading the energy; above the saturation threshold the low-gain branch
    // (which is still linear there) is used instead. The switch is decided on the
    // PEDESTAL-SUBTRACTED high-gain value (adc_high - ped_hg) -- the ADC quantity we
    // use everywhere else and are familiar with at user level -- not the raw ADC.
    const bool isMasked = m_calib.isMasked(slabId, chipId, channel);
    const bool saturated = m_calib.hasLowGain() && adcHighPedsub >= m_cfg.adcSaturationThreshold;

    // The energy, and the same energy the old way. They differ ONLY for saturated
    // hits: below the threshold both read the high gain and are identical.
    double energyMip = 0.0;
    double energyMipNoCalib = 0.0;
    const double hgEnergy = (mipHg > 0) ? (adcHighPedsub / mipHg) : 0.0;

    if (isMasked) {
      energyMip = 0.0;
      energyMipNoCalib = 0.0;
    } else if (!saturated) {
      energyMip = hgEnergy;
      energyMipNoCalib = hgEnergy;
    } else {
      // DEFAULT: anchor the low gain to the high gain. Invert the line the two
      // branches trace out, adc_low - ped_lg = k*(adc_high - ped_hg) + c, to get
      // the high-gain ADC an unsaturated preamp would have given, then divide by
      // MIP_hg. MIP_lg -- the number we cannot measure -- never enters.
      //
      // Only the low-gain PEDESTAL has to be good for this, not the low-gain MIP,
      // so it gates on isPedestalMaskedLg and not on isMaskedLg (which also masks
      // every channel whose MIP_lg fit failed -- and those fits fail precisely
      // because MIP_lg is unmeasurable, which is the whole reason we are here).
      if (mipHg > 0 && m_cfg.gainRatio > 0 && !m_calib.isPedestalMaskedLg(slabId, chipId, channel)) {
        const double adcHighEquiv = (adcLowPedsub - m_cfg.gainIntercept) / m_cfg.gainRatio;
        energyMip = adcHighEquiv / mipHg;
      }
      // The old way, for comparison: the low gain on its own inflated MIP scale.
      if (!m_calib.isMaskedLg(slabId, chipId, channel) && mipLg > 0) {
        energyMipNoCalib = adcLowPedsub / mipLg;
      }
    }
    // A saturated hit with no usable low-gain pedestal gets no energy at all (0),
    // rather than being silently under-read from the flattened high-gain ADC.

    float x = std::numeric_limits<float>::quiet_NaN();
    float y = std::numeric_limits<float>::quiet_NaN();
    if (m_padMap != nullptr) {
      // Geometric chip index (`chip`), NOT chipId -- see file header note 2.
      std::tie(x, y) = m_padMap->position(slab, chip, channel);
    }

    Hit hit;
    hit.slabPosition = slab;
    hit.chipId = chipId;
    hit.channel = channel;
    hit.sca = sca;
    hit.adcHighPedsub = static_cast<float>(adcHighPedsub);
    hit.adcLowPedsub = static_cast<float>(adcLowPedsub);
    hit.energyMip = static_cast<float>(energyMip);
    hit.energyMipNoCalib = static_cast<float>(energyMipNoCalib);
    hit.x = x;
    hit.y = y;
    hit.z = static_cast<float>(m_geom.slabZ(slab));
    hit.x0 = static_cast<float>(m_geom.slabX0(slab));
    // Sampling correction: layers do not all sit behind the same amount of
    // tungsten (2.8mm / 4.2mm / 5.6mm here), so a MIP deposited in a thick-
    // absorber layer stands for more incident energy than one in a thin layer.
    // Same weight EcalShowerVars.h::hitWeights applies when it forms `weighte`.
    hit.wEnergy = static_cast<float>(energyMip * m_geom.slabWOverX0(slab));
    hit.isMasked = isMasked;
    return hit;
  }

  const BuilderConfig& m_cfg;
  const CalibrationTables& m_calib;
  const SlabGeometry& m_geom;
  const PadMap* m_padMap;  // nullptr => hit_x/hit_y stay NaN (no mapping)
};

// Port of event_builder.py::EventBuilder.
class EventBuilder {
 public:
  EventBuilder(const BuilderConfig& cfg, const CalibrationTables& calib, const SlabGeometry& geom,
               const PadMap* padMap)
      : m_clusterer(cfg), m_hitCollector(cfg, calib, geom, padMap) {}

  // Port of build (event_builder.py:28-41). Windows whose hit collection
  // yields no surviving hits are dropped -- a window can have valid BCIDs
  // yet no channel passing the underflow cut.
  std::vector<ReconstructedEvent> build(const AcquisitionView& acq) const {
    std::vector<ReconstructedEvent> events;
    for (const auto& window : m_clusterer.findWindows(acq)) {
      auto hits = m_hitCollector.collect(acq, window);
      if (hits.empty()) continue;
      ReconstructedEvent ev;
      ev.bcid = window.bcidLabel;
      ev.hits = std::move(hits);
      events.push_back(std::move(ev));
    }
    return events;
  }

 private:
  BcidClusterer m_clusterer;
  HitCollector m_hitCollector;
};

}  // namespace k4siwecal
