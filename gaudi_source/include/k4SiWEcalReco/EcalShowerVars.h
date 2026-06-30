/*
 * Per-event particle-discrimination variables for the SiW-ECAL.
 *
 * C++ port of siwecal_validation/metrics.py (itself a port of
 * VarsJesus/analysis_template.C). metrics.py is the parity oracle, so where the
 * original .C and metrics.py differ (e.g. the energy>0 filter on transverse
 * moments, the shower thresholds), we follow metrics.py.
 *
 * Everything works on plain std::vector for one event, with no ROOT/Gaudi
 * dependency, so it is reusable and unit-testable. The canonical ordering of
 * the variables (scalarNames / perLayerNames) is the single source of truth for
 * the Cluster::shapeParameters layout written by EcalPidTransformer and read
 * back downstream.
 */
#ifndef K4SIWECALRECO_ECALSHOWERVARS_H
#define K4SIWECALRECO_ECALSHOWERVARS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace k4siwecal {

inline constexpr float NANF = std::numeric_limits<float>::quiet_NaN();

struct ShowerThresholds {
  float e_threshold = 5.0f;   ///< min per-layer activity to count as shower material
  float max_min = 10.0f;      ///< min profile peak for the event to qualify as shower
  float start_frac = 0.1f;    ///< fraction of the peak for the *_10 start/end variants
};

/// One event's worth of derived variables (mirrors EventData fields).
struct EventVars {
  // scalars
  float nhit = 0.f;
  float zbary = NANF;
  float energy = 0.f;
  float mip_likeness = 0.f;
  float weighte = 0.f;
  float bar_x = NANF;
  float bar_y = NANF;
  float bar_r = NANF;
  float moliere = NANF;
  float transverse_rms = NANF;
  float is_shower = 0.f;        // 0/1
  float shower_start = NANF;
  float shower_max = NANF;
  float shower_end = NANF;
  float shower_start_10 = NANF;
  float shower_end_10 = NANF;
  float shower_length = NANF;
  float first_layer = NANF;
  float last_layer = NANF;
  float n_layers_hit = 0.f;
  float e_over_nhit = NANF;
  // per-layer profiles (length n_layers)
  std::vector<float> hits_per_layer;
  std::vector<float> energy_per_layer;
  std::vector<float> weighte_per_layer;
};

/// Canonical, ordered names of the scalar variables (one float each).
inline const std::vector<std::string>& scalarNames() {
  static const std::vector<std::string> names = {
      "nhit", "zbary", "energy", "mip_likeness", "weighte", "bar_x", "bar_y", "bar_r",
      "moliere", "transverse_rms", "is_shower", "shower_start", "shower_max", "shower_end",
      "shower_start_10", "shower_end_10", "shower_length", "first_layer", "last_layer",
      "n_layers_hit", "e_over_nhit"};
  return names;
}

/// Canonical names of the per-layer blocks (each expands to [n_layers] floats).
inline const std::vector<std::string>& perLayerNames() {
  static const std::vector<std::string> names = {
      "hits_per_layer", "energy_per_layer", "weighte_per_layer"};
  return names;
}

// --------------------------------------------------------------------------
namespace detail {

// Accumulate in double (like numpy's bincount) and store the per-layer profile
// as float. An empty weights vector means "count hits".
template <typename W>
inline std::vector<float> bincount(const std::vector<int>& slab,
                                   const std::vector<W>& weights, int n_layers) {
  std::vector<double> acc(n_layers, 0.0);
  for (std::size_t i = 0; i < slab.size(); ++i) {
    const int s = slab[i];
    if (s >= 0 && s < n_layers) acc[s] += weights.empty() ? 1.0 : static_cast<double>(weights[i]);
  }
  return std::vector<float>(acc.begin(), acc.end());
}

}  // namespace detail

/// Tungsten-weighted per-hit energy: E_i * (W[slab_i] / X0). Kept in double to
/// match metrics.py (float64), which matters for the Molière cumulant boundary.
inline std::vector<double> hitWeights(const std::vector<int>& slab,
                                      const std::vector<float>& energy,
                                      const std::vector<double>& w_over_x0) {
  std::vector<double> w(energy.size(), 0.0);
  for (std::size_t i = 0; i < energy.size(); ++i) {
    const int s = slab[i];
    const double wx0 = (s >= 0 && s < static_cast<int>(w_over_x0.size())) ? w_over_x0[s] : 0.0;
    w[i] = static_cast<double>(energy[i]) * wx0;
  }
  return w;
}

/// _is_shower: a rising EM-like edge (3 consecutive increasing layers above
/// threshold) with a peak above max_min.
inline bool isShower(const std::vector<float>& profile, const ShowerThresholds& thr) {
  if (profile.empty()) return false;
  const float peak = *std::max_element(profile.begin(), profile.end());
  if (peak <= thr.max_min) return false;
  for (std::size_t i = 0; i + 2 < profile.size(); ++i) {
    if (profile[i] > thr.e_threshold && profile[i + 1] > profile[i] &&
        profile[i + 2] > profile[i + 1])
      return true;
  }
  return false;
}

/// Compute every derived variable for one hit set (already energy-cut if needed).
/// Mirrors siwecal_validation/event_data.py's per-event loop.
inline EventVars computeEventVars(const std::vector<int>& slab,
                                  const std::vector<float>& energy,
                                  const std::vector<float>& x,
                                  const std::vector<float>& y,
                                  const std::vector<double>& w_over_x0, int n_layers,
                                  const std::string& profile_kind,
                                  const ShowerThresholds& thr, double containment) {
  EventVars v;
  v.hits_per_layer.assign(n_layers, 0.f);
  v.energy_per_layer.assign(n_layers, 0.f);
  v.weighte_per_layer.assign(n_layers, 0.f);

  const std::size_t n = slab.size();
  if (n == 0) return v;

  // Longitudinal sums use signed energy over all hits.
  const std::vector<double> w = hitWeights(slab, energy, w_over_x0);
  v.hits_per_layer = detail::bincount(slab, std::vector<float>{}, n_layers);
  v.energy_per_layer = detail::bincount(slab, energy, n_layers);
  v.weighte_per_layer = detail::bincount(slab, w, n_layers);

  // Accumulate in double to match numpy (metrics.py) float64 sums.
  double esum = 0., wsum = 0., e_dot_slab = 0.;
  for (std::size_t i = 0; i < n; ++i) {
    esum += energy[i];
    wsum += w[i];
    e_dot_slab += static_cast<double>(energy[i]) * slab[i];
  }
  v.nhit = static_cast<float>(n);
  v.energy = static_cast<float>(esum);
  v.weighte = static_cast<float>(wsum);
  v.zbary = (esum > 0.) ? static_cast<float>(e_dot_slab / esum) : NANF;
  v.e_over_nhit = (n > 0) ? static_cast<float>(esum / static_cast<double>(n)) : NANF;

  // MIP-likeness from the per-layer hit counts.
  float score = 0.f;
  for (int i = 0; i < n_layers; ++i)
    if (v.hits_per_layer[i] > 0.f) score += 1.f / v.hits_per_layer[i];
  v.mip_likeness = score / static_cast<float>(n_layers);

  // Layer extent.
  int first = -1, last = -1, nhit_layers = 0;
  for (int i = 0; i < n_layers; ++i)
    if (v.hits_per_layer[i] > 0.f) {
      if (first < 0) first = i;
      last = i;
      ++nhit_layers;
    }
  v.first_layer = (first >= 0) ? static_cast<float>(first) : NANF;
  v.last_layer = (last >= 0) ? static_cast<float>(last) : NANF;
  v.n_layers_hit = static_cast<float>(nhit_layers);

  // Transverse moments use only hits with positive energy (pedestal-subtracted
  // noise has negative energy and would bias the centroid / Molière cumulant).
  std::vector<float> px, py;
  std::vector<double> pw;
  px.reserve(n);
  py.reserve(n);
  pw.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    if (energy[i] > 0.f) {
      px.push_back(x[i]);
      py.push_back(y[i]);
      pw.push_back(w[i]);
    }
  // Keep the barycenter in double for the internal r computations (Molière /
  // RMS), matching metrics.py; only the stored bar_* fields are cast to float.
  double ptot = std::accumulate(pw.begin(), pw.end(), 0.0);
  double bdx = 0., bdy = 0.;
  if (ptot > 0.) {
    double sx = 0., sy = 0.;
    for (std::size_t i = 0; i < px.size(); ++i) {
      sx += pw[i] * px[i];
      sy += pw[i] * py[i];
    }
    bdx = sx / ptot;
    bdy = sy / ptot;
  }
  v.bar_x = static_cast<float>(bdx);
  v.bar_y = static_cast<float>(bdy);
  v.bar_r = static_cast<float>(std::hypot(bdx, bdy));

  // Transverse RMS about the shower axis.
  if (ptot > 0.) {
    double acc = 0.;
    for (std::size_t i = 0; i < px.size(); ++i) {
      const double r2 = (px[i] - bdx) * (px[i] - bdx) + (py[i] - bdy) * (py[i] - bdy);
      acc += pw[i] * r2;
    }
    v.transverse_rms = static_cast<float>(std::sqrt(std::max(0.0, acc / ptot)));
  } else {
    v.transverse_rms = 0.f;
  }

  // Longitudinal shower descriptors from the configured profile.
  const std::vector<float>& profile =
      (profile_kind == "sume")      ? v.energy_per_layer
      : (profile_kind == "weighte") ? v.weighte_per_layer
                                    : v.hits_per_layer;
  const bool shower = isShower(profile, thr);
  v.is_shower = shower ? 1.f : 0.f;
  if (shower) {
    const float peak = *std::max_element(profile.begin(), profile.end());
    int max_layer = static_cast<int>(
        std::max_element(profile.begin(), profile.end()) - profile.begin());
    v.shower_max = static_cast<float>(max_layer);
    int length = 0;
    for (int i = 0; i < n_layers; ++i)
      if (profile[i] > thr.e_threshold) ++length;
    v.shower_length = static_cast<float>(length);

    auto firstAbove = [&](bool ten, int lo, int hi) -> float {
      for (int i = lo; i < hi; ++i)
        if (profile[i] > thr.e_threshold && (!ten || profile[i] > thr.start_frac * peak))
          return static_cast<float>(i);
      return NANF;
    };
    auto lastAbove = [&](bool ten, int lo, int hi) -> float {
      for (int i = hi - 1; i >= lo; --i)
        if (profile[i] > thr.e_threshold && (!ten || profile[i] > thr.start_frac * peak))
          return static_cast<float>(i);
      return NANF;
    };
    v.shower_start = firstAbove(false, 0, max_layer);
    v.shower_end = lastAbove(false, max_layer + 1, n_layers);
    v.shower_start_10 = firstAbove(true, 0, max_layer);
    v.shower_end_10 = lastAbove(true, max_layer + 1, n_layers);
  }

  // Molière radius (90% containment) only for showers.
  if (shower && ptot > 0.f && !px.empty()) {
    std::vector<std::size_t> order(px.size());
    std::iota(order.begin(), order.end(), 0);
    std::vector<double> r(px.size());  // double to match numpy's float64 hypot
    for (std::size_t i = 0; i < px.size(); ++i)
      r[i] = std::hypot(static_cast<double>(px[i]) - bdx, static_cast<double>(py[i]) - bdy);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return r[a] < r[b]; });
    const double target = containment * ptot;
    double cum = 0.;
    double mol = r[order.back()];
    for (std::size_t k = 0; k < order.size(); ++k) {
      cum += pw[order[k]];
      if (cum >= target) {
        mol = r[order[k]];
        break;
      }
    }
    v.moliere = static_cast<float>(mol);
  } else {
    v.moliere = shower ? 0.f : 0.f;  // metrics.py: 0.0 when not a shower
  }

  return v;
}

}  // namespace k4siwecal

#endif  // K4SIWECALRECO_ECALSHOWERVARS_H
