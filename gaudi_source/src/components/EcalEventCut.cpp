/*
 * EcalEventCut: keep only the events whose value of ONE shower variable falls
 * inside [Min, Max]. A FilterPredicate, so an event it rejects is never handed
 * to the writer and never reaches the output file.
 *
 * One instance = one cut. Chain instances in TopAlg to AND several cuts
 * together; each is independent and none of them rewrites the event, so the
 * surviving events keep every collection they arrived with, byte for byte.
 *
 *   cut1 = EcalEventCut("cut_moliere", Variable="moliere", Min=45)
 *   cut2 = EcalEventCut("cut_weighte", Variable="weighte", Min=5000)
 *   ApplicationMgr(TopAlg=[source, pid, cut1, cut2], ...)
 *
 * Variable is looked up in EcalShowerVars.h's scalarNames() -- the SAME list
 * EcalPidTransformer flattens into Cluster::shapeParameters, so the index is
 * derived from the one canonical order rather than from a second copy of it
 * living here. Only the base scalars (the first block) are addressable by name;
 * for anything further into the layout (per-layer profiles, MIP-variant
 * blocks), set ParamIndex directly.
 *
 * Bounds are INCLUSIVE and default to +-infinity, so a one-sided cut is just
 * Min (or just Max). An event whose value is NaN never passes: NaN compares
 * false against everything, which is the safe direction for a selection.
 */
#include "k4FWCore/FilterPredicate.h"

#include "edm4hep/ClusterCollection.h"

#include "Gaudi/Property.h"

#include "k4SiWEcalReco/EcalShowerVars.h"

#include <algorithm>
#include <atomic>
#include <limits>
#include <string>

struct EcalEventCut final : k4FWCore::FilterPredicate<bool(const edm4hep::ClusterCollection&)> {
  EcalEventCut(const std::string& name, ISvcLocator* svcLoc)
      // KeyValue, not KeyValues: k4FWCore aliases Gaudi's FilterPredicate
      // straight through, so it keeps the single-key signature (unlike
      // Transformer, which k4FWCore wraps to take a vector of keys).
      : FilterPredicate(name, svcLoc, KeyValue("InputClusters", "ECalPid")) {}

  StatusCode initialize() override {
    if (m_paramIndex.value() >= 0) {
      m_index = m_paramIndex.value();
    } else {
      const auto& names = k4siwecal::scalarNames();
      const auto it = std::find(names.begin(), names.end(), m_variable.value());
      if (it == names.end()) {
        error() << "Variable '" << m_variable.value() << "' is not one of the base scalars. "
                << "Use one of: ";
        for (const auto& n : names) error() << n << " ";
        error() << "-- or set ParamIndex to address the flat layout directly." << endmsg;
        return StatusCode::FAILURE;
      }
      m_index = static_cast<int>(std::distance(names.begin(), it));
    }
    if (m_min.value() > m_max.value())
      return error() << "Min (" << m_min.value() << ") > Max (" << m_max.value()
                     << "): this cut can never pass" << endmsg,
             StatusCode::FAILURE;

    info() << "cut: " << m_variable.value() << " (shapeParameters[" << m_index << "]) in ["
           << m_min.value() << ", " << m_max.value() << "]" << endmsg;
    return StatusCode::SUCCESS;
  }

  bool operator()(const edm4hep::ClusterCollection& clusters) const override {
    ++m_seen;
    // One Cluster per event by construction (EcalPidTransformer is 1->1). No
    // cluster means no variables to cut on, so nothing to keep.
    if (clusters.size() != 1) return false;

    const auto params = clusters[0].getShapeParameters();
    if (m_index >= static_cast<int>(params.size())) {
      // Not fatal per event, but it means the cut is not being applied at all.
      error() << "ParamIndex " << m_index << " out of range (" << params.size()
              << " shapeParameters): rejecting the event" << endmsg;
      return false;
    }

    const float value = params[m_index];
    const bool pass = value >= m_min.value() && value <= m_max.value();  // NaN -> false
    if (pass) ++m_passed;
    return pass;
  }

  StatusCode finalize() override {
    const auto seen = m_seen.load();
    const auto passed = m_passed.load();
    // Every cut in a chain sees EVERY event: k4FWCore wraps the algorithms in a
    // sequencer with ShortCircuit=False, so an earlier rejection does not stop
    // the later ones from running. This count is therefore THIS cut alone, not
    // the running total -- the events actually written are the ones that passed
    // all the cuts, which is <= the smallest of these numbers.
    info() << m_variable.value() << " in [" << m_min.value() << ", " << m_max.value()
           << "]: this cut alone keeps " << passed << " / " << seen << " events"
           << (seen ? " (" + std::to_string(100.0 * passed / seen) + "%)" : "")
           << "; events written = those passing every cut in the chain" << endmsg;
    return FilterPredicate::finalize();
  }

  Gaudi::Property<std::string> m_variable{this, "Variable", "moliere",
                                          "Base scalar to cut on (EcalShowerVars scalarNames)"};
  // +-1e30 and not +-infinity: Gaudi serialises an infinite default into the
  // generated configurable as the literal "-inf.0", which is not valid Python,
  // and the confdb2 merge step dies parsing it. 1e30 is unreachable by any of
  // these variables, so it means "no bound" just as well.
  Gaudi::Property<double> m_min{this, "Min", -1e30, "Lower bound, inclusive (-1e30 = no bound)"};
  Gaudi::Property<double> m_max{this, "Max", 1e30, "Upper bound, inclusive (1e30 = no bound)"};
  Gaudi::Property<int> m_paramIndex{
      this, "ParamIndex", -1,
      "Address shapeParameters by index instead of by name (-1 = use Variable). For entries beyond "
      "the base scalars, whose position depends on NLayers and MipThresholds."};

private:
  int m_index = -1;
  mutable std::atomic<long> m_seen{0};
  mutable std::atomic<long> m_passed{0};
};

DECLARE_COMPONENT(EcalEventCut)
