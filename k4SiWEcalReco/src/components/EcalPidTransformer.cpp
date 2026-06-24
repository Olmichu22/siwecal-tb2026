/*
 * EcalPidTransformer: per-event particle-discrimination variables for the
 * SiW-ECAL. Strict 1->1 functional Transformer: one input CalorimeterHit
 * collection (one event) -> one Cluster carrying the shower variables.
 *
 * The physics lives in EcalShowerVars.h (a C++ port of metrics.py). This
 * algorithm only marshals the EDM4hep hits into plain vectors, calls
 * computeEventVars, and flattens the result into Cluster::shapeParameters.
 *
 * shapeParameters layout (canonical, see EcalShowerVars.h):
 *   [ scalarNames()                       ]   base variables
 *   [ hits_per_layer[n] | energy_per_layer[n] | weighte_per_layer[n] ]
 *   [ <mipPrefix>_<scalarNames()> ... ]       one block per MIP-cut threshold
 * The matching shapeParameterNames vector is attached as collection metadata.
 */
#include "k4FWCore/MetaDataHandle.h"
#include "k4FWCore/Transformer.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"

#include "DDSegmentation/BitFieldCoder.h"

#include "Gaudi/Property.h"

#include "k4SiWEcalReco/EcalShowerVars.h"

#include <memory>
#include <string>
#include <vector>

using namespace k4siwecal;

namespace {
/// Branch-name-safe prefix for a MIP threshold: 0.5 -> "mip05", 1.0 -> "mip1".
std::string mipPrefix(float thr) {
  std::string s = std::to_string(thr);          // "0.500000"
  s = s.substr(0, s.find('.') + 2);             // "0.5"
  s.erase(std::remove(s.begin(), s.end(), '.'), s.end());  // "05"
  if (!s.empty() && s.back() == '0' && s.size() > 1) s.pop_back();  // "05"->"05","10"->"1"
  return "mip" + s;
}

void appendScalars(std::vector<float>& flat, const EventVars& v) {
  flat.insert(flat.end(),
              {v.nhit, v.zbary, v.energy, v.mip_likeness, v.weighte, v.bar_x, v.bar_y, v.bar_r,
               v.moliere, v.transverse_rms, v.is_shower, v.shower_start, v.shower_max, v.shower_end,
               v.shower_start_10, v.shower_end_10, v.shower_length, v.first_layer, v.last_layer,
               v.n_layers_hit, v.e_over_nhit});
}
}  // namespace

struct EcalPidTransformer final
    : k4FWCore::Transformer<edm4hep::ClusterCollection(const edm4hep::CalorimeterHitCollection&)> {
  EcalPidTransformer(const std::string& name, ISvcLocator* svcLoc)
      : Transformer(name, svcLoc, KeyValues("InputCaloHits", {"ECalHits"}),
                    KeyValues("OutputClusters", {"ECalPid"})) {}

  StatusCode initialize() override {
    m_coder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(m_cellIDEncoding.value());
    if (static_cast<int>(m_wThicknesses.value().size()) != m_nLayers.value())
      return error() << "WThicknesses length != NLayers" << endmsg, StatusCode::FAILURE;
    m_wOverX0.clear();
    for (double w : m_wThicknesses.value()) m_wOverX0.push_back(w / m_wX0.value());

    // Build and publish the shapeParameterNames metadata.
    m_paramNames = buildParamNames();
    m_shapeParamNames.put(m_paramNames);
    return StatusCode::SUCCESS;
  }

  edm4hep::ClusterCollection operator()(
      const edm4hep::CalorimeterHitCollection& hitColl) const override {
    edm4hep::ClusterCollection clusters;
    auto cluster = clusters.create();

    std::vector<int> slab;
    std::vector<float> energy, x, y;
    slab.reserve(hitColl.size());
    energy.reserve(hitColl.size());
    x.reserve(hitColl.size());
    y.reserve(hitColl.size());
    for (const auto& hit : hitColl) {
      slab.push_back(static_cast<int>(m_coder->get(hit.getCellID(), "slab")));
      energy.push_back(hit.getEnergy());
      x.push_back(hit.getPosition().x);
      y.push_back(hit.getPosition().y);
      cluster.addToHits(hit);
    }

    const ShowerThresholds thr{static_cast<float>(m_showerEThreshold.value()),
                               static_cast<float>(m_showerMaxMin.value()),
                               static_cast<float>(m_showerStartFrac.value())};
    const EventVars base = computeEventVars(slab, energy, x, y, m_wOverX0, m_nLayers.value(),
                                            m_showerProfile.value(), thr,
                                            m_moliereContainment.value());

    cluster.setEnergy(base.energy);
    cluster.setPosition({base.bar_x, base.bar_y, 0.f});

    std::vector<float> flat;
    appendScalars(flat, base);
    flat.insert(flat.end(), base.hits_per_layer.begin(), base.hits_per_layer.end());
    flat.insert(flat.end(), base.energy_per_layer.begin(), base.energy_per_layer.end());
    flat.insert(flat.end(), base.weighte_per_layer.begin(), base.weighte_per_layer.end());

    // MIP-cut variants: re-filter hits by energy >= threshold, scalars only.
    for (float mip : m_mipThresholds.value()) {
      std::vector<int> fs;
      std::vector<float> fe, fx, fy;
      for (std::size_t i = 0; i < slab.size(); ++i)
        if (energy[i] >= mip) {
          fs.push_back(slab[i]);
          fe.push_back(energy[i]);
          fx.push_back(x[i]);
          fy.push_back(y[i]);
        }
      const EventVars vv = computeEventVars(fs, fe, fx, fy, m_wOverX0, m_nLayers.value(),
                                            m_showerProfile.value(), thr,
                                            m_moliereContainment.value());
      appendScalars(flat, vv);
    }

    for (float value : flat) cluster.addToShapeParameters(value);
    return clusters;
  }

  /// Canonical shapeParameterNames matching the flat layout above.
  std::vector<std::string> buildParamNames() const {
    std::vector<std::string> names = scalarNames();
    for (const auto& block : perLayerNames())
      for (int i = 0; i < m_nLayers.value(); ++i) names.push_back(block + "_" + std::to_string(i));
    for (float mip : m_mipThresholds.value()) {
      const std::string pre = mipPrefix(mip);
      for (const auto& s : scalarNames()) names.push_back(pre + "_" + s);
    }
    return names;
  }

  Gaudi::Property<int> m_nLayers{this, "NLayers", 15, "Number of ECAL layers"};
  // double (not float) to bit-match metrics.py's float64 W/X0 weights.
  Gaudi::Property<std::vector<double>> m_wThicknesses{
      this, "WThicknesses",
      {2.8, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 4.2, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.6},
      "Per-slab tungsten thickness [mm]"};
  Gaudi::Property<double> m_wX0{this, "WX0", 3.5, "Tungsten radiation length [mm]"};
  Gaudi::Property<std::string> m_showerProfile{this, "ShowerProfile", "nhit",
                                               "Profile for shower flag: nhit|sume|weighte"};
  Gaudi::Property<float> m_showerEThreshold{this, "ShowerEThreshold", 5.0f, ""};
  Gaudi::Property<float> m_showerMaxMin{this, "ShowerMaxMin", 10.0f, ""};
  Gaudi::Property<float> m_showerStartFrac{this, "ShowerStartFrac", 0.1f, ""};
  Gaudi::Property<double> m_moliereContainment{this, "MoliereContainment", 0.90, ""};
  Gaudi::Property<std::vector<float>> m_mipThresholds{this, "MipThresholds", {0.5f, 1.0f}, ""};
  Gaudi::Property<std::string> m_cellIDEncoding{
      this, "CellIDEncoding", "system:8,slab:8,chip:16,channel:8,sca:8", ""};

private:
  std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_coder;
  std::vector<double> m_wOverX0;
  std::vector<std::string> m_paramNames;
  // Frame-level metadata: the names matching the flat shapeParameters layout.
  // Read downstream to interpret Cluster::shapeParameters by name.
  k4FWCore::MetaDataHandle<std::vector<std::string>> m_shapeParamNames{
      "ECalPid_shapeParameterNames", Gaudi::DataHandle::Writer};
};

DECLARE_COMPONENT(EcalPidTransformer)
