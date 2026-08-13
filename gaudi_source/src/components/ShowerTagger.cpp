/**
 * ShowerTagger — identify electromagnetic showers and keep their hits out of
 * the ACTS measurement pool.
 *
 * Ported from siwecal_k4sim/gaudi_source/ShowerTagger.cpp (same class name,
 * same properties, same shower-onset logic and outputs), adapted to run on
 * this repo's real reconstructed hits instead of simulated ones:
 *
 *   - Input is edm4hep::CalorimeterHitCollection (this repo's `ECalHits`,
 *     written by EcalToEDM4hep), not edm4hep::SimCalorimeterHitCollection.
 *   - Default BitField decodes this repo's own cellID convention
 *     ("system:8,slab:8,chip:16,channel:8,sca:8", the same one
 *     EcalToEDM4hep.cpp writes) and the decoded field name is "slab"
 *     (siwecal_k4sim's DD4hep bitfield calls the equivalent field "layer").
 *
 * Why this exists
 * ---------------
 * A track through an electromagnetic cascade is not a physical object. Once a
 * particle starts showering, the pads it lights up are secondaries spraying
 * transversely, not samples of a trajectory; fitting a helix (or a line)
 * through them produces a well-formed track object with a small chi2 and no
 * physical meaning. The right answer for an EM event is "this is a shower",
 * with its energy and barycentre — not "these are N tracks".
 *
 * This is sharper in a SiW-ECAL than in a tracker-plus-calorimeter setup:
 * every layer here is 1.2-2.0 X0, so a high-energy electron is already
 * showering in slab 0 or 1 and there is no incoming segment left to fit.
 *
 * What it does
 * ------------
 * Per slab, count the hits. A MIP crossing lights one or two pads per slab;
 * a shower lights tens. `ShowerNHitsThreshold` consecutive-dense slabs
 * (`ShowerMinConsecutive` of them, so a single delta-ray fluctuation does not
 * count) mark the shower onset.
 *
 * Outputs, per event:
 *   OutputFlags   podio::UserDataCollection<int32_t>, one entry per input hit:
 *                 0 = track-like (may enter ACTS), 1 = shower hit (may not).
 *                 Consumed by SiPadMeasConverter's InputFlags property.
 *   OutputShowers edm4hep::ClusterCollection, one cluster per identified
 *                 shower: total energy, energy-weighted barycentre, and shape
 *                 parameters (see below). Empty for a non-showering event.
 *
 * Hits in slabs at or beyond the onset are flagged as shower hits. If fewer
 * than `MinTrackLayers` slabs precede the onset, ALL hits are flagged: a stub
 * of one to three points cannot define a trajectory, and letting it through is
 * how a shower event ends up with a "track" anyway. So a high-energy electron
 * yields one shower and zero tracks, while a muon yields no shower and its
 * fifteen-slab track is untouched. A pion that punches through and then
 * showers keeps its incoming track AND gets a shower — which is exactly the
 * information PID wants.
 *
 * shapeParameters layout (all floats):
 *   [0] shower start slab
 *   [1] slab of maximum energy
 *   [2] number of slabs spanned by the shower
 *   [3] transverse RMS about the barycentre [mm]
 *   [4] number of hits in the shower
 *
 * Energy units are whatever the input carries (MIPs, since ECalHits is
 * already calibrated by EcalEventBuilder/EcalToEDM4hep).
 *
 * Properties
 * ----------
 *   InputCollection        CalorimeterHit collection          (ECalHits)
 *   OutputFlags            per-hit shower flags                (SiPadShowerFlags)
 *   OutputShowers          reconstructed showers                (EMShowers)
 *   BitField               cellID decoder (needs a 'slab' field)
 *   ShowerNHitsThreshold   hits/slab for a slab to be dense (default 4)
 *   ShowerMinConsecutive   consecutive dense slabs for onset (default 2)
 *   MinTrackLayers         pre-shower slabs needed to still try tracking (4)
 *   Enabled                false disables tagging entirely (flags all zero)
 *   DebugFrequency         print every Nth event
 */

#include "k4FWCore/Transformer.h"

#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "podio/UserDataCollection.h"

#include "DDSegmentation/BitFieldCoder.h"

#include "Gaudi/Property.h"

#include <atomic>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

using OutputType = std::tuple<podio::UserDataCollection<std::int32_t>,
                              edm4hep::ClusterCollection>;

struct ShowerTagger final
    : k4FWCore::MultiTransformer<OutputType(const edm4hep::CalorimeterHitCollection&)> {

  ShowerTagger(const std::string& name, ISvcLocator* svc)
      : MultiTransformer(
            name, svc,
            {KeyValues("InputCollection", {"ECalHits"})},
            {KeyValues("OutputFlags", {"SiPadShowerFlags"}),
             KeyValues("OutputShowers", {"EMShowers"})}) {}

  StatusCode initialize() override {
    StatusCode sc = MultiTransformer::initialize();
    if (sc.isFailure()) return sc;

    if (m_bitField.value().empty()) {
      error() << "[ShowerTagger] BitField must be set." << endmsg;
      return StatusCode::FAILURE;
    }
    m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
        m_bitField.value());

    if (m_showerMinConsecutive.value() < 1) {
      error() << "[ShowerTagger] ShowerMinConsecutive must be >= 1." << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "[ShowerTagger] " << (m_enabled.value() ? "enabled" : "DISABLED")
           << ": a slab is dense at >= " << m_showerNHitsThreshold.value()
           << " hits, onset needs " << m_showerMinConsecutive.value()
           << " consecutive dense slabs, tracking needs "
           << m_minTrackLayers.value() << " slabs before the onset." << endmsg;
    return StatusCode::SUCCESS;
  }

  OutputType operator()(
      const edm4hep::CalorimeterHitCollection& input) const override {
    podio::UserDataCollection<std::int32_t> flags;
    edm4hep::ClusterCollection showers;

    const long long evtNum  = m_eventCount.fetch_add(1);
    const bool      doPrint = (m_debugFreq.value() > 0) &&
                              (evtNum % m_debugFreq.value() == 0);

    const std::size_t nHits = input.size();
    flags.vec().assign(nHits, 0);

    if (!m_enabled.value() || nHits == 0) {
      return std::make_tuple(std::move(flags), std::move(showers));
    }

    // ---- per-hit slab, and hits per slab --------------------------------
    std::vector<int> hitLayer(nHits, -1);
    std::map<int, int> nPerLayer;
    for (std::size_t i = 0; i < nHits; ++i) {
      const int slab = static_cast<int>(
          m_decoder->get(input[i].getCellID(), "slab"));
      hitLayer[i] = slab;
      ++nPerLayer[slab];
    }

    // ---- shower onset: first run of consecutive dense slabs -------------
    // "Consecutive" means consecutive slab NUMBERS, not consecutive entries
    // of nPerLayer: a slab with no hits at all breaks the run, which is what
    // we want (a shower does not skip slabs).
    const int  nThr  = m_showerNHitsThreshold.value();
    const int  runNeeded = m_showerMinConsecutive.value();
    const int  firstLayer = nPerLayer.begin()->first;
    const int  lastLayer  = nPerLayer.rbegin()->first;

    int startLayer = -1;
    int run = 0;
    for (int l = firstLayer; l <= lastLayer; ++l) {
      auto it = nPerLayer.find(l);
      const int n = (it == nPerLayer.end()) ? 0 : it->second;
      if (n >= nThr) {
        if (run == 0) run = 1; else ++run;
        if (run >= runNeeded) { startLayer = l - (run - 1); break; }
      } else {
        run = 0;
      }
    }

    if (startLayer < 0) {
      // No shower: every hit stays available to the tracker.
      if (doPrint) {
        debug() << "[ShowerTagger] evt=" << evtNum << " hits=" << nHits
                << " no shower found." << endmsg;
      }
      return std::make_tuple(std::move(flags), std::move(showers));
    }

    // ---- how many slabs precede the onset --------------------------------
    int nLayersBefore = 0;
    for (const auto& [l, n] : nPerLayer) {
      if (l < startLayer && n > 0) ++nLayersBefore;
    }
    // Too short a stub to be a trajectory -> nothing goes to ACTS at all.
    const bool vetoAll = (nLayersBefore < m_minTrackLayers.value());

    // ---- flag the shower hits, and accumulate the cluster ---------------
    double eTot = 0.0, sx = 0.0, sy = 0.0, sz = 0.0;
    std::map<int, double> ePerLayer;
    int nShowerHits = 0;

    for (std::size_t i = 0; i < nHits; ++i) {
      const bool inShower = (hitLayer[i] >= startLayer);
      if (inShower || vetoAll) flags.vec()[i] = 1;
      if (!inShower) continue;

      const auto&  hit = input[i];
      const double e   = hit.getEnergy();
      const auto&  p   = hit.getPosition();
      eTot += e;
      sx   += e * p.x;
      sy   += e * p.y;
      sz   += e * p.z;
      ePerLayer[hitLayer[i]] += e;
      ++nShowerHits;
    }

    if (eTot <= 0.0) {
      // Dense but with no energy: nothing sane to report as a shower, and the
      // flags already keep those hits out of the fit.
      return std::make_tuple(std::move(flags), std::move(showers));
    }

    const double bx = sx / eTot, by = sy / eTot, bz = sz / eTot;

    // Transverse RMS about the barycentre, energy weighted.
    double s2 = 0.0;
    for (std::size_t i = 0; i < nHits; ++i) {
      if (hitLayer[i] < startLayer) continue;
      const auto&  hit = input[i];
      const double e   = hit.getEnergy();
      const auto&  p   = hit.getPosition();
      s2 += e * ((p.x - bx) * (p.x - bx) + (p.y - by) * (p.y - by));
    }
    const double rmsT = std::sqrt(s2 / eTot);

    int layerMax = startLayer;
    double eMax  = -1.0;
    for (const auto& [l, e] : ePerLayer) {
      if (e > eMax) { eMax = e; layerMax = l; }
    }

    auto cluster = showers.create();
    cluster.setType(1);  // 1 = electromagnetic shower
    cluster.setEnergy(static_cast<float>(eTot));
    cluster.setPosition({static_cast<float>(bx),
                         static_cast<float>(by),
                         static_cast<float>(bz)});
    cluster.addToShapeParameters(static_cast<float>(startLayer));
    cluster.addToShapeParameters(static_cast<float>(layerMax));
    cluster.addToShapeParameters(static_cast<float>(ePerLayer.size()));
    cluster.addToShapeParameters(static_cast<float>(rmsT));
    cluster.addToShapeParameters(static_cast<float>(nShowerHits));

    if (doPrint) {
      debug() << "[ShowerTagger] evt=" << evtNum << " hits=" << nHits
              << " shower starts at slab " << startLayer
              << " (" << nLayersBefore << " slabs before it)"
              << " nShowerHits=" << nShowerHits
              << " E=" << eTot
              << " layerMax=" << layerMax
              << " rmsT=" << rmsT << " mm"
              << (vetoAll ? " -- pre-shower stub too short, no tracking"
                          : " -- keeping the pre-shower segment for tracking")
              << endmsg;
    }

    return std::make_tuple(std::move(flags), std::move(showers));
  }

private:
  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,slab:8,chip:16,channel:8,sca:8",
      "DD4hep BitField string for the input cellIDs; needs a 'slab' field. "
      "Defaults to the same convention EcalToEDM4hep writes (see its "
      "CellIDEncoding property)."};
  Gaudi::Property<int> m_showerNHitsThreshold{
      this, "ShowerNHitsThreshold", 4,
      "Hits in one slab for that slab to count as dense. A MIP gives 1-2."};
  Gaudi::Property<int> m_showerMinConsecutive{
      this, "ShowerMinConsecutive", 2,
      "Consecutive dense slabs required to declare a shower onset. 2 rejects "
      "a single-slab delta-ray fluctuation."};
  Gaudi::Property<int> m_minTrackLayers{
      this, "MinTrackLayers", 4,
      "Slabs that must precede the shower onset for the pre-shower segment to "
      "be offered to the tracker. Below this, every hit is flagged, so the "
      "event yields a shower and no track."};
  Gaudi::Property<bool> m_enabled{
      this, "Enabled", true,
      "If false, no tagging: all flags 0 and no showers written."};
  Gaudi::Property<int> m_debugFreq{
      this, "DebugFrequency", 500, "Print every Nth event (0 = never)."};

  mutable std::atomic<long long> m_eventCount{0};
  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
};

DECLARE_COMPONENT(ShowerTagger)
