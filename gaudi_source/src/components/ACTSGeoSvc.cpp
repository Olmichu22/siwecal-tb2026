#include "k4SiWEcalReco/ACTSGeoSvc.h"
#include "k4SiWEcalReco/PadMapGeometry.h"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/PlaneLayer.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/CuboidVolumeBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include <algorithm>
#include <cmath>
#include <string>

// ---------------------------------------------------------------------------
// SNDDetectorElement
// A minimal DetectorElementBase so the CKF actor classifies each plane as
// sensitive (CKF checks associatedDetectorElement() != nullptr at line 442
// of CombinatorialKalmanFilter.hpp; without this, every surface is "passive").
// Identical to siwecal_k4sim/gaudi_source/ACTSGeoSvc.cpp.
// ---------------------------------------------------------------------------

class SNDDetectorElement : public Acts::DetectorElementBase {
public:
  SNDDetectorElement(std::shared_ptr<const Acts::PlanarBounds> bounds,
                     Acts::Transform3 transform, double thickness)
      : m_transform(std::move(transform)), m_thickness(thickness) {
    // PlaneSurface(bounds, detElement) sets m_associatedDetectorElement = this
    m_surface = Acts::Surface::makeShared<Acts::PlaneSurface>(bounds, *this);
  }

  const Acts::Transform3& transform(const Acts::GeometryContext&) const override {
    return m_transform;
  }
  const Acts::Surface& surface() const override { return *m_surface; }
  Acts::Surface&       surface() override       { return *m_surface; }
  double thickness() const override { return m_thickness; }

private:
  Acts::Transform3 m_transform;
  double m_thickness;
  std::shared_ptr<Acts::PlaneSurface> m_surface;
};

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------

ACTSGeoSvc::ACTSGeoSvc(const std::string& name, ISvcLocator* svcLoc)
    : extends(name, svcLoc) {}

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

namespace {
// Tungsten material constants, consistent with
// k4SiWEcalReco/PadMapGeometry.h's kWX0Mm (3.5 mm radiation length) -- the
// same constant EcalEventBuilder uses for hit_X0. L0/A/Z/rho are the PDG
// values for W; L0 only affects hadronic-interaction bookkeeping, not the
// Highland multiple-scattering estimate the CKF actually uses (that one only
// reads X0), so approximate precision there is fine.
constexpr float kWL0Mm  = 99.4f;   // nuclear interaction length [mm]
constexpr float kWA     = 183.84f;
constexpr float kWZ     = 74.0f;
constexpr float kWRhoG  = 19.3f;   // g/cm^3

Acts::MaterialSlab makeTungstenSlab(double thicknessMm) {
  return Acts::MaterialSlab(
      Acts::Material::fromMassDensity(
          static_cast<float>(k4siwecal::kWX0Mm), kWL0Mm, kWA, kWZ, kWRhoG),
      static_cast<float>(thicknessMm));
}
} // namespace

StatusCode ACTSGeoSvc::initialize() {
  StatusCode sc = Service::initialize();
  if (sc.isFailure()) return sc;

  if (m_slabZFile.value().empty()) {
    error() << "[ACTSGeoSvc] SlabZFile property must be set." << endmsg;
    return StatusCode::FAILURE;
  }

  k4siwecal::SlabGeometry slabGeo = k4siwecal::SlabGeometry::fromYamlFile(m_slabZFile.value());
  info() << "[ACTSGeoSvc] Slab geometry loaded from: " << m_slabZFile.value() << endmsg;

  // Struct holding extracted plane info -- same fields siwecal_k4sim's
  // ACTSGeoSvc fills from the DD4hep DetElement tree; here they come from
  // SlabGeometry (mappings/slab_z_positions.yml) plus the configured pad
  // pitch / envelope instead.
  struct PlaneInfo {
    double z;           // global Z center [mm]  (beam axis, before rot90Y)
    double halfX;        // half-size in X [mm]   (transverse)
    double halfY;        // half-size in Y [mm]   (transverse)
    double thickness;    // half-thickness in Z [mm]
    int    plane;        // -1 = pixel (SiPad)
    int    detID;        // 1 = SiPad
    int    station;      // -1 (SiPad has no stations)
    int    layerInDet;   // hit_slab index
    Acts::MaterialSlab material;  // W absorber only (dominates the budget)
  };

  const double halfX = m_halfSizeX.value();
  const double halfY = m_halfSizeY.value();
  // Half-thickness placeholder [mm]: this repo has no per-slab Si/PCB/Cu
  // slice thicknesses to sum (unlike siwecal_k4sim's DD4hep slice stack), so
  // the plane is given a thin nominal thickness -- what matters for the CKF
  // is the *surface material* (set below from the W budget), not this value.
  constexpr double kNominalHalfThicknessMm = 0.5;

  std::vector<PlaneInfo> allPlanes;
  constexpr int kNumSlabs = 15;
  for (int slab = 0; slab < kNumSlabs; ++slab) {
    const double z = slabGeo.slabZ(slab);
    if (!std::isfinite(z)) {
      error() << "[ACTSGeoSvc] hit_slab " << slab
              << " has no z position in " << m_slabZFile.value() << endmsg;
      return StatusCode::FAILURE;
    }
    PlaneInfo pi{};
    pi.z          = z;
    pi.halfX      = halfX;
    pi.halfY      = halfY;
    pi.thickness  = kNominalHalfThicknessMm;
    pi.plane      = -1;
    pi.detID      = 1;   // SiPad -- same address convention as siwecal_k4sim
    pi.station    = -1;
    pi.layerInDet = slab;
    // slabWOverX0() = w_thickness_mm[slab] / kWX0Mm, so the thickness this
    // slab's W absorber has is slabWOverX0() * kWX0Mm -- recovers the raw
    // mm value from the file without adding a new getter to PadMapGeometry.h.
    const double wThicknessMm = slabGeo.slabWOverX0(slab) * k4siwecal::kWX0Mm;
    pi.material = makeTungstenSlab(wThicknessMm);
    allPlanes.push_back(pi);
  }

  // Already built in hit_slab order == beam-z order for this detector
  // (monotonic downstream), but sort defensively as siwecal_k4sim does.
  std::sort(allPlanes.begin(), allPlanes.end(),
            [](const PlaneInfo& a, const PlaneInfo& b) {
              return a.z < b.z;
            });

  info() << "[ACTSGeoSvc] Found " << allPlanes.size() << " SiPad planes." << endmsg;
  for (const auto& pi : allPlanes) {
    info() << "[ACTSGeoSvc]   slab=" << pi.layerInDet
           << " z=" << pi.z
           << " halfX=" << pi.halfX << " halfY=" << pi.halfY
           << " | material t=" << pi.material.thickness()
           << " mm X0=" << pi.material.material().X0()
           << " mm t/X0=" << pi.material.thicknessInX0() << endmsg;
  }

  // =========================================================================
  // Build TrackingGeometry using PlaneLayer directly (same scaffolding as
  // siwecal_k4sim's ACTSGeoSvc: avoids CuboidVolumeBuilder's position offset
  // bug). Each plane becomes a PlaneSurface + PlaneLayer with absolute
  // position; all layers sit in one CuboidVolumeBounds TrackingVolume.
  // =========================================================================

  std::vector<Acts::LayerPtr> allLayers;
  allLayers.reserve(allPlanes.size());

  // Rotation 90 deg around Y: maps local Z (surface normal) to global X, so
  // tracks run along ACTS X and avoid the theta=0 bound-coordinate
  // singularity. Same convention as siwecal_k4sim -- local x maps to minus
  // global Z, see docs/acts_integration.md's coordinate section.
  Acts::RotationMatrix3 rot90Y =
      Eigen::AngleAxisd(M_PI / 2.0,
                        Acts::Vector3::UnitY()).toRotationMatrix();

  for (std::size_t i = 0; i < allPlanes.size(); ++i) {
    const auto& pi = allPlanes[i];

    Acts::Transform3 transform = Acts::Transform3::Identity();
    transform.rotate(rot90Y);
    transform.pretranslate(Acts::Vector3(pi.z, 0.0, 0.0));

    auto bounds = std::make_shared<Acts::RectangleBounds>(pi.halfX, pi.halfY);
    const double layerThickness = 2.0 * pi.thickness + 0.01;

    auto detElem = std::make_shared<SNDDetectorElement>(bounds, transform, pi.thickness);
    m_detectorElements.push_back(detElem);

    auto surfaceArray = std::make_unique<Acts::SurfaceArray>(
        detElem->surface().getSharedPtr());

    auto layer = Acts::PlaneLayer::create(
        transform, bounds, std::move(surfaceArray), layerThickness,
        nullptr, Acts::active);

    if (pi.material.thickness() > 0.0) {
      auto surfMat = std::make_shared<Acts::HomogeneousSurfaceMaterial>(pi.material);
      detElem->surface().assignSurfaceMaterial(surfMat);
    } else {
      warning() << "[ACTSGeoSvc] hit_slab " << pi.layerInDet
                << " has no material budget -- surface left material-free."
                << endmsg;
    }

    allLayers.push_back(layer);
  }

  // ---- Navigation layers at the beam-Z extremes -----------------------
  const double zMin = allPlanes.front().z - allPlanes.front().thickness - 5.0;
  const double zMax = allPlanes.back().z  + allPlanes.back().thickness  + 5.0;

  {
    Acts::Transform3 t = Acts::Transform3::Identity();
    t.rotate(rot90Y);
    t.pretranslate(Acts::Vector3(zMin, 0.0, 0.0));
    auto navBounds = std::make_shared<Acts::RectangleBounds>(halfX + 5.0, halfY + 5.0);
    auto navLayer = Acts::PlaneLayer::create(t, navBounds, nullptr, 0.0, nullptr, Acts::navigation);
    allLayers.insert(allLayers.begin(), navLayer);
  }
  {
    Acts::Transform3 t = Acts::Transform3::Identity();
    t.rotate(rot90Y);
    t.pretranslate(Acts::Vector3(zMax, 0.0, 0.0));
    auto navBounds = std::make_shared<Acts::RectangleBounds>(halfX + 5.0, halfY + 5.0);
    auto navLayer = Acts::PlaneLayer::create(t, navBounds, nullptr, 0.0, nullptr, Acts::navigation);
    allLayers.push_back(navLayer);
  }

  // ---- LayerArray, binned along X (beam axis) --------------------------
  Acts::LayerArrayCreator::Config lacCfg{};
  Acts::LayerArrayCreator lac(
      lacCfg, Acts::getDefaultLogger("LayerArrayCreator", Acts::Logging::WARNING));

  auto layerArray = lac.layerArray(
      m_gctx, allLayers, zMin - 1.0, zMax + 1.0,
      Acts::arbitrary, Acts::AxisDirection::AxisX);

  // ---- TrackingVolume ----------------------------------------------------
  const double volumeHalfZ = halfX + 15.0;
  const double volumeHalfY = halfY + 15.0;
  const double volumeHalfX = std::max(std::abs(zMin), std::abs(zMax)) + 2.0;

  auto volumeBounds = std::make_shared<Acts::CuboidVolumeBounds>(
      volumeHalfX, volumeHalfY, volumeHalfZ);

  Acts::Transform3 volumeTransform = Acts::Transform3::Identity();

  auto trackingVolume = std::make_shared<Acts::TrackingVolume>(
    volumeTransform, volumeBounds, nullptr, std::move(layerArray),
    nullptr, Acts::MutableTrackingVolumeVector{}, "SNDVolume");

  // ---- TrackingGeometry ----------------------------------------------------
  Acts::GeometryIdentifierHook hook{};
  m_trackingGeometry = std::make_shared<Acts::TrackingGeometry>(
    trackingVolume, nullptr, hook, Acts::getDummyLogger());

  if (!m_trackingGeometry) {
    error() << "[ACTSGeoSvc] Failed to build ACTS TrackingGeometry." << endmsg;
    return StatusCode::FAILURE;
  }
  info() << "[ACTSGeoSvc] TrackingGeometry built with PlaneLayer." << endmsg;

  // Collect sensitive surfaces from confined layers -- same idiom as
  // siwecal_k4sim: use the module surfaces stored in each layer's
  // SurfaceArray (the ones the Navigator actually visits), not
  // surfaceRepresentation() (that would give the layer envelope's geoID).
  const Acts::TrackingVolume* world = m_trackingGeometry->highestTrackingVolume();
  if (world && world->confinedLayers()) {
    for (const auto& layer : world->confinedLayers()->arrayObjects()) {
      if (layer && layer->layerType() == Acts::active) {
        const Acts::SurfaceArray* sa = layer->surfaceArray();
        if (sa) {
          for (const Acts::Surface* sf : sa->surfaces()) {
            if (sf) m_allSurfaces.push_back(sf);
          }
        } else {
          m_allSurfaces.push_back(&layer->surfaceRepresentation());
        }
      }
    }
  }
  std::sort(m_allSurfaces.begin(), m_allSurfaces.end(),
            [&](const Acts::Surface* a, const Acts::Surface* b) {
              return a->center(m_gctx).x() < b->center(m_gctx).x();
            });

  if (allPlanes.size() == m_allSurfaces.size()) {
    for (std::size_t i = 0; i < allPlanes.size(); ++i) {
      const auto& pi = allPlanes[i];
      m_surfaceByAddressMap[{pi.detID, pi.station, pi.layerInDet, pi.plane}]
          = m_allSurfaces[i];
    }
  } else {
    warning() << "[ACTSGeoSvc] allPlanes.size()=" << allPlanes.size()
              << " != m_allSurfaces.size()=" << m_allSurfaces.size()
              << " -- surfaceByAddress map not populated." << endmsg;
  }

  info() << "[ACTSGeoSvc] TrackingGeometry built. Total surfaces: "
         << m_allSurfaces.size() << endmsg;

  return StatusCode::SUCCESS;
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode ACTSGeoSvc::finalize() {
  return Service::finalize();
}

// ---------------------------------------------------------------------------
// Interface method implementations
// ---------------------------------------------------------------------------

const Acts::TrackingGeometry& ACTSGeoSvc::trackingGeometry() const {
  return *m_trackingGeometry;
}

const std::vector<const Acts::Surface*>& ACTSGeoSvc::allSurfaces() const {
  return m_allSurfaces;
}

const Acts::Surface* ACTSGeoSvc::surfaceByAddress(
    int detID, int station, int layer, int plane) const {
  auto it = m_surfaceByAddressMap.find({detID, station, layer, plane});
  return (it != m_surfaceByAddressMap.end()) ? it->second : nullptr;
}

const Acts::GeometryContext& ACTSGeoSvc::geometryContext() const {
  return m_gctx;
}

DECLARE_COMPONENT(ACTSGeoSvc)
