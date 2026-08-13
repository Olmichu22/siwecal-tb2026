#pragma once

#include "k4SiWEcalReco/ISNDGeoSvc.h"
#include "GaudiKernel/Service.h"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

// ACTSGeoSvc — same public interface, same address convention (detID=1
// "SiPad", station=-1, layer=hit_slab 0..14, plane=-1) and same rot90Y
// beam-axis convention as siwecal_k4sim's gaudi_source/ACTSGeoSvc, so
// ACTSProtoTracker is a verbatim port between the two repos.
//
// The one real difference: siwecal_k4sim builds its 15 plane surfaces by
// walking a DD4hep compact-XML detector description. This repo has no such
// description -- real geometry here is mappings/slab_z_positions.yml (per
// hit_slab z and W thickness), already ported into C++ as
// k4siwecal::SlabGeometry (include/k4SiWEcalReco/PadMapGeometry.h). So this
// service builds the same 15 PlaneSurface/PlaneLayer scaffolding directly
// from SlabGeometry instead of a DetElement tree.
class ACTSGeoSvc : public extends<Service, ISNDGeoSvc> {
public:
  ACTSGeoSvc(const std::string& name, ISvcLocator* svcLoc);

  StatusCode initialize() override;
  StatusCode finalize()   override;

  // IActsGeoSvc interface
  const Acts::TrackingGeometry& trackingGeometry() const override;

  // ISNDGeoSvc extensions
  const std::vector<const Acts::Surface*>& allSurfaces() const override;
  const Acts::GeometryContext& geometryContext() const override;
  const Acts::Surface* surfaceByAddress(
      int detID, int station, int layer, int plane) const override;

private:
  Gaudi::Property<std::string> m_slabZFile{
      this, "SlabZFile", "",
      "Path to mappings/slab_z_positions.yml (absolute -- the job may run "
      "from any directory). Replaces siwecal_k4sim's CompactFile: this repo "
      "has no DD4hep detector description, so the 15 surfaces are built "
      "directly from this file instead of a compact XML."};

  Gaudi::Property<double> m_padPitchMm{
      this, "PadPitchMm", 5.53,
      "Pad pitch [mm] -- same physical hardware as siwecal_k4sim "
      "(event_viewer/model/detector.py:PAD_PITCH_MM there)."};

  Gaudi::Property<double> m_halfSizeX{
      this, "HalfSizeX", 90.0,
      "Transverse half-size in X [mm] of each plane surface (detector envelope)."};

  Gaudi::Property<double> m_halfSizeY{
      this, "HalfSizeY", 90.0,
      "Transverse half-size in Y [mm] of each plane surface (detector envelope)."};

  std::shared_ptr<const Acts::TrackingGeometry> m_trackingGeometry;
  std::vector<const Acts::Surface*>             m_allSurfaces;
  Acts::GeometryContext                         m_gctx;
  // Detector elements must outlive all surfaces (surfaces hold raw pointers back).
  std::vector<std::shared_ptr<Acts::DetectorElementBase>> m_detectorElements;

  // key = (detID, station, layerInDet, plane)
  std::map<std::tuple<int,int,int,int>, const Acts::Surface*> m_surfaceByAddressMap;
};
