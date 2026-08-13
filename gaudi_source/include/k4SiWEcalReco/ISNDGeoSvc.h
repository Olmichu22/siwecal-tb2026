#pragma once
#include "k4ActsTracking/IActsGeoSvc.h"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include <vector>

class ISNDGeoSvc : virtual public IActsGeoSvc {
public:
  DeclareInterfaceID(ISNDGeoSvc, 1, 0);

  // Returns all rectangle surfaces sorted by X, for all detectors combined.
  virtual const std::vector<const Acts::Surface*>& allSurfaces() const = 0;

  // Returns the ACTS geometry context (needed by algorithms in execute()).
  virtual const Acts::GeometryContext& geometryContext() const = 0;

  // Returns the ACTS surface for the given detector address, or nullptr.
  // detID: 1=SiPad
  // station: -1 (SiPad has no stations)
  // layer:   layer index within detector (= hit_slab, 0..14)
  // plane:   -1 (SiPad pixels)
  virtual const Acts::Surface* surfaceByAddress(
      int detID, int station, int layer, int plane) const = 0;

  virtual ~ISNDGeoSvc() {}
};
