// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/DataHandle.h"

// edm4hep input — this repo's real reconstructed hits (EcalToEDM4hep's
// ECalHits), not simulated ones. siwecal_k4sim's SiPadMeasConverter reads
// edm4hep::SimCalorimeterHitCollection instead.
#include "edm4hep/CalorimeterHitCollection.h"
#include "podio/UserDataCollection.h"

// edm4hep output
#include "edm4hep/TrackerHit3DCollection.h"
#include "edm4hep/MutableTrackerHit3D.h"
#include "edm4hep/CovMatrix3f.h"
#include "edm4hep/Vector3d.h"

// DD4hep segmentation
#include "DDSegmentation/BitFieldCoder.h"

// ACTS
#include "Acts/Definitions/Units.hpp"

// SND geometry service
#include "k4SiWEcalReco/ISNDGeoSvc.h"

// Standard
#include <atomic>
#include <cmath>
#include <memory>
#include <string>

// ---------------------------------------------------------------------------
// SiPadMeasConverter
//
// Ported from siwecal_k4sim/gaudi_source/SiPadMeasConverter.cpp -- same class
// name, same properties, same output schema (TrackerHit3DCollection,
// quality=slab, covariance=pitch^2/12). Reads edm4hep::CalorimeterHitCollection
// (this repo's ECalHits) instead of edm4hep::SimCalorimeterHitCollection, and
// takes the hit time directly from hit.getTime() -- a real field on
// CalorimeterHit -- instead of the MC-truth-only
// hit.getContributions(0).getTime() the simulation uses.
// ---------------------------------------------------------------------------

class SiPadMeasConverter : public Gaudi::Algorithm {
public:
  SiPadMeasConverter(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  // ---- Properties ----
  Gaudi::Property<std::string> m_inputCollection{
      this, "InputCollection", "ECalHits",
      "Input CalorimeterHitCollection name"};

  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "SiPadMeasurements",
      "Output TrackerHit3DCollection name"};

  Gaudi::Property<std::string> m_bitField{
      this, "BitField", "system:8,slab:8,chip:16,channel:8,sca:8",
      "DD4hep BitField string for the input cellIDs; needs a 'slab' field. "
      "Defaults to the same convention EcalToEDM4hep writes."};

  Gaudi::Property<std::string> m_inputFlags{
      this, "InputFlags", "",
      "Optional per-hit veto flags (podio::UserDataCollection<int32_t>, one "
      "entry per input hit, non-zero = drop). Set to ShowerTagger's "
      "OutputFlags so that hits belonging to an electromagnetic cascade never "
      "become ACTS measurements: a track through a shower is not a physical "
      "object. Empty disables the veto."};

  Gaudi::Property<double> m_pixelSizeX{
      this, "PixelSizeX", 5.53,
      "Pixel pitch in X [mm] — sets the measurement variance"};
  Gaudi::Property<double> m_pixelSizeY{
      this, "PixelSizeY", 5.53,
      "Pixel pitch in Y [mm] — sets the measurement variance"};

  // ---- DataHandles ----
  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection>>
      m_inputHandle;

  mutable std::unique_ptr<
      k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_outputHandle;

  mutable std::unique_ptr<
      k4FWCore::DataHandle<podio::UserDataCollection<std::int32_t>>>
      m_flagsHandle;

  // ---- Services ----
  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  // ---- Internal state ----
  mutable std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  mutable std::atomic<long long> m_eventCount{0};
};

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode SiPadMeasConverter::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    // Build DataHandles
    m_inputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection>>(
        m_inputCollection.value(), Gaudi::DataHandle::Reader, this);

    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    if (!m_inputFlags.value().empty()) {
      m_flagsHandle = std::make_unique<
          k4FWCore::DataHandle<podio::UserDataCollection<std::int32_t>>>(
          m_inputFlags.value(), Gaudi::DataHandle::Reader, this);
      info() << "[SiPadMeasConverter] Vetoing hits flagged in '"
             << m_inputFlags.value() << "'." << endmsg;
    }

    // Build BitField decoder
    m_decoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(
        m_bitField.value());

    // Retrieve geometry service
    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[SiPadMeasConverter] Failed to retrieve ACTSGeoSvc."
              << endmsg;
      return StatusCode::FAILURE;
    }

    const std::size_t nSurfaces = m_geoSvc->allSurfaces().size();
    info() << "[SiPadMeasConverter] Initialized. GeoSvc has "
           << nSurfaces << " surfaces." << endmsg;

    return sc;
  } catch (const std::exception& e) {
    error() << "[SiPadMeasConverter] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[SiPadMeasConverter] Unknown exception in initialize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode SiPadMeasConverter::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    auto* output = m_outputHandle->createAndPut();

    const auto* hits = m_inputHandle->get();
    if (!hits || hits->size() == 0) {
      debug() << "[SiPadMeasConverter] evt=" << evtNum
              << " no SiPad hits." << endmsg;
      return StatusCode::SUCCESS;
    }

    // SiPad is the only detector; surfaces are indexed directly by slab.
    const auto& allSurfaces       = m_geoSvc->allSurfaces();
    const std::size_t sipadOffset = 0;

    // Precompute pixel variances
    // Z=beam: SiPad measures X and Y (transverse to beam)
    const double varX = (m_pixelSizeX.value() * m_pixelSizeX.value()) / 12.0;
    const double varY = (m_pixelSizeY.value() * m_pixelSizeY.value()) / 12.0;

    // Build CovMatrix3f: xx and yy diagonal elements
    // Upper triangle storage: [xx=0, xy=1, xz=2, yy=3, yz=4, zz=5]
    edm4hep::CovMatrix3f cov{};
    cov[0] = static_cast<float>(varX);  // xx
    cov[3] = static_cast<float>(varY);  // yy
    // Same covariance for all SiPad hits — computed once outside the loop

    // Optional shower veto. The flags are one-per-input-hit and positional, so
    // a size mismatch means the flags were produced from a different
    // collection; refusing to guess is safer than silently vetoing the wrong
    // hits (or none).
    const podio::UserDataCollection<std::int32_t>* vetoFlags = nullptr;
    if (m_flagsHandle) {
      vetoFlags = m_flagsHandle->get();
      if (vetoFlags && vetoFlags->size() != hits->size()) {
        error() << "[SiPadMeasConverter] evt=" << evtNum << " flag collection '"
                << m_inputFlags.value() << "' has " << vetoFlags->size()
                << " entries but the input has " << hits->size()
                << " hits — they must be parallel." << endmsg;
        return StatusCode::FAILURE;
      }
    }

    std::size_t nConverted = 0;
    std::size_t nSkipped   = 0;
    std::size_t nVetoed    = 0;

    for (std::size_t i = 0; i < hits->size(); ++i) {
      if (vetoFlags && (*vetoFlags)[i] != 0) { ++nVetoed; continue; }
      const auto&    hit = (*hits)[i];
      const uint64_t cid = hit.getCellID();

      // Decode slab from cellID (no plane field for SiPad)
      int slab = static_cast<int>(m_decoder->get(cid, "slab"));

      const std::size_t surfIdx = sipadOffset +
          static_cast<std::size_t>(slab);

      if (surfIdx >= allSurfaces.size()) {
        warning() << "[SiPadMeasConverter] evt=" << evtNum
                  << " hit=" << i << " slab=" << slab
                  << " surfIdx=" << surfIdx << " out of range, skipping."
                  << endmsg;
        ++nSkipped;
        continue;
      }

      const auto& pos = hit.getPosition();

      // Create TrackerHit3D
      auto mhit = output->create();
      mhit.setCellID(cid);
      mhit.setType(1);      // detector ID: 1 = SiPad
      mhit.setQuality(slab);  // SiPad slab (used by ACTSProtoTracker for surfaceByAddress)
      mhit.setTime(hit.getTime());
      mhit.setEDep(hit.getEnergy());
      mhit.setEDepError(0.0f);
      mhit.setPosition(edm4hep::Vector3d{pos.x, pos.y, pos.z});
      mhit.setCovMatrix(cov);

      ++nConverted;
    }

    debug() << "[SiPadMeasConverter] evt=" << evtNum
            << " input=" << hits->size()
            << " converted=" << nConverted
            << " skipped=" << nSkipped
            << " vetoed(shower)=" << nVetoed << endmsg;

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[SiPadMeasConverter] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[SiPadMeasConverter] Unknown exception in execute()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode SiPadMeasConverter::finalize() {
  try {
    m_inputHandle.reset();
    m_outputHandle.reset();
    info() << "[SiPadMeasConverter] Done. Events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[SiPadMeasConverter] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[SiPadMeasConverter] Unknown exception in finalize()."
            << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(SiPadMeasConverter)
