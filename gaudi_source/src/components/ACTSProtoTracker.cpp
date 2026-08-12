// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4FWCore/DataHandle.h"

// edm4hep input
#include "edm4hep/TrackerHit3DCollection.h"

// edm4hep output
#include "edm4hep/TrackCollection.h"
#include "edm4hep/MutableTrack.h"

// ACTS
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"

// ACTS propagator
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"

// ACTS track fitting + finding
#include "Acts/TrackFitting/GainMatrixUpdater.hpp"
#include "Acts/TrackFitting/GainMatrixSmoother.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/TrackFinding/CombinatorialKalmanFilter.hpp"
#include "Acts/TrackFinding/MeasurementSelector.hpp"
#include "Acts/TrackFinding/TrackStateCreator.hpp"
#include "Acts/EventData/Types.hpp"

// ACTS track containers
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/TrackParameters.hpp"

// ACTS magnetic field
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

// ACTS calibration context
#include "Acts/Utilities/CalibrationContext.hpp"

// ACTS measurement helpers
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/SubspaceHelpers.hpp"
#include "Acts/EventData/SourceLink.hpp"

// SND geometry service
#include "k4SiWEcalReco/ISNDGeoSvc.h"

// Standard
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// ---------------------------------------------------------------------------
// SNDMeasurement
// ---------------------------------------------------------------------------

struct SNDMeasurement {
  const Acts::Surface* surface     = nullptr;
  double               localCoord  = 0.0;
  double               localCoord2 = 0.0;
  double               variance    = 0.0;
  double               variance2   = 0.0;
  bool                 is2D        = false;
  int                  detectorID  = -1;
  int                  plane       = -1;
  float                time        = 0.0f;
  float                eDep        = 0.0f;

  SNDMeasurement(const Acts::Surface* sf,
                 double lc, double lc2,
                 double var, double var2,
                 bool twod, int detID, int pl,
                 float t, float e)
      : surface(sf), localCoord(lc), localCoord2(lc2),
        variance(var), variance2(var2), is2D(twod),
        detectorID(detID), plane(pl), time(t), eDep(e) {}
};

// ---------------------------------------------------------------------------
// SNDSourceLink
// ---------------------------------------------------------------------------

struct SNDSourceLink {
  Acts::GeometryIdentifier geometryId() const { return m_geometryId; }
  std::size_t index = 0;
  void setGeometryId(Acts::GeometryIdentifier gid) { m_geometryId = gid; }
private:
  Acts::GeometryIdentifier m_geometryId;
};

// ---------------------------------------------------------------------------
// SNDSourceLinkAccessor
// Returns the range [begin, end) of SourceLinks associated with a surface.
// The slinks vector must be sorted by geometryId() before use.
// ---------------------------------------------------------------------------

struct SNDSourceLinkAccessor {
  const std::vector<Acts::SourceLink>* slinks = nullptr;
  std::pair<std::vector<Acts::SourceLink>::const_iterator,
            std::vector<Acts::SourceLink>::const_iterator>
  operator()(const Acts::Surface& surface) const {
    const auto geoId = surface.geometryId();
    // std::equal_range needs the comparator in both directions, so use
    // lower_bound + upper_bound with direction-specific lambdas instead.
    auto lo = std::lower_bound(
        slinks->begin(), slinks->end(), geoId,
        [](const Acts::SourceLink& sl, const Acts::GeometryIdentifier& id) {
          return sl.get<SNDSourceLink>().geometryId() < id;
        });
    auto hi = std::upper_bound(
        lo, slinks->end(), geoId,
        [](const Acts::GeometryIdentifier& id, const Acts::SourceLink& sl) {
          return id < sl.get<SNDSourceLink>().geometryId();
        });
    return {lo, hi};
  }
};

// ---------------------------------------------------------------------------
// SNDSurfaceAccessor — resolves a SourceLink back to its Acts::Surface.
// Required by KalmanFitterExtensions: the KF groups input source links by
// surface before propagation (SNDSourceLink only stores the geometryId).
// ---------------------------------------------------------------------------

struct SNDSurfaceAccessor {
  const std::vector<SNDMeasurement>* meas = nullptr;
  const Acts::Surface* operator()(const Acts::SourceLink& sl) const {
    return (*meas)[sl.get<SNDSourceLink>().index].surface;
  }
};

// ---------------------------------------------------------------------------
// SNDFixedNavigator — wraps DirectNavigator; injects the surface list from
// Config at makeState() time so the CKF's setPlainOptions() (which only
// copies NavigatorPlainOptions base fields) cannot accidentally erase it.
// ---------------------------------------------------------------------------

struct SNDFixedNavigator {
  // Surface sequence set once at construction; never touched by CKF options.
  std::vector<const Acts::Surface*> surfaces;

  // Options: minimal wrapper around NavigatorPlainOptions.
  // The CKF calls setPlainOptions() on this — it only copies base fields,
  // which is exactly what we want (surfaces stay in Config, not Options).
  struct Options : public Acts::NavigatorPlainOptions {
    explicit Options(const Acts::GeometryContext& gctx)
        : Acts::NavigatorPlainOptions(gctx) {}
    void setPlainOptions(const Acts::NavigatorPlainOptions& opts) {
      static_cast<Acts::NavigatorPlainOptions&>(*this) = opts;
    }
  };

  using State = Acts::DirectNavigator::State;

  // makeState: inject our surface list before creating the DirectNavigator state.
  State makeState(const Options& opts) const {
    // Diagnostic: confirm this navigator is actually being used.
    static std::atomic<int> s_callCount{0};
    if (s_callCount.fetch_add(1) < 3) {
      std::fprintf(stderr,
          "[SNDFixedNavigator] makeState called: %zu surfaces in Config\n",
          surfaces.size());
    }
    Acts::DirectNavigator::Options dirOpts(opts.geoContext);
    static_cast<Acts::NavigatorPlainOptions&>(dirOpts) =
        static_cast<const Acts::NavigatorPlainOptions&>(opts);
    dirOpts.surfaces = surfaces;
    Acts::DirectNavigator inner;
    return inner.makeState(dirOpts);
  }

  // Delegate all navigator interface methods to DirectNavigator.
  Acts::Result<void> initialize(State& state, const Acts::Vector3& pos,
                                const Acts::Vector3& dir,
                                Acts::Direction propDir) const {
    Acts::DirectNavigator inner;
    return inner.initialize(state, pos, dir, propDir);
  }
  Acts::NavigationTarget nextTarget(State& state, const Acts::Vector3& pos,
                                    const Acts::Vector3& dir) const {
    Acts::DirectNavigator inner;
    return inner.nextTarget(state, pos, dir);
  }
  bool checkTargetValid(const State& state, const Acts::Vector3& pos,
                        const Acts::Vector3& dir) const {
    Acts::DirectNavigator inner;
    return inner.checkTargetValid(state, pos, dir);
  }
  void handleSurfaceReached(State& state, const Acts::Vector3& pos,
                            const Acts::Vector3& dir,
                            const Acts::Surface& sf) const {
    Acts::DirectNavigator inner;
    inner.handleSurfaceReached(state, pos, dir, sf);
  }
  const Acts::Surface* currentSurface(const State& s) const {
    return s.currentSurface;
  }
  const Acts::TrackingVolume* currentVolume(const State&) const {
    return nullptr;
  }
  const Acts::IVolumeMaterial* currentVolumeMaterial(const State&) const {
    return nullptr;
  }
  const Acts::Surface* startSurface(const State& s) const {
    return s.options.startSurface;
  }
  const Acts::Surface* targetSurface(const State& s) const {
    return s.options.targetSurface;
  }
  // For a FIXED surface sequence, running out of surfaces IS the end of the
  // world, and it must be reported as such: Propagator::propagate() only
  // leaves its stepping loop when an aborter fires, and the CKF's only
  // geometric aborter is EndOfWorldReached (which calls exactly this method).
  // Returning a hardcoded false left the propagator free-stepping at
  // maxStepSize after DirectNavigator had set navigationBreak, until it hit
  // maxSteps and returned PropagatorError::StepCountLimitReached — the CKF
  // then discarded the branch, so EVERY event produced zero tracks.
  bool endOfWorldReached(State& s) const { return s.navigationBreak; }
  bool navigationBreak(const State& s) const { return s.navigationBreak; }
};

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

using SNDStepper       = Acts::EigenStepper<>;
using SNDCKFNavigator  = SNDFixedNavigator;
using SNDCKFPropagator = Acts::Propagator<SNDStepper, SNDCKFNavigator>;

using SNDTrackContainer = Acts::TrackContainer<
    Acts::VectorTrackContainer,
    Acts::VectorMultiTrajectory,
    std::shared_ptr>;

using SNDCKF = Acts::CombinatorialKalmanFilter<SNDCKFPropagator,
                                               SNDTrackContainer>;

// Final-refit KalmanFitter: the surface-sequence fit() overload requires a
// plain Acts::DirectNavigator (the sequence is passed per-fit call).
using SNDKFPropagator = Acts::Propagator<SNDStepper, Acts::DirectNavigator>;
using SNDKF           = Acts::KalmanFitter<SNDKFPropagator,
                                           Acts::VectorMultiTrajectory>;

// ---------------------------------------------------------------------------
// IronSlabBField — By inside registered slabs, zero everywhere else.
// Unlike Acts::MultiRangeBField, never returns MagneticFieldError for
// positions outside all registered ranges (returns zero instead).
// ---------------------------------------------------------------------------

class IronSlabBField : public Acts::MagneticFieldProvider {
public:
  struct Cache {
    explicit Cache(const Acts::MagneticFieldContext&) {}
  };
  struct Slab { double xlo, xhi, ylo, yhi, zlo, zhi, by; };

  explicit IronSlabBField(std::vector<Slab> slabs) : m_slabs(std::move(slabs)) {}

  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override {
    return Acts::MagneticFieldProvider::Cache(std::in_place_type<Cache>, mctx);
  }

  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& pos,
      Acts::MagneticFieldProvider::Cache&) const override {
    for (const auto& s : m_slabs) {
      if (pos.x() >= s.xlo && pos.x() <= s.xhi &&
          pos.y() >= s.ylo && pos.y() <= s.yhi &&
          pos.z() >= s.zlo && pos.z() <= s.zhi) {
        return Acts::Result<Acts::Vector3>::success(Acts::Vector3(0.0, s.by, 0.0));
      }
    }
    return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
  }

private:
  std::vector<Slab> m_slabs;
};

// ---------------------------------------------------------------------------
// ACTSProtoTracker
// ---------------------------------------------------------------------------

class ACTSProtoTracker : public Gaudi::Algorithm {
public:
  ACTSProtoTracker(const std::string& name, ISvcLocator* svcLoc)
      : Gaudi::Algorithm(name, svcLoc) {}

  StatusCode initialize() override;
  StatusCode execute(const EventContext&) const override;
  StatusCode finalize() override;

private:
  Gaudi::Property<std::string> m_inputSiPad{
      this, "InputSiPad", "SiPadMeasurements",
      "SiPad TrackerHit3DCollection from SiPadMeasConverter"};
  Gaudi::Property<std::string> m_outputCollection{
      this, "OutputCollection", "ACTSTracks",
      "Output edm4hep::TrackCollection name"};
  Gaudi::Property<double> m_bFieldX{this, "BFieldX", 0.0, "BField X [T]"};
  Gaudi::Property<double> m_bFieldY{this, "BFieldY", 0.0, "BField Y [T]"};
  Gaudi::Property<double> m_bFieldZ{this, "BFieldZ", 0.0, "BField Z [T]"};
  Gaudi::Property<std::vector<double>> m_ironFieldRanges{
      this, "IronFieldRanges", {},
      "Per-slab field: [xlo,xhi, ylo,yhi, zlo,zhi, by] x N slabs in ACTS coords [mm, T]"};

  // ---- Seed configuration (DD4hep convention: Z=beam) ----
  // SeedPositions: flat list of (x, y, z) triplets in mm.
  //   x = transverse X [mm]
  //   y = transverse Y [mm]
  //   z = beam position [mm] — used only to select starting surface
  // Default: one seed at the beam center.
  Gaudi::Property<std::vector<double>> m_seedPositions{
      this, "SeedPositions", {0.0, 0.0, 0.0},
      "Seed positions in DD4hep convention (Z=beam) as flat (x,y,z) triplets [mm]."};

  // SeedDirections: flat list of (dx, dy, dz) unit vectors.
  //   dz = beam direction component (dominant for forward tracks)
  // The algorithm internally swaps to ACTS convention before fitting.
  // Default: one seed along +Z (beam direction).
  Gaudi::Property<std::vector<double>> m_seedDirections{
      this, "SeedDirections", {0.0, 0.0, 1.0},
      "Seed directions in DD4hep convention (Z=beam) as flat (dx,dy,dz) triplets."};

  // Seed momentum magnitude [GeV], same for all seeds.
  Gaudi::Property<double> m_seedMomentum{
      this, "SeedMomentum", 10.0,
      "Seed momentum magnitude [GeV]."};

  // Enable automatic Hough Transform seeding (ignores SeedPositions/SeedDirections).
  Gaudi::Property<bool> m_autoSeed{
      this, "AutoSeed", false,
      "If true, use Hough Transform automatic seeding instead of manual seeds."};

  // Maximum number of seeds from auto-seeding.
  Gaudi::Property<int> m_maxSeeds{
      this, "MaxSeeds", 5,
      "Maximum number of seeds returned by auto-seeding."};

  // Hough Transform histogram bin size [mm]
  Gaudi::Property<double> m_houghBinSize{
      this, "HoughBinSize", 5.0,
      "Bin size for Hough Transform 2D histogram [mm]. "
      "Should be ~1 SiPad pixel size."};

  // Hough Transform detector half-size [mm] — sets histogram range
  Gaudi::Property<double> m_houghHalfSize{
      this, "HoughHalfSize", 200.0,
      "Half-size of transverse detector for Hough histogram [mm]."};

  // Minimum number of votes for a Hough peak to be considered a seed
  Gaudi::Property<int> m_houghMinVotes{
      this, "HoughMinVotes", 3,
      "Minimum number of hits voting for a Hough peak to create a seed."};

  // Compatibility radius: hit is compatible with seed if within this [mm]
  Gaudi::Property<double> m_seedCompatRadius{
      this, "SeedCompatRadius", 10.0,
      "Radius [mm] within which a hit is considered compatible with a seed."};

  // Strip pitch used for seed position refinement [mm]
  Gaudi::Property<double> m_stripPitch{
      this, "SeedStripPitch", 0.0755,
      "Strip pitch [mm] used for most-frequent-strip seed refinement."};

  Gaudi::Property<double> m_maxChi2PerNdf{
    this, "MaxChi2PerNdf", 500.0,
    "Maximum chi2/ndf threshold for track acceptance, where "
    "ndf = Σ calibratedSize − 5 helix params. Tracks above this "
    "threshold are rejected as false seeds (≈1 = ideal, >>1 = bad fit)."};

  Gaudi::Property<double> m_chi2CutOff{
      this, "Chi2CutOff", 15.0,
      "chi2 cut for MeasurementSelector: maximum local chi2 to accept a "
      "measurement on a surface during CKF track finding."};

  Gaudi::Property<int> m_numMeasCutOff{
      this, "NumMeasCutOff", 1,
      "Maximum number of measurements accepted per surface by "
      "MeasurementSelector during CKF track finding."};

  Gaudi::Property<int> m_maxPropSteps{
      this, "MaxPropSteps", 100000,
      "Maximum number of propagation steps allowed by the CKF EigenStepper. "
      "Increase if PropagatorError:2 occurs (step-count limit hit)."};

  Gaudi::Property<double> m_maxStepSize{
      this, "MaxStepSize", 100.0,
      "Maximum single-step size [mm] for the CKF EigenStepper. "
      "100 mm gives ~0.06 mm position error per step at 1.7 T / 10 GeV."};

  // Track/shower classification: maximum crossing multiplicity per station.
  // For each Hough peak, multiplicity = nCompatCrossings / nCompatStations.
  //   Track:  1 StripX x 1 StripY per station → multiplicity ≈ 1-4
  //   Shower: many strips per station          → multiplicity ≈ 100+
  // Set to 1e9 to disable (default: disabled).
  // Tune by running with OutputLevel=DEBUG and observing "multiplicity=" values.
  Gaudi::Property<double> m_houghMaxMultiplicity{
      this, "HoughMaxMultiplicity", 1e9,
      "Maximum crossing multiplicity per station for a Hough peak to be "
      "considered a track candidate. Peaks above this are classified as "
      "showers and skipped. Set to 1e9 to disable (default: disabled)." };

  // Strip isolation filter: window size [mm] for counting strip neighbors.
  // A strip is considered isolated if fewer than IsolationMaxNeighbors
  // other strips of the same type (StripX or StripY) fall within this
  // window in the same station. Set to 0.0 to disable (default: disabled).
  Gaudi::Property<double> m_isolationWindow{
      this, "IsolationWindow", 0.0,
      "2D distance [mm] for crossing isolation filter. A crossing is kept "
      "only if fewer than IsolationMaxNeighbors OTHER crossings in the same "
      "station fall within this 2D distance. Set to 0.0 to disable."};

  // Maximum number of crossing neighbors within IsolationWindow to be
  // considered isolated (track-like).
  // Muon:  0 neighbors (exactly 1 crossing per station) → passes filter
  // Shower: hundreds of neighbors → fails filter → discarded
  Gaudi::Property<int> m_isolationMaxNeighbors{
      this, "IsolationMaxNeighbors", 2,
      "Maximum number of crossing neighbors within IsolationWindow for a "
      "crossing to be considered isolated (track-like). "
      "Muon: 0 neighbors. Shower: hundreds of neighbors."};

  Gaudi::Property<bool> m_finalRefit{
      this, "FinalRefit", true,
      "If true, refit the hits selected by the CKF with an Acts::KalmanFitter "
      "on the frozen hit set (no re-selection possible), giving unbiased "
      "parameters and covariances. Falls back to the CKF track on failure."};

  Gaudi::Property<bool> m_seedCleaning{
      this, "SeedCleaning", true,
      "If true, the measurements assigned to an accepted track are removed "
      "from the source-link pool before processing the next seed. Nearby "
      "seeds then cannot converge onto an already-used particle, which "
      "recovers the second track at small opening angles."};

  // ---- Shower-hit purge (measurement-pool level) ---------------------------
  Gaudi::Property<double> m_hitPurgeWindow{
      this, "HitPurgeWindow", 0.0,
      "Distance [mm] for the density-based shower-hit purge of the "
      "measurement pool (applied per surface, before seeding and CKF). A "
      "measurement with more than HitPurgeMaxNeighbors same-surface "
      "neighbors within this distance is dropped. 0.0 = disabled."};

  Gaudi::Property<int> m_hitPurgeMaxNeighbors{
      this, "HitPurgeMaxNeighbors", 4,
      "Maximum number of same-surface neighbors within HitPurgeWindow for a "
      "measurement to be kept. MIP + delta rays: a few. Shower core: tens."};

  // ---- Duplicate-track filter ---------------------------------------------
  // Two accepted tracks are duplicates if they share more than
  // DuplicateOverlapFraction of the SMALLER track's hits. The event's tracks
  // are deduplicated at the END of the event, best-first (lowest chi2/ndf), so
  // the survivor of each overlap group is always the best-quality fit —
  // independent of the order seeds were processed in.
  Gaudi::Property<double> m_dupOverlapFraction{
      this, "DuplicateOverlapFraction", 0.7,
      "Fraction of the smaller track's hits that must be shared for two tracks "
      "to be considered duplicates (0..1). The higher-chi2/ndf one is dropped."};

  // ---- Hough-based auto-seeding -------------------------------------------
  struct SeedCandidate {
    double x;            // transverse X [mm], extrapolated to the first surface
    double y;            // transverse Y [mm], extrapolated to the first surface
    double dx;           // dx/dz slope [-] from the straight-line regression
    double dy;           // dy/dz slope [-] from the straight-line regression
    double z_start;      // beam Z of first compatible layer [mm]
    int    nVotes;       // number of crossing points supporting this seed
    double multiplicity; // nCompatCrossings / nCompatStations — track/shower discriminant
  };

  std::vector<SeedCandidate> findSeeds(
      const std::vector<SNDMeasurement>& measurements,
      const Acts::GeometryContext& gctx) const;

  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>
      m_siPadHandle;
  mutable std::unique_ptr<k4FWCore::DataHandle<edm4hep::TrackCollection>>
      m_outputHandle;

  ServiceHandle<ISNDGeoSvc> m_geoSvc{
      this, "GeoSvc", "ACTSGeoSvc", "ACTS geometry service"};

  mutable std::atomic<long long> m_eventCount{0};
  mutable Acts::MagneticFieldContext m_mctx;
  mutable Acts::CalibrationContext   m_cctx;


};

// ---------------------------------------------------------------------------
// findSeeds() — Hough Transform automatic seeder
// ---------------------------------------------------------------------------

std::vector<ACTSProtoTracker::SeedCandidate> ACTSProtoTracker::findSeeds(
    const std::vector<SNDMeasurement>& measurements,
    const Acts::GeometryContext& gctx) const
{
  std::vector<SeedCandidate> seeds;
  if (measurements.empty()) return seeds;

  // =========================================================================
  // STEP 1: Build list of unambiguous 2D points from SiPad hits
  // =========================================================================

  struct Point2D {
    double x;      // transverse X [mm]
    double y;      // transverse Y [mm]
    double z;      // beam Z [mm] (ACTS X coordinate of surface)
    int    weight; // 1 for SiPad
  };
  std::vector<Point2D> points2D;

  // Isolation filter settings
  const double isolWin      = m_isolationWindow.value();
  const int    isolMaxNeigh = m_isolationMaxNeighbors.value();
  const bool   doIsolation  = (isolWin > 0.0);
  const double isolWin2     = isolWin * isolWin;

  // SiPad 2D hits with position-level isolation per layer.
  std::map<const Acts::Surface*, std::vector<const SNDMeasurement*>> sipadBySurface;
  for (const auto& m : measurements) {
    if (m.is2D) sipadBySurface[m.surface].push_back(&m);
  }

  for (const auto& [surf, layerHits] : sipadBySurface) {
    double beamZ = surf->center(gctx).x();
    for (std::size_t i = 0; i < layerHits.size(); ++i) {
      const auto* mp = layerHits[i];

      if (doIsolation) {
        int nNeigh = 0;
        for (std::size_t j = 0; j < layerHits.size(); ++j) {
          if (j == i) continue;
          double dx = layerHits[j]->localCoord  - mp->localCoord;
          double dy = layerHits[j]->localCoord2 - mp->localCoord2;
          if (dx*dx + dy*dy < isolWin2) ++nNeigh;
          if (nNeigh > isolMaxNeigh) break;
        }
        if (nNeigh > isolMaxNeigh) continue;
      }

      points2D.push_back({mp->localCoord, mp->localCoord2, beamZ, 1});
    }
  }

  if (points2D.empty()) {
    SeedCandidate sc{};
    sc.x = 0.0; sc.y = 0.0; sc.dx = 0.0; sc.dy = 0.0; sc.nVotes = 0;
    sc.z_start = measurements.empty() ? 0.0
               : measurements.front().surface->center(gctx).x();
    seeds.push_back(sc);
    return seeds;
  }

  // =========================================================================
  // STEP 2: Hough Transform on 2D points only
  // =========================================================================
  // Each point2D votes for one (x,y) bin — no strip ambiguity!

  const double halfSize = m_houghHalfSize.value();
  const double binSize  = m_houghBinSize.value();
  const int    nBins    = static_cast<int>(2.0 * halfSize / binSize) + 1;

  std::vector<std::vector<int>> histo(nBins, std::vector<int>(nBins, 0));

  auto toBin = [&](double coord) -> int {
    int bin = static_cast<int>((coord + halfSize) / binSize);
    return std::max(0, std::min(nBins - 1, bin));
  };
  auto fromBin = [&](int bin) -> double {
    return -halfSize + (bin + 0.5) * binSize;
  };

  for (const auto& p : points2D) {
    int ix = toBin(p.x);
    int iy = toBin(p.y);
    histo[ix][iy] += p.weight;
  }

  // =========================================================================
  // STEP 3: Find local maxima
  // =========================================================================
  const int minVotes    = m_houghMinVotes.value();
  const int suppressRad = static_cast<int>(std::ceil(15.0 / binSize));

  struct Peak { int ix; int iy; int votes; };
  std::vector<Peak> peaks;

  for (int ix = 0; ix < nBins; ++ix) {
    for (int iy = 0; iy < nBins; ++iy) {
      int v = histo[ix][iy];
      if (v < minVotes) continue;
      bool isMax = true;
      for (int dx = -1; dx <= 1 && isMax; ++dx) {
        for (int dy = -1; dy <= 1 && isMax; ++dy) {
          if (dx == 0 && dy == 0) continue;
          int nx = ix + dx, ny = iy + dy;
          if (nx < 0 || nx >= nBins || ny < 0 || ny >= nBins) continue;
          if (histo[nx][ny] > v) isMax = false;
        }
      }
      if (isMax) peaks.push_back({ix, iy, v});
    }
  }

  std::sort(peaks.begin(), peaks.end(),
            [](const Peak& a, const Peak& b) { return a.votes > b.votes; });

  // =========================================================================
  // STEP 3b: Compute crossing multiplicity per station for each peak
  // =========================================================================
  // multiplicity = nCompatCrossings / nCompatStations
  // Track:  1 crossing per station → multiplicity ≈ 1
  // Shower: many crossings per station → multiplicity ≈ 100+
  // A Point2D from a crossing has a z coordinate (station beam Z).
  // Round z to the nearest stationTolerance to identify unique stations.

  const double compatR         = m_seedCompatRadius.value();
  const double compatR2        = compatR * compatR;
  const double maxMult         = m_houghMaxMultiplicity.value();
  // Layer pitch ~11 mm; use half-pitch so each SiPad layer maps to one station key.
  const double stationTolerance = 6.0;  // mm

  struct PeakWithMult {
    int    ix, iy, votes;
    double multiplicity;
  };
  std::vector<PeakWithMult> peaksWithMult;

  for (const auto& pk : peaks) {
    const double peakX = fromBin(pk.ix);
    const double peakY = fromBin(pk.iy);

    int nCompatCrossings = 0;
    std::set<int> compatStations;  // unique station keys within compatR

    for (const auto& p : points2D) {
      double dx = p.x - peakX, dy = p.y - peakY;
      if (dx*dx + dy*dy < compatR2) {
        nCompatCrossings += p.weight;
        // Station key: round beam Z to nearest stationTolerance mm
        int stationKey = static_cast<int>(
            std::round(p.z / stationTolerance));
        compatStations.insert(stationKey);
      }
    }

    const int nStations = static_cast<int>(compatStations.size());
    const double mult   = (nStations > 0)
        ? static_cast<double>(nCompatCrossings) / nStations
        : 0.0;

    peaksWithMult.push_back({pk.ix, pk.iy, pk.votes, mult});
  }

  // Non-maximum suppression (same as before, operating on peaksWithMult)
  std::vector<bool> suppressed(peaksWithMult.size(), false);
  for (std::size_t i = 0; i < peaksWithMult.size(); ++i) {
    if (suppressed[i]) continue;
    for (std::size_t j = i + 1; j < peaksWithMult.size(); ++j) {
      if (suppressed[j]) continue;
      int dx = std::abs(peaksWithMult[i].ix - peaksWithMult[j].ix);
      int dy = std::abs(peaksWithMult[i].iy - peaksWithMult[j].iy);
      if (dx <= suppressRad && dy <= suppressRad) suppressed[j] = true;
    }
  }

  // =========================================================================
  // STEP 4: Classify peaks and refine seed position using most-frequent strip
  // =========================================================================
  const int    maxS       = m_maxSeeds.value();
  const double stripPitch = m_stripPitch.value();
  // Beam coordinate of the surface the seed is built on (execute() always uses
  // allSurfaces.front()), so the regression extrapolates to the right plane.
  const auto&  geoSurfaces = m_geoSvc->allSurfaces();
  const double zSeed = geoSurfaces.empty()
      ? measurements.front().surface->center(gctx).x()
      : geoSurfaces.front()->center(gctx).x();

  for (std::size_t pi = 0; pi < peaksWithMult.size() && (int)seeds.size() < maxS; ++pi) {
    if (suppressed[pi]) continue;

    const auto& pk     = peaksWithMult[pi];
    const double peakX = fromBin(pk.ix);
    const double peakY = fromBin(pk.iy);

    // Log all peaks at DEBUG level for tuning HoughMaxMultiplicity
    debug() << "[ACTSProtoTracker] Hough peak: x=" << peakX
            << " y=" << peakY
            << " votes=" << pk.votes
            << " multiplicity=" << pk.multiplicity
            << " (crossings/station)" << endmsg;

    // Track/shower classification: skip high-multiplicity peaks
    if (pk.multiplicity > maxMult) {
      debug() << "[ACTSProtoTracker]   -> SHOWER (multiplicity=" << pk.multiplicity
              << " > HoughMaxMultiplicity=" << maxMult << ") -- skipped." << endmsg;
      continue;
    }

    // Straight-line regression over the compatible points, plus the
    // most-frequent-position fallback.
    //
    // The Hough peak is only the (x,y) CENTROID over the layers: seeding with
    // it and a direction fixed along the beam is right only for a track
    // exactly parallel to the beam. For a track crossing pads (a tilt of just
    // atan(5.53/15) ~ 20 deg is a full pad per layer, and beam divergence
    // alone moves the entry point over 225 mm) the prediction on the first
    // layer lands pads away from the hit, the MeasurementSelector rejects it,
    // and the trajectory is never corrected — the CKF then picks up one or two
    // hits and the candidate is discarded. Fitting x(z) and y(z) gives the
    // seed both the correct entry point AND a direction estimate.
    std::map<int,int> xFreq, yFreq;
    double firstZ = std::numeric_limits<double>::max();
    int nPts = 0;
    // Regression accumulators (z = beam coordinate).
    double sw = 0.0, sz = 0.0, szz = 0.0;
    double sx = 0.0, szx = 0.0, sy = 0.0, szy = 0.0;
    std::set<int> zKeys;

    for (const auto& p : points2D) {
      double ddx = p.x - peakX, ddy = p.y - peakY;
      if (ddx*ddx + ddy*ddy < compatR2) {
        xFreq[static_cast<int>(std::round(p.x/stripPitch))] += p.weight;
        yFreq[static_cast<int>(std::round(p.y/stripPitch))] += p.weight;
        nPts += p.weight;
        if (p.z < firstZ) firstZ = p.z;
        const double w = static_cast<double>(p.weight);
        sw  += w;        sz  += w * p.z;   szz += w * p.z * p.z;
        sx  += w * p.x;  szx += w * p.z * p.x;
        sy  += w * p.y;  szy += w * p.z * p.y;
        zKeys.insert(static_cast<int>(std::round(p.z / stationTolerance)));
      }
    }

    if (nPts == 0) continue;

    double refinedX = peakX;
    int maxXF = 0;
    for (const auto& [strip, freq] : xFreq)
      if (freq > maxXF) { maxXF = freq; refinedX = strip * stripPitch; }

    double refinedY = peakY;
    int maxYF = 0;
    for (const auto& [strip, freq] : yFreq)
      if (freq > maxYF) { maxYF = freq; refinedY = strip * stripPitch; }

    // The seed is built on the first surface, so extrapolate there. Require
    // at least 3 distinct layers and a non-degenerate z lever arm; otherwise
    // fall back to the parallel-to-beam assumption.
    const double det = sw * szz - sz * sz;
    double slopeX = 0.0, slopeY = 0.0;
    if (zKeys.size() >= 3 && std::abs(det) > 1e-6) {
      slopeX = (sw * szx - sz * sx) / det;
      slopeY = (sw * szy - sz * sy) / det;
      const double interX = (sx - slopeX * sz) / sw;
      const double interY = (sy - slopeY * sz) / sw;
      refinedX = interX + slopeX * zSeed;
      refinedY = interY + slopeY * zSeed;
    }

    SeedCandidate sc{};
    sc.x            = refinedX;
    sc.y            = refinedY;
    sc.dx           = slopeX;
    sc.dy           = slopeY;
    sc.z_start      = firstZ;
    sc.nVotes       = pk.votes;
    sc.multiplicity = pk.multiplicity;
    seeds.push_back(sc);
  }

  return seeds;
}

// ---------------------------------------------------------------------------
// initialize()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::initialize() {
  try {
    StatusCode sc = Gaudi::Algorithm::initialize();
    if (sc.isFailure()) return sc;

    m_siPadHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackerHit3DCollection>>(
        m_inputSiPad.value(), Gaudi::DataHandle::Reader, this);
    m_outputHandle = std::make_unique<
        k4FWCore::DataHandle<edm4hep::TrackCollection>>(
        m_outputCollection.value(), Gaudi::DataHandle::Writer, this);

    if (!m_geoSvc.retrieve().isSuccess()) {
      error() << "[ACTSProtoTracker] Failed to retrieve ACTSGeoSvc." << endmsg;
      return StatusCode::FAILURE;
    }

    // Verify surfaces are available.
    const auto& allSurf = m_geoSvc->allSurfaces();
    const auto& gctx    = m_geoSvc->geometryContext();

    if (allSurf.empty()) {
      error() << "[ACTSProtoTracker] No surfaces in geometry." << endmsg;
      return StatusCode::FAILURE;
    }

    info() << "[ACTSProtoTracker] Initialized. GeoSvc has "
           << allSurf.size() << " SiPad surfaces." << endmsg;

    // Validate seed properties
    const std::size_t nPosVals = m_seedPositions.value().size();
    const std::size_t nDirVals = m_seedDirections.value().size();
    if (nPosVals % 3 != 0 || nDirVals % 3 != 0) {
      error() << "[ACTSProtoTracker] SeedPositions and SeedDirections must "
              << "have sizes that are multiples of 3." << endmsg;
      return StatusCode::FAILURE;
    }
    if (nPosVals / 3 != nDirVals / 3) {
      error() << "[ACTSProtoTracker] SeedPositions and SeedDirections must "
              << "have the same number of (x,y,z) triplets." << endmsg;
      return StatusCode::FAILURE;
    }
    info() << "[ACTSProtoTracker] Configured with "
           << nPosVals / 3 << " seed(s)." << endmsg;

    return sc;
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in initialize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in initialize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// execute()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::execute(const EventContext&) const {
  try {
    const long long evtNum = m_eventCount.fetch_add(1);

    auto* output = m_outputHandle->createAndPut();

    const auto& allSurfaces = m_geoSvc->allSurfaces();
    const auto& gctx        = m_geoSvc->geometryContext();

    // =========================================================================
    // STEP 1: Collect measurements
    // =========================================================================
    std::vector<SNDMeasurement> measurements;

    // ---- SiPad (2D pixels) --------------------------------------------------
    // Surface lookup is address-based via ISNDGeoSvc::surfaceByAddress (detID=1,
    // station=-1, layer=quality, plane=-1). The layer index is written by
    // SiPadMeasConverter into TrackerHit3D::quality. This replaces the previous
    // nearest-z search, which depended on a fragile largest-gap heuristic for the
    // Surface lookup is address-based via ISNDGeoSvc::surfaceByAddress.
    const auto* spHits = m_siPadHandle->get();
    if (spHits) {
      for (std::size_t i = 0; i < spHits->size(); ++i) {
        const auto& hit   = (*spHits)[i];
        const auto& pos   = hit.getPosition();
        const int   layer = hit.getQuality();

        const Acts::Surface* surf =
            m_geoSvc->surfaceByAddress(1, -1, layer, -1);
        if (!surf) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " SiPad hit " << i << " layer=" << layer
                    << " has no matching surface." << endmsg;
          continue;
        }

        const auto& cov = hit.getCovMatrix();
        measurements.emplace_back(surf, pos.x, pos.y,
                                   cov[0], cov[3], true, 1, -1,
                                   hit.getTime(), hit.getEDep());
      }
    }

    if (measurements.empty()) {
      debug() << "[ACTSProtoTracker] evt=" << evtNum
              << " no measurements." << endmsg;
      return StatusCode::SUCCESS;
    }

    std::sort(measurements.begin(), measurements.end(),
              [&](const SNDMeasurement& a, const SNDMeasurement& b) {
                return a.surface->center(gctx).x() <
                       b.surface->center(gctx).x();
              });

    // ---- Shower-hit purge ---------------------------------------------------
    // Density-based cleaning of the measurement pool BEFORE seeding and the
    // CKF: on each surface, a measurement with more than HitPurgeMaxNeighbors
    // other measurements within HitPurgeWindow is shower-like and is dropped
    // entirely (a MIP leaves 1-2 hits per layer; a shower leaves tens to
    // hundreds). Unlike the Hough isolation filter (which only steers the
    // seeding), this removes the hits from the fit itself, so shower cores
    // can no longer contaminate track candidates.
    const std::size_t nBeforePurge = measurements.size();
    const double purgeWin  = m_hitPurgeWindow.value();
    const int    purgeMaxN = m_hitPurgeMaxNeighbors.value();
    if (purgeWin > 0.0) {
      std::unordered_map<const Acts::Surface*, std::vector<std::size_t>> bySurf;
      for (std::size_t i = 0; i < measurements.size(); ++i) {
        bySurf[measurements[i].surface].push_back(i);
      }

      std::vector<bool> keep(measurements.size(), true);
      for (const auto& [surf, idxs] : bySurf) {
        if (static_cast<int>(idxs.size()) <= purgeMaxN) continue;  // cannot exceed
        for (std::size_t i : idxs) {
          const auto& a = measurements[i];
          int n = 0;
          for (std::size_t j : idxs) {
            if (j == i) continue;
            const auto& b = measurements[j];
            const double d = std::hypot(a.localCoord  - b.localCoord,
                                        a.localCoord2 - b.localCoord2);
            if (d < purgeWin && ++n > purgeMaxN) break;
          }
          if (n > purgeMaxN) keep[i] = false;
        }
      }

      std::vector<SNDMeasurement> purged;
      purged.reserve(measurements.size());
      for (std::size_t i = 0; i < measurements.size(); ++i) {
        if (keep[i]) purged.push_back(measurements[i]);
      }
      measurements.swap(purged);

      if (measurements.empty()) {
        debug() << "[ACTSProtoTracker] evt=" << evtNum
                << " all measurements purged as shower-like." << endmsg;
        return StatusCode::SUCCESS;
      }
    }

    debug() << "[ACTSProtoTracker] evt=" << evtNum
            << " SiPad=" << (spHits ? spHits->size() : 0)
            << " total measurements=" << measurements.size()
            << " (purged " << (nBeforePurge - measurements.size())
            << " shower-like)" << endmsg;

    // =========================================================================
    // STEP 4: Build shared KF components (once per event)
    // =========================================================================

    std::shared_ptr<Acts::MagneticFieldProvider> bField;
    const auto& ironRanges = m_ironFieldRanges.value();
    if (!ironRanges.empty() && ironRanges.size() % 7 == 0) {
      // Build IronSlabBField: each entry covers one outer iron slab.
      // Format per entry: [xlo,xhi, ylo,yhi, zlo,zhi, by] in [mm, T].
      // ACTS coords: x = beam axis (= DD4hep Z), y = DD4hep Y, z = DD4hep X.
      std::vector<IronSlabBField::Slab> slabs;
      slabs.reserve(ironRanges.size() / 7);
      for (std::size_t i = 0; i < ironRanges.size(); i += 7) {
        slabs.push_back({
          ironRanges[i+0] * Acts::UnitConstants::mm,
          ironRanges[i+1] * Acts::UnitConstants::mm,
          ironRanges[i+2] * Acts::UnitConstants::mm,
          ironRanges[i+3] * Acts::UnitConstants::mm,
          ironRanges[i+4] * Acts::UnitConstants::mm,
          ironRanges[i+5] * Acts::UnitConstants::mm,
          ironRanges[i+6] * Acts::UnitConstants::T
        });
      }
      bField = std::make_shared<IronSlabBField>(std::move(slabs));
    } else {
      bField = std::make_shared<Acts::ConstantBField>(
          Acts::Vector3(m_bFieldX.value() * Acts::UnitConstants::T,
                        m_bFieldY.value() * Acts::UnitConstants::T,
                        m_bFieldZ.value() * Acts::UnitConstants::T));
    }

    // ---- Calibrator ---------------------------------------------------------
    struct SNDCalibrator {
      const std::vector<SNDMeasurement>* meas = nullptr;
      void operator()(const Acts::GeometryContext& /*gctx*/,
                      const Acts::CalibrationContext& /*cctx*/,
                      const Acts::SourceLink& sl,
                      Acts::VectorMultiTrajectory::TrackStateProxy ts) const {
        const auto& ssl = sl.get<SNDSourceLink>();
        const auto& m   = (*meas)[ssl.index];
        // Store the source link on the track state: TrackStateCreator leaves
        // this to the calibrator (ACTS convention). Without it the hit
        // fingerprints (duplicate rejection, frozen-hit refit) are empty.
        ts.setUncalibratedSourceLink(Acts::SourceLink{sl});
        // SiPad: 2D measurement (loc0,loc1) = (DD4hep X, DD4hep Y)
        constexpr std::array<Acts::BoundIndices, 2> indices = {
            Acts::eBoundLoc0, Acts::eBoundLoc1};
        ts.allocateCalibrated(2);
        ts.template calibrated<2>() =
            Acts::ActsVector<2>(m.localCoord, m.localCoord2);
        ts.template calibratedCovariance<2>() =
            Acts::ActsSquareMatrix<2>{{m.variance, 0.0},
                                      {0.0, m.variance2}};
        ts.setProjectorSubspaceIndices(indices);
      }
    };

    SNDCalibrator calibrator;
    calibrator.meas = &measurements;

    // ---- Build sorted SourceLink list (CKF accessor needs sorted-by-geoId) --
    // One SourceLink per measurement; includes ALL measurements (not pre-filtered
    // per seed — the CKF's MeasurementSelector handles per-surface selection).
    std::vector<Acts::SourceLink> sortedSLinks;
    sortedSLinks.reserve(measurements.size());
    for (std::size_t i = 0; i < measurements.size(); ++i) {
      SNDSourceLink ssl;
      ssl.index = i;
      ssl.setGeometryId(measurements[i].surface->geometryId());
      sortedSLinks.push_back(Acts::SourceLink(ssl));
    }
    std::sort(sortedSLinks.begin(), sortedSLinks.end(),
              [](const Acts::SourceLink& a, const Acts::SourceLink& b) {
                return a.get<SNDSourceLink>().geometryId() <
                       b.get<SNDSourceLink>().geometryId();
              });

    // ---- CKF infrastructure (built once per event, reused across seeds) -----
    Acts::GainMatrixUpdater gainMatrixUpdater;
    SNDSourceLinkAccessor slAccessor;
    slAccessor.slinks = &sortedSLinks;

    Acts::MeasurementSelectorCuts measCuts;
    measCuts.chi2CutOff            = {m_chi2CutOff.value()};
    measCuts.numMeasurementsCutOff = {static_cast<std::size_t>(m_numMeasCutOff.value())};
    Acts::MeasurementSelector measSelector(measCuts);

    using SLinkIter = std::vector<Acts::SourceLink>::const_iterator;
    Acts::TrackStateCreator<SLinkIter, SNDTrackContainer> tsc;
    tsc.sourceLinkAccessor
        .connect<&SNDSourceLinkAccessor::operator()>(&slAccessor);
    tsc.calibrator
        .template connect<&SNDCalibrator::operator()>(&calibrator);
    tsc.measurementSelector
        .connect<&Acts::MeasurementSelector::select<Acts::VectorMultiTrajectory>>(
            &measSelector);

    Acts::CombinatorialKalmanFilterExtensions<SNDTrackContainer> ckfExtensions;
    ckfExtensions.updater
        .connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
            &gainMatrixUpdater);
    ckfExtensions.createTrackStates
        .connect<&Acts::TrackStateCreator<SLinkIter, SNDTrackContainer>::createTrackStates>(
            &tsc);

    // SNDFixedNavigator: surfaces injected from allSurfaces at makeState()
    // time. The CKF's setPlainOptions() only copies NavigatorPlainOptions base
    // fields and cannot erase the surface list stored in the Config.
    SNDFixedNavigator ckfNavigator;
    ckfNavigator.surfaces.assign(allSurfaces.begin(), allSurfaces.end());

    if (evtNum < 1) {
      debug() << "[ACTSProtoTracker] DIAG: ckfNavigator has "
              << ckfNavigator.surfaces.size() << " surfaces" << endmsg;
      if (!ckfNavigator.surfaces.empty()) {
        debug() << "[ACTSProtoTracker] DIAG: first surf geoId="
                << ckfNavigator.surfaces.front()->geometryId() << endmsg;
        debug() << "[ACTSProtoTracker] DIAG: last surf geoId="
                << ckfNavigator.surfaces.back()->geometryId() << endmsg;
      }
    }

    SNDStepper       ckfStepper(bField);
    SNDCKFPropagator ckfPropagator(std::move(ckfStepper), std::move(ckfNavigator));
    SNDCKF ckf(std::move(ckfPropagator),
               Acts::getDefaultLogger("CKF", Acts::Logging::WARNING));

    // ---- KalmanFitter infrastructure for the final refit --------------------
    // Reuses the same updater and calibrator as the CKF; the smoother and the
    // surface accessor are KF-specific. The surface sequence is passed per
    // fit() call, which is why a plain DirectNavigator is used.
    Acts::GainMatrixSmoother kfSmoother;
    SNDSurfaceAccessor kfSurfaceAccessor;
    kfSurfaceAccessor.meas = &measurements;

    Acts::KalmanFitterExtensions<Acts::VectorMultiTrajectory> kfExtensions;
    kfExtensions.calibrator
        .template connect<&SNDCalibrator::operator()>(&calibrator);
    kfExtensions.updater
        .connect<&Acts::GainMatrixUpdater::operator()<Acts::VectorMultiTrajectory>>(
            &gainMatrixUpdater);
    kfExtensions.smoother
        .connect<&Acts::GainMatrixSmoother::operator()<Acts::VectorMultiTrajectory>>(
            &kfSmoother);
    kfExtensions.surfaceAccessor
        .connect<&SNDSurfaceAccessor::operator()>(&kfSurfaceAccessor);

    SNDStepper kfStepper(bField);
    Acts::DirectNavigator kfNavigator;
    SNDKFPropagator kfPropagator(std::move(kfStepper), std::move(kfNavigator));
    SNDKF kf(std::move(kfPropagator),
             Acts::getDefaultLogger("KF", Acts::Logging::WARNING));

    // Track container accumulates tracks from all seeds this event.
    auto trackBackend = std::make_shared<Acts::VectorTrackContainer>();
    auto trajBackend  = std::make_shared<Acts::VectorMultiTrajectory>();
    SNDTrackContainer ckfTracks(trackBackend, trajBackend);

    // Shared seed covariance (loose — same for all seeds)
    Acts::BoundSquareMatrix seedCov = Acts::BoundSquareMatrix::Zero();
    seedCov(Acts::eBoundLoc0,   Acts::eBoundLoc0)   = 1e2;
    seedCov(Acts::eBoundLoc1,   Acts::eBoundLoc1)   = 1e2;
    seedCov(Acts::eBoundPhi,    Acts::eBoundPhi)     = 1.0;
    seedCov(Acts::eBoundTheta,  Acts::eBoundTheta)   = 1.0;
    seedCov(Acts::eBoundQOverP, Acts::eBoundQOverP)  = 0.04;
    seedCov(Acts::eBoundTime,   Acts::eBoundTime)    = 1e9;

    const double seedQoverP =
        -1.0 / (m_seedMomentum.value() * Acts::UnitConstants::GeV);

    // ---- Build seed list (auto or manual) -----------------------------------
    // Each entry is (dd_x, dd_y, dd_z) in DD4hep convention.
    std::vector<std::array<double, 3>> seedPositions3D;
    std::vector<std::array<double, 3>> seedDirections3D;

    if (m_autoSeed.value()) {
      // Hough Transform seeding using all measurements
      auto autoSeeds = findSeeds(measurements, gctx);

      for (const auto& sc : autoSeeds) {
        seedPositions3D.push_back({sc.x, sc.y, sc.z_start});
        // Direction from the seed's straight-line regression (dz = 1 by
        // construction); the caller normalises before creating the parameters.
        seedDirections3D.push_back({sc.dx, sc.dy, 1.0});
      }

      if (evtNum < 3) {
        debug() << "[ACTSProtoTracker] evt=" << evtNum
                << " Hough seeding: " << autoSeeds.size()
                << " seed(s) found" << endmsg;
        for (std::size_t si = 0; si < autoSeeds.size(); ++si) {
          debug() << "[ACTSProtoTracker]   seed " << si
                  << ": x=" << autoSeeds[si].x
                  << " y=" << autoSeeds[si].y
                  << " z_start=" << autoSeeds[si].z_start
                  << " votes=" << autoSeeds[si].nVotes
                  << " multiplicity=" << autoSeeds[si].multiplicity << endmsg;
        }
      }

      if (autoSeeds.empty()) {
        // Fallback: single seed at first surface
        double zFirst = allSurfaces.empty() ? -370.0
                      : allSurfaces.front()->center(gctx).x();
        seedPositions3D.push_back({0.0, 0.0, zFirst});
        seedDirections3D.push_back({0.0, 0.0, 1.0});
        if (evtNum < 3) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " Hough seeding found no seeds, using fallback." << endmsg;
        }
      }
    } else {
      // Manual seeding from Gaudi properties
      const auto& seedPosCfg = m_seedPositions.value();
      const auto& seedDirCfg = m_seedDirections.value();
      const std::size_t nManual = seedPosCfg.size() / 3;
      for (std::size_t i = 0; i < nManual; ++i) {
        seedPositions3D.push_back({seedPosCfg[3*i+0],
                                   seedPosCfg[3*i+1],
                                   seedPosCfg[3*i+2]});
        seedDirections3D.push_back({seedDirCfg[3*i+0],
                                    seedDirCfg[3*i+1],
                                    seedDirCfg[3*i+2]});
      }
    }

    const std::size_t nSeeds = seedPositions3D.size();

    std::size_t nTracks = 0;

    // Provisional per-seed winners, deduplicated at the end of the event.
    // Fingerprints store sets of hit indices (not surface geoIDs) for accurate
    // duplicate detection across seeds with different hit selection.
    struct AcceptedTrack {
      Acts::TrackIndexType   idx;    // into ckfTracks (states read at write time)
      std::set<std::size_t>  fp;     // hit-index fingerprint
      double                 chi2;
      int                    ndf;
      double                 ddx, ddy;  // seed transverse position (AtIP state)
      std::size_t            iSeed;
    };
    std::vector<AcceptedTrack> provisional;

    // Standard goodness-of-fit convention:
    //   chi² = Σ per-state innovation chi² (filled by ACTS during filtering)
    //   ndf  = Σ calibratedSize − n_fit_params   (5 helix params: loc0, loc1,
    //          phi, theta, q/p; ACTS's nDoF() returns Σ calibratedSize only).
    // For a well-fit track, chi²/ndf is centered at ≈ 1.
    constexpr int kHelixParams = 5;

    // ---- Fit options (seed-independent, shared by all seeds) ----------------
    Acts::PropagatorPlainOptions pOptions(gctx, m_mctx);
    pOptions.direction = Acts::Direction::Forward();
    pOptions.stepping.maxStepSize = m_maxStepSize.value();
    pOptions.maxSteps = static_cast<std::size_t>(m_maxPropSteps.value());

    Acts::CombinatorialKalmanFilterOptions<SNDTrackContainer> ckfOptions(
        gctx, m_mctx, std::cref(m_cctx),
        ckfExtensions,
        pOptions,
        true,
        true);

    Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory> kfOptions(
        gctx, m_mctx, std::cref(m_cctx), kfExtensions, pOptions);

    // ---- Per-seed fit helpers -----------------------------------------------
    auto chi2NdfOf = [&](const auto& track) -> double {
      const int ndf =
          std::max(1, static_cast<int>(track.nDoF()) - kHelixParams);
      return track.chi2() / static_cast<double>(ndf);
    };

    // NOTE: TrackStateCreator sets MeasurementFlag on outlier states too (plus
    // OutlierFlag), so outliers must be excluded explicitly — they are not
    // part of the fit and must not be fed to the frozen-hit KF refit.
    auto fingerprintOf = [](const auto& track) -> std::set<std::size_t> {
      std::set<std::size_t> fp;
      for (const auto& ts : track.trackStatesReversed()) {
        if (ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag) &&
            !ts.typeFlags().test(Acts::TrackStateFlag::OutlierFlag) &&
            ts.hasUncalibratedSourceLink()) {
          fp.insert(ts.getUncalibratedSourceLink()
                        .template get<SNDSourceLink>().index);
        }
      }
      return fp;
    };

    // Run one CKF pass from the given seed parameters; return the index (into
    // ckfTracks) of the best candidate: nMeas >= 3, lowest chi2/ndf. The
    // MaxChi2PerNdf cut and duplicate rejection are applied to the FINAL
    // track only, after the optional KF refit.
    auto runCKFPass = [&](const Acts::BoundTrackParameters& sp,
                          std::size_t iSeed)
        -> std::optional<Acts::TrackIndexType> {
      auto ckfResult = ckf.findTracks(sp, ckfOptions, ckfTracks);
      if (!ckfResult.ok()) {
        warning() << "[ACTSProtoTracker] evt=" << evtNum
                  << " seed=" << iSeed
                  << " CKF failed: " << ckfResult.error() << endmsg;
        return std::nullopt;
      }
      std::optional<Acts::TrackIndexType> bestIdx;
      double bestChi2Ndf = std::numeric_limits<double>::max();
      for (const auto& t : *ckfResult) {
        const std::size_t nMeas = t.nMeasurements();
        const double      c2n   = chi2NdfOf(t);
        debug() << "[ACTSProtoTracker] evt=" << evtNum
                << " seed=" << iSeed
                << " CKF nMeas=" << nMeas
                << " nHoles=" << t.nHoles()
                << " chi2=" << t.chi2()
                << " chi2/ndf=" << c2n << endmsg;
        if (nMeas < 3) continue;
        if (c2n < bestChi2Ndf) {
          bestChi2Ndf = c2n;
          bestIdx     = t.index();
        }
      }
      return bestIdx;
    };

    // Re-seed parameters from a fitted track: smoothed state at the most
    // upstream measurement, extrapolated back to the first surface with a
    // straight line (the setup is field-free), at the given loose covariance
    // so the refit stays measurement-dominated.
    auto reseedFrom = [&](const auto& track, const Acts::BoundSquareMatrix& cov)
        -> std::optional<Acts::BoundTrackParameters> {
      // trackStatesReversed walks tip -> stem, so the LAST measurement state
      // visited is the most upstream one.
      std::optional<Acts::BoundVector> bound;
      const Acts::Surface* surf = nullptr;
      for (const auto& ts : track.trackStatesReversed()) {
        if (!ts.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          continue;
        }
        if (ts.hasSmoothed()) {
          bound = ts.smoothed();
          surf  = &ts.referenceSurface();
        } else if (ts.hasFiltered()) {
          bound = ts.filtered();
          surf  = &ts.referenceSurface();
        }
      }
      if (!bound || surf == nullptr) return std::nullopt;
      const double phi   = (*bound)[Acts::eBoundPhi];
      const double theta = (*bound)[Acts::eBoundTheta];
      Acts::Vector3 dir(std::sin(theta) * std::cos(phi),
                        std::sin(theta) * std::sin(phi),
                        std::cos(theta));
      if (std::abs(dir.x()) < 1e-6) return std::nullopt;
      Acts::Vector3 pos = surf->localToGlobal(
          gctx,
          Acts::Vector2((*bound)[Acts::eBoundLoc0], (*bound)[Acts::eBoundLoc1]),
          dir);
      const Acts::Surface* front = allSurfaces.front();
      const double xFront = front->center(gctx).x();
      pos += dir * ((xFront - pos.x()) / dir.x());
      Acts::Vector4 pos4(pos.x(), pos.y(), pos.z(), 0.0);
      auto res = Acts::BoundTrackParameters::create(
          gctx, front->getSharedPtr(), pos4, dir,
          (*bound)[Acts::eBoundQOverP], cov, Acts::ParticleHypothesis::muon());
      if (!res.ok()) return std::nullopt;
      return *res;
    };

    // =========================================================================
    // STEPS 5-7: Loop over seeds
    // =========================================================================
    for (std::size_t iSeed = 0; iSeed < nSeeds; ++iSeed) {


      // ---- Seed position (DD4hep → ACTS coordinate swap) -------------------
      // DD4hep convention: x=transverse X, y=transverse Y, z=beam.
      const double dd_x = seedPositions3D[iSeed][0];  // transverse X [mm]
      const double dd_y = seedPositions3D[iSeed][1];  // transverse Y [mm]
      // dd_z (beam) is used to find the starting surface
      const double dd_z = seedPositions3D[iSeed][2];  // beam Z [mm]

      // Always use the first surface as the KF seed reference surface.
      // The Hough z_start is used only to inform hit selection, not the KF seed.
      // Starting from the first surface is safe — the KF handles empty layers
      // (holes) gracefully with loose covariance.
      const Acts::Surface* sfSeed = allSurfaces.front();

      Acts::Vector4 seedPos4;
      // Beam coordinate goes into ePos0 (ACTS X, after X<->Z swap in geometry).
      //
      // Surface frame: the geometry applies rot90Y = R_Y(pi/2), whose columns
      // are the images of the local axes:
      //     local x -> (0, 0, -1)     local y -> (0, 1, 0)     local z -> (1, 0, 0)
      // so local z (the normal) is the beam, as intended, but local x maps to
      // MINUS global Z. Inverting, a global offset (0, dy, dz) from the surface
      // centre has bound coordinates (loc0, loc1) = (-dz, dy).
      //
      // The calibrator defines the measurement as (loc0, loc1) = (dd_x, dd_y),
      // so consistency requires global Z = -dd_x and global Y = dd_y. Filling
      // ePos2 with +dd_x instead put the seed at -dd_x in loc0, i.e. 2*|dd_x|
      // away from its own track: harmless at dd_x = 0, but 88 mm for a muon at
      // dd_x = -44 mm, whose first-layer chi2 (88/10.1)^2 = 76 then exceeds
      // Chi2CutOff = 70 and the CKF loses the track entirely.
      seedPos4[Acts::ePos0] = sfSeed->center(gctx).x();  // beam coord (ACTS X)
      seedPos4[Acts::ePos1] =  dd_y;  // DD4hep Y → global Y  → eBoundLoc1
      seedPos4[Acts::ePos2] = -dd_x;  // DD4hep X → -global Z → eBoundLoc0
      seedPos4[Acts::eTime] = 0.0;

      // ---- Seed direction (DD4hep → ACTS coordinate swap) ------------------
      // Same mapping as the position above: beam → ACTS X, dd_y → ACTS Y,
      // dd_x → MINUS ACTS Z. The old (dd_dz, dd_dx, dd_dy) swap was wrong in
      // both transverse components; it stayed invisible only because the seed
      // direction used to be hardcoded to (0, 0, 1), which maps to (1, 0, 0)
      // either way. It matters now that the Hough seeds carry real slopes.
      const double dd_dx = seedDirections3D[iSeed][0];
      const double dd_dy = seedDirections3D[iSeed][1];
      const double dd_dz = seedDirections3D[iSeed][2];

      Acts::Vector3 seedDir(dd_dz, dd_dy, -dd_dx);
      if (seedDir.norm() < 1e-6) seedDir = Acts::Vector3(1.0, 0.001, 0.001);
      seedDir = seedDir.normalized();

      // ---- Create seed parameters -------------------------------------------
      // Hit selection is delegated to the CKF's MeasurementSelector (chi2 cut).
      auto seedParamsResult = Acts::BoundTrackParameters::create(
          gctx, sfSeed->getSharedPtr(), seedPos4, seedDir, seedQoverP,
          seedCov, Acts::ParticleHypothesis::muon());

      if (!seedParamsResult.ok()) {
        warning() << "[ACTSProtoTracker] evt=" << evtNum
                  << " seed=" << iSeed
                  << " seed creation failed: "
                  << seedParamsResult.error() << endmsg;
        continue;
      }
      const auto& seedParams = *seedParamsResult;

      // ---- Seed surface diagnostic ------------------------------------------
      if (evtNum < 1 && iSeed == 0) {
        const Acts::Surface* sf = &seedParams.referenceSurface();
        debug() << "[ACTSProtoTracker] DIAG evt=0 seed=0: seedSurf geoId="
                << sf->geometryId()
                << " center=" << sf->center(gctx).transpose() << endmsg;
        // allSurfaces is the same list injected into the navigator
        bool found = std::any_of(
            allSurfaces.begin(), allSurfaces.end(),
            [sf](const Acts::Surface* s){ return s == sf; });
        debug() << "[ACTSProtoTracker] DIAG: seedSurf in allSurfaces=" << found
                << " allSurfaces.size()=" << allSurfaces.size() << endmsg;
      }

      // ---- Run CKF: best chi2/ndf candidate for this seed --------------------
      auto bestIdxOpt = runCKFPass(seedParams, iSeed);
      if (!bestIdxOpt) continue;
      std::set<std::size_t> bestFp =
          fingerprintOf(ckfTracks.getTrack(*bestIdxOpt));

      // ---- Final KalmanFitter refit on the frozen hit set --------------------
      // The hit set selected by the CKF is refit as-is (no re-selection
      // possible), giving unbiased parameters and covariances.
      if (m_finalRefit.value() && bestFp.size() >= 3) {
        std::vector<Acts::SourceLink> fitSLinks;
        fitSLinks.reserve(bestFp.size());
        for (std::size_t mIdx : bestFp) {
          SNDSourceLink ssl;
          ssl.index = mIdx;
          ssl.setGeometryId(measurements[mIdx].surface->geometryId());
          fitSLinks.emplace_back(ssl);
        }
        // Start from the CKF solution; the loose covariance keeps the prior
        // uninformative, so the refit stays unbiased.
        Acts::BoundTrackParameters kfStart = seedParams;
        if (auto rs = reseedFrom(ckfTracks.getTrack(*bestIdxOpt), seedCov)) {
          kfStart = *rs;
        }
        auto kfResult = kf.fit(fitSLinks.begin(), fitSLinks.end(), kfStart,
                               kfOptions, allSurfaces, ckfTracks);
        if (evtNum < 3 && kfResult.ok()) {
          debug() << "[ACTSProtoTracker][DIAG-KF] evt=" << evtNum
                  << " seed=" << iSeed
                  << " fitSLinks=" << fitSLinks.size()
                  << " nMeas=" << (*kfResult).nMeasurements()
                  << " nHoles=" << (*kfResult).nHoles()
                  << " chi2=" << (*kfResult).chi2() << endmsg;
        }
        if (kfResult.ok() && (*kfResult).nMeasurements() >= 3) {
          bestIdxOpt = (*kfResult).index();
          bestFp     = fingerprintOf(*kfResult);
        } else if (kfResult.ok()) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " seed=" << iSeed
                    << " KF refit kept too few measurements, keeping CKF track"
                    << endmsg;
        } else {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " seed=" << iSeed
                    << " KF refit failed: " << kfResult.error()
                    << ", keeping CKF track" << endmsg;
        }
      }

      // ---- Track-level acceptance (final track only) -------------------------
      auto finalCand = ckfTracks.getTrack(*bestIdxOpt);
      const int    ndf     = std::max(
          1, static_cast<int>(finalCand.nDoF()) - kHelixParams);
      const double chi2    = finalCand.chi2();
      const double chi2Ndf = chi2 / static_cast<double>(ndf);

      debug() << "[ACTSProtoTracker] evt=" << evtNum
              << " seed=" << iSeed
              << " final nMeas=" << finalCand.nMeasurements()
              << " nHoles=" << finalCand.nHoles()
              << " chi2=" << chi2
              << " ndf=" << ndf
              << " chi2/ndf=" << chi2Ndf << endmsg;

      if (chi2Ndf > m_maxChi2PerNdf.value()) {
        debug() << "[ACTSProtoTracker] evt=" << evtNum
                << " seed=" << iSeed
                << " rejected: chi2/ndf=" << chi2Ndf << endmsg;
        continue;
      }

      // ---- Seed cleaning ----------------------------------------------------
      // Remove the accepted track's hits from the source-link pool so the
      // remaining seeds can only build tracks from unused measurements.
      // sortedSLinks stays sorted (erase preserves order) and the accessor
      // holds a pointer to the vector, so no re-wiring is needed.
      if (m_seedCleaning.value()) {
        sortedSLinks.erase(
            std::remove_if(sortedSLinks.begin(), sortedSLinks.end(),
                           [&bestFp](const Acts::SourceLink& sl) {
                             return bestFp.count(
                                        sl.get<SNDSourceLink>().index) > 0;
                           }),
            sortedSLinks.end());
      }

      // Defer the accept/write decision: collect this seed's winner and
      // deduplicate the whole event at the end, best-first (see below).
      provisional.push_back({*bestIdxOpt, bestFp, chi2, ndf, dd_x, dd_y, iSeed});

    }  // end loop over seeds

    // ---- Write helper -------------------------------------------------------
    auto writeTrack = [&](Acts::TrackIndexType bestIdx, double chi2, int ndf,
                          double dd_x, double dd_y, std::size_t iSeed) {
      auto finalTrack = ckfTracks.getTrack(bestIdx);
        // ---- Write output --------------------------------------------------
        auto track = output->create();
        track.setType(1);
        track.setChi2(static_cast<float>(chi2));
        track.setNdf(ndf);  // already Σ-dim − 5; downstream uses chi2/ndf directly

        // Store seed transverse position as the first TrackState (location=AtIP).
        {
          edm4hep::TrackState seedState{};
          seedState.location = edm4hep::TrackState::AtIP;
          seedState.D0       = static_cast<float>(dd_x);
          seedState.Z0       = static_cast<float>(dd_y);
          track.addToTrackStates(seedState);
        }

        // Write per-surface filtered states for event display / analysis.
        try {
          auto tipIdx = finalTrack.tipIndex();
          auto& mutableTraj = ckfTracks.trackStateContainer();
          std::vector<edm4hep::TrackState> collected;
          while (true) {
            auto ts = mutableTraj.getTrackState(tipIdx);
            if (ts.hasCalibrated()) {
              edm4hep::TrackState edm4ts;
              edm4ts.location  = edm4hep::TrackState::AtOther;
              edm4ts.D0        = 0.0f;
              // Z0 is unused on AtOther states (only AtIP uses it for seed Y).
              // Carry the per-state innovation chi² here so downstream tools can
              // see WHICH surface drives the fit quality — track-total chi² is
              // useless for that.
              edm4ts.Z0        = static_cast<float>(ts.chi2());
              const auto& surf  = ts.referenceSurface();
              const float beamZ = static_cast<float>(surf.center(gctx).x());

              // stereoTiltZ ≈ 0 for SiPad (no stereo rotation).
              const Acts::Vector3 localYinGlobal =
                  surf.transform(gctx).rotation() * Acts::Vector3::UnitY();
              const double stereoTiltZ = localYinGlobal.z();  // sin(+α) or sin(−α)
              const double cosAlpha    = localYinGlobal.y();  // cos(α), > 0 always
              // omega≈0 for SiPad (no stereo tilt).
              edm4ts.omega = static_cast<float>(stereoTiltZ);

              // Measurement convention on stereo surfaces:
              //   loc0 = cos(α)·dd_x − sin(α)·dd_y   (= −eBoundLoc0_geometric)
              //   loc1 = cos(α)·dd_y + sin(α)·dd_x   (=  eBoundLoc1_geometric)
              // Inverse rotation recovers DD4hep coordinates:
              //   dd_x =  loc0·cosAlpha + loc1·stereoTiltZ
              //   dd_y = −loc0·stereoTiltZ + loc1·cosAlpha
              // For non-stereo surfaces: stereoTiltZ=0, cosAlpha=1 → dd_x=loc0, dd_y=loc1.
              auto toGlobalX = [&](double loc0, double loc1) -> float {
                return static_cast<float>(loc0 * cosAlpha + loc1 * stereoTiltZ);
              };
              auto toGlobalY = [&](double loc0, double loc1) -> float {
                return static_cast<float>(-loc0 * stereoTiltZ + loc1 * cosAlpha);
              };

              if (ts.hasSmoothed()) {
                edm4ts.phi       = static_cast<float>(
                    ts.smoothed()[Acts::eBoundPhi]);
                edm4ts.tanLambda = static_cast<float>(
                    std::tan(M_PI / 2.0 - ts.smoothed()[Acts::eBoundTheta]));
                const double loc0 = ts.smoothed()[Acts::eBoundLoc0];
                const double loc1 = ts.smoothed()[Acts::eBoundLoc1];
                edm4ts.D0 = static_cast<float>(loc0);  // raw strip measurement for stereo pairing
                edm4ts.referencePoint = edm4hep::Vector3f{
                    toGlobalX(loc0, loc1),
                    toGlobalY(loc0, loc1),
                    beamZ};
              } else if (ts.hasFiltered()) {
                edm4ts.phi       = static_cast<float>(
                    ts.filtered()[Acts::eBoundPhi]);
                edm4ts.tanLambda = 0.0f;
                const double loc0 = ts.filtered()[Acts::eBoundLoc0];
                const double loc1 = ts.filtered()[Acts::eBoundLoc1];
                edm4ts.D0 = static_cast<float>(loc0);  // raw strip measurement for stereo pairing
                edm4ts.referencePoint = edm4hep::Vector3f{
                    toGlobalX(loc0, loc1),
                    toGlobalY(loc0, loc1),
                    beamZ};
              } else {
                edm4ts.phi       = 0.0f;
                edm4ts.tanLambda = 0.0f;
                edm4ts.referencePoint = edm4hep::Vector3f{0.f, 0.f, beamZ};
              }

              collected.push_back(edm4ts);
            }

            // TEMP DIAG (Issue 2/3) — fires for EVERY state (Meas+Outlier+Hole)
            // for first 3 events. For SiPad investigation: a SiPad surface
            // should appear here with isHole=1 (or isOutl=1 if chi2 fails);
            // pred(loc0,loc1) vs meas(loc0,loc1) shows whether the chi2 is
            // huge (rejected) or small (would be accepted).
            if (evtNum < 3) {
              const auto& surfD = ts.referenceSurface();
              const float zD = static_cast<float>(surfD.center(gctx).x());
              const auto tf = ts.typeFlags();
              const bool isMeas = tf.test(Acts::TrackStateFlag::MeasurementFlag);
              const bool isOutl = tf.test(Acts::TrackStateFlag::OutlierFlag);
              const bool isHole = tf.test(Acts::TrackStateFlag::HoleFlag);
              const char* region = "SiP";
              double pL0 = 0., pL1 = 0., sL0 = -1., sL1 = -1.;
              if (ts.hasPredicted()) {
                pL0 = ts.predicted()[Acts::eBoundLoc0];
                pL1 = ts.predicted()[Acts::eBoundLoc1];
                const auto& pc = ts.predictedCovariance();
                sL0 = std::sqrt(pc(Acts::eBoundLoc0, Acts::eBoundLoc0));
                sL1 = std::sqrt(pc(Acts::eBoundLoc1, Acts::eBoundLoc1));
              }
              std::ostringstream ms;
              if (ts.hasCalibrated()) {
                const unsigned cs = ts.calibratedSize();
                if (cs == 2) {
                  auto c2 = ts.template calibrated<2>();
                  ms << " meas(" << c2[0] << "," << c2[1] << ") cs=2";
                } else if (cs == 1) {
                  auto c1 = ts.template calibrated<1>();
                  ms << " meas(" << c1[0] << ",-) cs=1";
                }
              }
              info() << "[ACTSProtoTracker][DIAG-TS] evt=" << evtNum
                     << " seed=" << iSeed
                     << " z=" << zD << " " << region
                     << " M=" << isMeas << " O=" << isOutl << " H=" << isHole
                     << " pred(" << pL0 << "," << pL1 << ")"
                     << " sig(" << sL0 << "," << sL1 << ")"
                     << ms.str()
                     << endmsg;
            }

            if (!ts.hasPrevious()) break;
            tipIdx = ts.previous();
          }  // end while (per-surface collect)

          // ---- Pair-average U/V stereo partners ---------------------------
          // Sort collected states by global z (beam axis), then walk pairs.
          // Pairing condition: both have nonzero omega (stereo), opposite
          // sign, and |Δz| < 5 mm.  Output omega is set to 0 to signal
          // "already paired" to downstream consumers (event display, diag).
          std::sort(collected.begin(), collected.end(),
                    [](const edm4hep::TrackState& a,
                       const edm4hep::TrackState& b) {
                      return a.referencePoint.z < b.referencePoint.z;
                    });
          for (std::size_t i = 0; i < collected.size(); ) {
            bool didPair = false;
            if (i + 1 < collected.size()) {
              const auto& a = collected[i];
              const auto& b = collected[i + 1];
              const float dz = std::abs(a.referencePoint.z - b.referencePoint.z);
              if (std::abs(a.omega) > 0.01f &&
                  std::abs(b.omega) > 0.01f &&
                  a.omega * b.omega < 0.0f &&
                  dz < 5.0f) {
                edm4hep::TrackState avg;
                avg.location  = edm4hep::TrackState::AtOther;
                avg.D0        = 0.0f;     // loc0 not meaningful after pairing
                // Z0 carries per-state chi²; merged pair gets the sum (each
                // partner contributed one innovation update with E[chi²]=1).
                avg.Z0        = a.Z0 + b.Z0;
                avg.omega     = 0.0f;     // signals "paired/non-stereo"
                avg.phi       = 0.5f * (a.phi + b.phi);
                avg.tanLambda = 0.5f * (a.tanLambda + b.tanLambda);
                avg.referencePoint = edm4hep::Vector3f{
                    0.5f * (a.referencePoint.x + b.referencePoint.x),
                    0.5f * (a.referencePoint.y + b.referencePoint.y),
                    0.5f * (a.referencePoint.z + b.referencePoint.z)};
                track.addToTrackStates(avg);
                i += 2;
                didPair = true;
              }
            }
            if (!didPair) {
              track.addToTrackStates(collected[i]);
              ++i;
            }
          }
        } catch (const std::exception& e) {
          warning() << "[ACTSProtoTracker] evt=" << evtNum
                    << " seed=" << iSeed
                    << " trackStates iteration failed: " << e.what() << endmsg;
        }
        ++nTracks;
    };  // end writeTrack

    // ---- Deduplicate the event's tracks, best-first ------------------------
    // Sort accepted candidates by chi2/ndf ascending, then greedily keep a
    // candidate only if it does not overlap an already-kept (better) track by
    // more than DuplicateOverlapFraction of the smaller hit set. The survivor
    // of each overlap group is therefore always the best-quality fit,
    // independent of the order the seeds were processed in.
    std::sort(provisional.begin(), provisional.end(),
              [](const AcceptedTrack& a, const AcceptedTrack& b) {
                return a.chi2 / std::max(1, a.ndf) <
                       b.chi2 / std::max(1, b.ndf);
              });
    const double dupFrac = m_dupOverlapFraction.value();
    std::vector<const AcceptedTrack*> survivors;
    for (const auto& cand : provisional) {
      bool dup = false;
      for (const auto* kept : survivors) {
        std::size_t nShared = 0;
        for (std::size_t h : cand.fp) if (kept->fp.count(h)) ++nShared;
        const double smaller =
            static_cast<double>(std::min(cand.fp.size(), kept->fp.size()));
        if (smaller > 0 && nShared / smaller > dupFrac) { dup = true; break; }
      }
      if (dup) {
        debug() << "[ACTSProtoTracker] evt=" << evtNum << " seed=" << cand.iSeed
                << " dropped as duplicate (kept a better-chi2/ndf track)" << endmsg;
        continue;
      }
      survivors.push_back(&cand);
    }
    for (const auto* sv : survivors) {
      writeTrack(sv->idx, sv->chi2, sv->ndf, sv->ddx, sv->ddy, sv->iSeed);
    }

    info() << "[ACTSProtoTracker] evt=" << evtNum
           << " measurements=" << measurements.size()
           << " CKF tracks=" << nTracks << endmsg;

    return StatusCode::SUCCESS;
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in execute(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in execute()." << endmsg;
    return StatusCode::FAILURE;
  }
}

// ---------------------------------------------------------------------------
// finalize()
// ---------------------------------------------------------------------------

StatusCode ACTSProtoTracker::finalize() {
  try {
    m_siPadHandle.reset();
    m_outputHandle.reset();
    info() << "[ACTSProtoTracker] Done. Total events processed: "
           << m_eventCount.load() << endmsg;
    return Gaudi::Algorithm::finalize();
  } catch (const std::exception& e) {
    error() << "[ACTSProtoTracker] Exception in finalize(): "
            << e.what() << endmsg;
    return StatusCode::FAILURE;
  } catch (...) {
    error() << "[ACTSProtoTracker] Unknown exception in finalize()." << endmsg;
    return StatusCode::FAILURE;
  }
}

DECLARE_COMPONENT(ACTSProtoTracker)