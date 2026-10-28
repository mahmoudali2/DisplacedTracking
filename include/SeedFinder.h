// SeedFinder — N-plet seeding + helix extension, for events too dense for an exhaustive
// C(N,k) search over all layers.
//
// 1. Seed layer sets are derived once from LayerGeometry: straight rays from the origin are
//    swept in theta and every run of NumSeedLayers consecutive crossed layers (optionally
//    allowing skipped layers) is a candidate seed layer set. Policy picks all runs
//    ("consecutive", displacement-friendly), only the first ("innermost") or only the last
//    ("outermost"), or user-given sets ("explicit").
// 2. Seeds: one hit per seed layer, grown hit by hit with cheap pairwise cuts (azimuth
//    difference, 3D distance), then a circle fit (via the component's own fitCircleNHits,
//    passed in as a callback) and a z-vs-arc-length consistency check. A 3-hit circle fit
//    has chi2 == 0 by construction, so the z residual is the real triplet discriminant.
// 3. Extension: the seed helix is propagated layer by layer (HelixExtrapolator) in the
//    configured direction; on each crossed layer the nearest compatible hits (measured on
//    the hit's own sensor plane) are tried, refitting each time; a hit is kept only if the
//    refit chi2/ndf stays below the cut. The track stops after MaxMissedLayers consecutive
//    layers without an accepted hit.
// 4. Candidates are de-duplicated and sorted (more hits first, then lower chi2/ndf). The
//    caller takes them greedily, skipping any that share hits with an accepted track.
#ifndef TRACKING_SEEDFINDER_H
#define TRACKING_SEEDFINDER_H

#include "HelixExtrapolator.h"
#include "LayerGeometry.h"
#include "TrackQuality.h"

#include "edm4hep/TrackerHitPlane.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include <Eigen/Dense>

#include <cstdint>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <string>
#include <utility>
#include <vector>

namespace k4rt {

struct SeedConfig {
  int    numSeedLayers   = 3;              ///< hits per seed (>= 3), when `passes` is empty
  /// Iterative seeding: hits per seed for each pass, e.g. {5, 4, 3}. Later passes seed only on
  /// hits not claimed by an earlier pass's candidates. Empty = one pass with numSeedLayers.
  std::vector<int> passes;
  int    claimMinHits    = 4;     ///< a candidate with >= this many hits claims its hits
  std::string policy     = "consecutive";  ///< consecutive | innermost | outermost | explicit
  std::vector<std::vector<int>> explicitSets;  ///< compositeIDs, used when policy == explicit
  int    layerSkip       = 0;     ///< layers that may be skipped inside one seed window
  double maxDeltaPhi     = 0.3;   ///< rad, azimuth difference between consecutive seed hits
  double maxHitDist_cm   = 100.;  ///< 3D distance between consecutive seed hits
  double minRadius_cm    = 0.;    ///< seed circle radius floor (a pT floor); 0 = off
  double maxZResid_cm    = 0.5;   ///< max |z - line| of seed hits about the z-vs-s fit
  double seedMaxChi2NDF  = 10.;   ///< only bites for numSeedLayers > 3
  std::string direction  = "Outward";  ///< Outward | Inward | Bidirectional
  double window_cm       = 1.0;   ///< hit search radius at the crossing ...
  double windowSlope     = 0.02;  ///< ... growing by this many cm per cm of extrapolation
  double extMaxChi2NDF   = 10.;   ///< refit chi2/ndf a new hit must keep
  int    maxMissedLayers = 1;
  int    maxCandPerLayer = 3;     ///< hits tried per crossed layer, nearest first
  int    minTrackHits    = 3;
  long   maxSeeds        = 200000;  ///< per event safety cap on seeds passing the cheap pre-check

  // ── BestFirst mode ────────────────────────────────────────────────────────────
  // Seed windows inside-out; within a window the seeds best in kink chi2 are extended first;
  // an accepted track claims its hits at once, so later seeds that use them are never
  // extended. Hits are attached by a multiple-scattering-aware kink chi2 (TrackQuality.h),
  // each track gets a confirmation pass, and sequential pieces of one trajectory are merged.
  std::string mode       = "Classic";  ///< Classic | BestFirst
  double kinkMaxChi2Triplet  = 16.0;   ///< a hit is attached only if the refit's kink chi2 grows by less than this (2 dof)
  double seedKinkMaxChi2NDF  = 9.0;    ///< seed acceptance (kink chi2/ndf)
  double trackKinkMaxChi2NDF = 5.0;    ///< a finished track is kept (and claims its hits) below this
  int    confirmIterations   = 2;      ///< confirmation passes after extension
  bool   mergePieces         = true;   ///< merge sequential pieces of one trajectory
  /// Hit-pattern confirmation of short tracks: a track with <= shortTrackMaxHits hits is kept
  /// only if it has <= shortTrackMaxHoles holes, i.e. layers its helix crosses (between its
  /// first and last hit, and outward until it leaves the tracker) with no hit near the crossing.
  /// Inward holes are not counted: a displaced track legitimately has none inside its origin.
  int    shortTrackMaxHits   = 4;
  int    shortTrackMaxHoles  = 1;
  /// Short tracks (<= shortTrackMaxHits hits) must also have fitted pT >= this (GeV); 0 = off.
  /// A 3-hit track has a single constraint, and at low pT its expected scattering is large,
  /// so random hit triplets pass any kink cut; in dense events their momentum is what
  /// separates them from real tracks.
  double shortTrackMinPt     = 0.0;
  /// After a track is accepted, unclaimed hits lying on its trajectory within this distance
  /// (cm, on their own sensor plane) are claimed too, without being added: overlapping
  /// sensors give a particle a second hit on a layer, which would otherwise seed a parallel
  /// duplicate. 0 = off.
  double shadowHitWindow_cm  = 0.1;
  double bField_T            = 2.0;
  double mass_GeV            = 0.10566;  ///< mass hypothesis for the scattering angle
};

/// A found track candidate. Hits are ordered along the direction of travel.
struct SECandidate {
  std::vector<size_t> gIdx;
  double chi2ndf = 0, x0 = 0, y0 = 0, R = 0;
  Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
  std::vector<double> residuals;
  int nSeedHits = 0;
  // BestFirst only
  double kinkChi2 = -1.0;  ///< multiple-scattering-aware kink chi2 (TrackQuality.h)
  int    kinkNdf  = 0;
  bool   looper   = false; ///< turns by more than half a circle, or overlaps another arc of the same circle
  bool   merged   = false; ///< absorbed into another candidate (internal)
};

struct SeedStats {
  long seedLayerSets = 0;
  long rawTriplets = 0;  ///< hit combinations enumerated, before any check
  long seeds = 0;        ///< combinations passing the cheap geometric pre-check (this is what MaxSeedsPerEvent caps)
  long seedsAfterZ = 0, extended = 0, candidates = 0;
  // BestFirst only
  long skippedClaimed = 0, rejectedTracks = 0, merged = 0, loopers = 0;
  long confirmAdded = 0, confirmSwapped = 0, confirmRemoved = 0, rejectedHoles = 0, rejectedShortPt = 0;
  long shadowClaimed = 0;
  bool hitSeedCap = false;
};

class SeedFinder {
public:
  /// Circle fit: (ordered hits) -> x0, y0, R (cm), chi2/ndf, 3x3 cov, residuals. This is the
  /// component's fitCircleNHits, so there is one fit, not two.
  using FitFn = std::function<bool(const std::vector<edm4hep::TrackerHitPlane>&,
                                   double&, double&, double&, double&,
                                   Eigen::Matrix3d&, std::vector<double>&)>;
  /// Sensor plane of a hit: origin (cm) and unit normal. false if unknown.
  using PlaneFn = std::function<bool(uint64_t cellID, Eigen::Vector3d&, Eigen::Vector3d&)>;

  SeedFinder(const LayerGeometry& geo, SeedConfig cfg, FitFn fit, PlaneFn plane);

  /// Derive the seed layer sets from the geometry. Call once (initialize()).
  /// Returns problems found (e.g. unknown compositeIDs in explicit sets).
  std::vector<std::string> prepare();
  const std::vector<std::vector<int>>& seedLayerSets() const { return m_sets; }  ///< first pass
  const std::vector<int>& passes() const { return m_passes; }
  size_t numSets(int n) const { auto it = m_setsByN.find(n); return it == m_setsByN.end() ? 0 : it->second.size(); }
  const SeedConfig& config() const { return m_cfg; }

  /// Find candidates among `pool` = (global hit index, compositeID) of usable hits.
  std::vector<SECandidate> run(const edm4hep::TrackerHitPlaneCollection& hits,
                               const std::vector<std::pair<size_t, int>>& pool,
                               SeedStats& stats) const;

private:
  struct LayerHits {          // hits of one layer, sorted by azimuth for window lookups
    std::vector<size_t> gIdx;
    std::vector<double> phi;
  };

  std::vector<SECandidate> runBestFirst(const edm4hep::TrackerHitPlaneCollection& hits,
                                        const std::vector<std::pair<size_t, int>>& pool,
                                        SeedStats& stats) const;

  /// Seed layer sets of size n from the ray sweep (policy consecutive/innermost/outermost).
  std::vector<std::vector<int>> buildSets(int n, std::vector<std::string>& problems) const;

  bool extend(const edm4hep::TrackerHitPlaneCollection& hits,
              const std::map<int, LayerHits>& byLayer,
              const std::unordered_map<size_t, int>& cidOf,
              std::vector<size_t>& track, int dir,
              std::vector<int>& visited, double& x0, double& y0, double& R, double& chi2ndf,
              Eigen::Matrix3d& cov, std::vector<double>& res) const;

  const LayerGeometry&          m_geo;
  SeedConfig                    m_cfg;
  FitFn                         m_fit;
  PlaneFn                       m_plane;
  std::vector<std::vector<int>> m_sets;                    ///< sets of the first pass
  std::map<int, std::vector<std::vector<int>>> m_setsByN;  ///< seed layer sets per seed size
  std::vector<int>              m_passes;                  ///< seed size of each pass
};

}  // namespace k4rt

#endif  // TRACKING_SEEDFINDER_H
