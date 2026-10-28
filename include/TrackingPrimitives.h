// TrackingPrimitives — the pieces of CombinatorialTrackFinder's pattern recognition that
// do not depend on the per-event search state, lifted out of findTracks() so that a seed
// finder can reuse them instead of carrying a second copy of every cut.
//
// Everything here is a pure function of its arguments. The bodies are moved verbatim from
// findTracks(); only the access to inputs changed (explicit parameters instead of [&]
// capture), so the arithmetic and its order are unchanged.
//
// Units follow the call sites they came from: hit positions in mm (EDM4hep), circle
// parameters and distance cuts in cm, B in tesla.
#ifndef TRACKING_TRACKINGPRIMITIVES_H
#define TRACKING_TRACKINGPRIMITIVES_H

#include "edm4hep/TrackerHitPlane.h"
#include "edm4hep/TrackerHitPlaneCollection.h"

#include <Eigen/Dense>

#include <cstddef>
#include <limits>
#include <unordered_map>
#include <vector>

class MsgStream;

namespace k4rt {

/// One unused hit offered to the search, with its decoded layer identity.
struct HitInfo {
  size_t index;  ///< index into the full input TrackerHitPlaneCollection
  edm4hep::TrackerHitPlane hit;
  int layerID;
  int typeID;
  int compositeID;

  HitInfo(size_t idx, const edm4hep::TrackerHitPlane& h, int layer, int type, int composite)
      : index(idx), hit(h), layerID(layer), typeID(type), compositeID(composite) {}
};

/// A fitted hit combination and the track parameters derived from it.
struct CombResult {
  std::vector<size_t> hiiVec;   // indices into allHitInfo
  std::vector<size_t> gIdxVec;  // indices into the full hit collection
  double chi2ndf   = std::numeric_limits<double>::max();
  double avgProxCm = std::numeric_limits<double>::max();  // avg edep-deviation score (keV)
  double x0 = 0, y0 = 0, radius = 0;
  Eigen::Matrix3d fitCov = Eigen::Matrix3d::Zero();
  std::vector<double> residuals;
  int nHits = 0;
  // Track parameters stored after N-hit fit (used to build seed state)
  double charge    = 0.0;
  double d0_cm     = 0.0;
  double phi_rad   = 0.0;
  double omega_mm  = 0.0;
  double z0_cm     = 0.0;
  double tanLambda = 0.0;
  double pT        = 0.0;
  // SeedAndExtend BestFirst only: multiple-scattering-aware kink chi2, looper flag
  double kinkChi2  = -1.0;
  int    kinkNdf   = 0;
  bool   looper    = false;
};

/// Cut values for the hit-set predicates. A value <= 0 disables that cut, as before.
struct RoadCutConfig {
  double maxConsecDeltaPhi     = -1.0;  ///< rad
  double maxComboPhiSpread     = -1.0;  ///< rad
  double maxConsecHitDist      = -1.0;  ///< cm
  double maxPairHitDist        = -1.0;  ///< cm
  double maxComboTanLambdaDev  = -1.0;
};

/// Region of a compositeID: 0 barrel, +1 +endcap, -1 -endcap.
///
/// The encoding is barrel 100*system+layer, endcap +-(1000+100*system+layer), so the old
/// test "cid >= 1000 means endcap" is WRONG as soon as a barrel system ID reaches 10
/// (IDEA o1v04: Si-wrapper barrel = system 23 -> cid 2300, outer vertex barrel = system 31
/// -> cid 3103 were treated as +endcap). The authoritative answer is the readout's region
/// field, collected per compositeID at initialize(); pass that map here. Without a map
/// (or for an unknown cid) the old range rule is used, which is exact when every system
/// ID is < 10.
int regionOfCid(int cid, const std::unordered_map<int, int>* cidRegion);

/// Everything the predicates read besides the candidate itself.
/// `dbg` may be null; when set, rejection reasons are written to it at DEBUG level.
struct RoadCutContext {
  const edm4hep::TrackerHitPlaneCollection& hits;
  const std::unordered_map<size_t, int>&    gIdxToComposite;
  const RoadCutConfig&                      cfg;
  MsgStream*                                dbg = nullptr;
  const std::unordered_map<int, int>*       cidRegion = nullptr;  ///< see regionOfCid()
};

/// phi road: consecutive-layer |dphi|, total phi spread, endcap |z|-r monotonicity.
/// All phi cuts are bypassed for near-beam combos (theta < 15 deg or > 160 deg).
bool passesRoadCuts(const std::vector<size_t>& gIdxVec, const RoadCutContext& ctx);

/// 3D distance between consecutive hits and between any two hits.
bool passesDistanceCuts(const std::vector<size_t>& gIdxVec, const RoadCutContext& ctx);

/// Per-pair dz/dr_xy must stay within maxComboTanLambdaDev of its median.
/// Mixed barrel+endcap combos always pass.
bool passesTanLambdaCut(const std::vector<size_t>& gIdxVec, const RoadCutContext& ctx);

// ── Circle -> helix pieces ────────────────────────────────────────────────────────────
// Shared by the N-hit track-state block and the AtLastHit block, which used to carry
// separate copies of the same arithmetic.

/// Mean normalised cross-product sin(bend) over consecutive triplets of ordered hits.
struct BendSense {
  double sinBendMean = 0.0;
  double avgChord    = 0.0;  ///< mm
  int    nTriplets   = 0;
};
BendSense bendSense(const std::vector<edm4hep::TrackerHitPlane>& orderedHits);

/// Charge sign from the bend and the field sign: bSign * (sinBend < 0 ? +1 : -1).
double chargeFromBend(double sinBendMean, double bFieldT);

/// Transverse momentum of a circle of radius r_cm in field B (GeV/c).
double ptFromRadius(double bFieldT, double r_cm);

/// d0 (cm), phi0 (rad, wrapped to [-pi, pi]) and omega (1/mm) from a circle and a charge.
struct Perigee2D {
  double d0_cm    = 0.0;
  double phi_rad  = 0.0;
  double omega_mm = 0.0;
};
Perigee2D perigeeFromCircle(double cx_cm, double cy_cm, double r_cm, double charge);

/// Linear fit z = z0 + tanLambda * s over the ordered hits, s = charge * r * dAngle
/// measured from the first hit. NOTE: s omits the 1/cos(lambda) factor, which is harmless
/// for this 2-parameter regression but is not a 3D arc length.
struct ZFit {
  double z0_cm     = 0.0;
  double tanLambda = 0.0;
};
ZFit fitZvsArc(const std::vector<edm4hep::TrackerHitPlane>& orderedHits,
               double cx_cm, double cy_cm, double r_cm, double charge);

}  // namespace k4rt

#endif  // TRACKING_TRACKINGPRIMITIVES_H
