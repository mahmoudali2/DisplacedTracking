// TrackQuality — a multiple-scattering-aware chi2 for track finding, built from kinks.
//
// The circle fit weights hits by their measurement error only (micrometres). At pT ~ 0.5 GeV
// multiple scattering moves hits by 0.1-1 mm, so real low-pT tracks fail any chi2 cut tuned on
// it and extension stops part-way. Inflating every hit's error with the scattering of all
// upstream material over-counts by ~3x, because those displacements are correlated.
//
// Here the scattering enters where it is independent: as KINKS. For three consecutive hits
// a-b-c the change of direction at b, in the bending plane and in dip, is compared with what
// the fitted helix predicts. Each kink has an expected spread
//     sigma^2 = theta0(b)^2  (Highland, x/X0 of b's layer / cos(incidence), at the track's p)
//             + the measurement errors of a, b, c propagated to the angle.
// Kinks at different layers are statistically independent, so
//     chi2 = sum over triplets of (dphi/sigma_phi)^2 + (dlambda/sigma_lambda)^2,
//     ndf  = 2 (n - 2) - 1          (the radius is taken from the fit)
// has chi2/ndf ~ 1 for a real track at any momentum. This is the idea behind General Broken
// Lines, used here only as a gate and a ranking, not for the final track parameters.
//
// Units: cm, GeV, tesla.
#ifndef TRACKING_TRACKQUALITY_H
#define TRACKING_TRACKQUALITY_H

#include <Eigen/Dense>

#include <vector>

namespace k4rt {

struct KinkConfig {
  double bField_T = 2.0;
  double mass_GeV = 0.10566;  ///< mass hypothesis for beta in Highland
};

/// One hit as the kink chi2 sees it.
struct KinkPoint {
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();     ///< cm
  double sigT_cm = 1e-3;                             ///< measurement error, transverse
  double sigZ_cm = 1e-3;                             ///< measurement error, along z / dip
  double xOverX0 = 0.0;                              ///< material of the hit's layer at normal incidence
  Eigen::Vector3d normal = Eigen::Vector3d::UnitX(); ///< sensor normal, for the incidence angle
};

struct KinkResult {
  double chi2 = 0.0;
  int    ndf = 0;
  std::vector<double> perTriplet;  ///< chi2 of the triplet centred at hit i (0 at both ends)
  double chi2ndf() const { return ndf > 0 ? chi2 / ndf : 0.0; }
};

/// Highland multiple-scattering angle (plane), with the log term.
double highlandTheta0(double xOverX0, double p_GeV, double mass_GeV);

/// Kink chi2 (2 dof) of the triplet a-b-c for a helix of radius R_cm and rotation sense rot
/// (+1 = anticlockwise seen from +z, the order of travel being a -> b -> c).
/// sigmaR_cm: uncertainty of R, propagated to the expected turn. Pass it when R comes from hits
/// other than a, b, c (e.g. testing a new hit against the current track's curvature); leave it
/// 0 when R was fitted to these hits.
double tripletKinkChi2(const KinkPoint& a, const KinkPoint& b, const KinkPoint& c,
                       double R_cm, int rot, const KinkConfig& cfg, double sigmaR_cm = 0.0);

/// Kink chi2 of hits ordered along the track.
KinkResult kinkChi2(const std::vector<KinkPoint>& pts, double R_cm, int rot, const KinkConfig& cfg);

}  // namespace k4rt

#endif  // TRACKING_TRACKQUALITY_H
