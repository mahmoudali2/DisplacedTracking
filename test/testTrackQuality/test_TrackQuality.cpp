// Toy calibration test of the kink chi2 (TrackQuality.h): muons through 9 barrel layers with
// Highland scattering at each layer and Gaussian measurement noise. With the true radius the
// expected chi2 per kink pair is 2, so <chi2> / (2 (n-2)) must be ~1 at every momentum.
#include "HelixExtrapolator.h"
#include "TrackQuality.h"

#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

using namespace k4rt;
using V = Eigen::Vector3d;

namespace {
constexpr double kKappa = 0.299792458e-2;
Helix helixFrom(const V& x0, const V& dir, double pT, double q, double B) {
  Helix h; const double tx = dir.x() / std::hypot(dir.x(), dir.y()), ty = dir.y() / std::hypot(dir.x(), dir.y());
  h.R = pT / (kKappa * std::abs(q * B)); h.rot = (q * B > 0) ? -1 : +1;
  h.cx = x0.x() + h.rot * h.R * (-ty); h.cy = x0.y() + h.rot * h.R * tx;
  h.phiRef = std::atan2(x0.y() - h.cy, x0.x() - h.cx); h.zRef = x0.z();
  h.tanL = dir.z() / std::hypot(dir.x(), dir.y()); return h;
}
}  // namespace

int main() {
  const std::vector<double> radii = {1.4, 2.4, 3.5, 13.0, 40.0, 67.0, 100.0, 157.0, 204.0};
  const double x0 = 0.01, sig = 5e-4, B = 2.0, m = 0.10566;
  std::mt19937 rng(7); std::normal_distribution<double> g(0, 1); std::uniform_real_distribution<double> u(0, 2 * M_PI);
  const KinkConfig kc{B, m};
  int nFail = 0;
  std::printf("  p [GeV]   tracks   <chi2>/(2(n-2))   expected 1.00 +- 0.05\n");
  for (double p : {0.3, 0.5, 1.0, 2.0, 5.0, 20.0}) {
    double sum = 0; int nTrk = 0;
    for (int t = 0; t < 2000; ++t) {
      const double theta = 70.0 * M_PI / 180.0, phi = u(rng), q = (t % 2) ? 1.0 : -1.0;
      V pos(0, 0, 0), dir(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
      const double pT = p * std::sin(theta);
      std::vector<KinkPoint> pts;
      for (double R : radii) {
        const Helix h = helixFrom(pos, dir, pT, q, B);
        const Crossing c = crossCylinder(h, R, 0.0, +1);
        if (!c.ok) break;
        pos = c.pos; dir = h.dirAt(c.s);
        KinkPoint k; k.normal = V(pos.x(), pos.y(), 0).normalized(); k.xOverX0 = x0; k.sigT_cm = sig; k.sigZ_cm = sig;
        const V tphi = V(-k.normal.y(), k.normal.x(), 0);                     // measurement directions
        k.pos = pos + sig * g(rng) * tphi + sig * g(rng) * V::UnitZ();
        pts.push_back(k);
        // scattering in two planes orthogonal to the direction
        const double cosInc = std::abs(k.normal.dot(dir));
        const double th0 = highlandTheta0(x0 / cosInc, p, m);
        const V e1 = dir.cross(V::UnitZ()).normalized(), e2 = dir.cross(e1).normalized();
        dir = (dir + th0 * g(rng) * e1 + th0 * g(rng) * e2).normalized();
      }
      if (pts.size() < 4) continue;
      const Helix h0 = helixFrom(V(0, 0, 0), V(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)), pT, q, B);
      const KinkResult r = kinkChi2(pts, h0.R, h0.rot, kc);
      sum += r.chi2 / (2.0 * (pts.size() - 2)); ++nTrk;
    }
    const double mean = sum / nTrk; const bool ok = std::abs(mean - 1.0) < 0.05;
    std::printf("  %6.1f   %6d   %8.3f            %s\n", p, nTrk, mean, ok ? "ok" : "FAIL");
    nFail += !ok;
  }
  std::printf("%s: %d failure(s)\n", nFail ? "FAILED" : "PASSED", nFail);
  return nFail ? 1 : 0;
}
