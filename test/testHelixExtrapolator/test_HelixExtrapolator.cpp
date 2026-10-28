// Unit test for k4rt::HelixExtrapolator — no Gaudi, no DD4hep, no input files.
//
// Every analytic crossing is checked against an INDEPENDENT reference: a 4th-order
// Runge-Kutta integration of the Lorentz force in a uniform field, starting from the same
// physical initial condition (position, momentum, charge, B). Agreement must be < 10 um.
#include "HelixExtrapolator.h"

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

using k4rt::Helix;
using Vec = Eigen::Vector3d;

namespace {
constexpr double kKappa = 0.299792458e-2;  // GeV / (T cm)
int nFail = 0;

void check(bool ok, const std::string& what, double val, double lim) {
  std::printf("  %-58s %12.3e  (limit %8.1e)  %s\n", what.c_str(), val, lim, ok ? "ok" : "FAIL");
  if (!ok) ++nFail;
}

struct State { Vec x; Vec p; };

State rk4Step(const State& st, double h, double q, double B) {
  auto f = [&](const State& s) {
    const Vec t = s.p.normalized();
    return State{t, kKappa * q * B * t.cross(Vec(0, 0, 1)) * 1.0};
  };
  // dp/dl = kappa q (t x B) with B along +z
  auto add = [](const State& a, const State& b, double k) { return State{a.x + k * b.x, a.p + k * b.p}; };
  const State k1 = f(st);
  const State k2 = f(add(st, k1, h / 2));
  const State k3 = f(add(st, k2, h / 2));
  const State k4 = f(add(st, k3, h));
  return State{st.x + h / 6 * (k1.x + 2 * k2.x + 2 * k3.x + k4.x),
               st.p + h / 6 * (k1.p + 2 * k2.p + 2 * k3.p + k4.p)};
}

// Integrate (forward if dir=+1, backward if -1) until g(x) changes sign, then bisect.
Vec rk4Until(State st, double q, double B, int dir, const std::function<double(const Vec&)>& g) {
  const double h0 = 0.05 * dir;  // cm of 3D path
  if (dir < 0) st.p = -st.p, q = -q;  // backward = forward with reversed momentum and charge
  double g0 = g(st.x);
  for (int i = 0; i < 400000; ++i) {
    State nx = rk4Step(st, std::abs(h0), q, B);
    const double g1 = g(nx.x);
    if ((g0 < 0) != (g1 < 0)) {
      double lo = 0, hi = std::abs(h0);
      for (int k = 0; k < 80; ++k) {
        const double mid = 0.5 * (lo + hi);
        if ((g(rk4Step(st, mid, q, B).x) < 0) == (g0 < 0)) lo = mid; else hi = mid;
      }
      return rk4Step(st, 0.5 * (lo + hi), q, B).x;
    }
    st = nx; g0 = g1;
  }
  return Vec(NAN, NAN, NAN);
}

// Helix for a particle at x0 with momentum p0, charge q, in B (T) along +z.
Helix helixFromPhysics(const Vec& x0, const Vec& p0, double q, double B) {
  Helix h;
  const double pT = std::hypot(p0.x(), p0.y());
  const double tx = p0.x() / pT, ty = p0.y() / pT;
  h.R   = pT / (kKappa * std::abs(q * B));
  h.rot = (q * B > 0) ? -1 : +1;  // positive charge in +Bz turns clockwise
  h.cx  = x0.x() + h.rot * h.R * (-ty);
  h.cy  = x0.y() + h.rot * h.R * (tx);
  h.phiRef = std::atan2(x0.y() - h.cy, x0.x() - h.cx);
  h.zRef = x0.z();
  h.tanL = p0.z() / pT;
  return h;
}
}  // namespace

int main() {
  struct Case { const char* name; Vec x0; Vec p0; double q; double B; };
  const std::vector<Case> cases = {
    {"prompt mu+  pT=1 GeV  B=+2T", Vec(0, 0, 0),      Vec(0.8, 0.6, 0.5),   +1, 2.0},
    {"prompt mu-  pT=1 GeV  B=+2T", Vec(0, 0, 0),      Vec(0.6, -0.8, -0.3), -1, 2.0},
    {"displaced mu+ pT=3   B=-1.7", Vec(10, 5, 3),     Vec(2.0, 2.2, 1.5),   +1, -1.7},
    {"low-pT mu-  pT=0.4  B=+2T",   Vec(1, -2, 0.5),   Vec(-0.3, 0.26, 0.2), -1, 2.0},
  };
  const double tol = 1e-3;  // cm = 10 um

  for (const auto& c : cases) {
    std::printf("%s\n", c.name);
    const Helix h = helixFromPhysics(c.x0, c.p0, c.q, c.B);
    const State st{c.x0, c.p0};

    // 1. the helix starts where it should and points the right way
    check((h.at(0) - c.x0).norm() < 1e-12, "helix passes through x0", (h.at(0) - c.x0).norm(), 1e-12);
    check((h.dirAt(0) - c.p0.normalized()).norm() < 1e-12, "direction at x0 = p0 direction",
          (h.dirAt(0) - c.p0.normalized()).norm(), 1e-12);

    // 2. forward barrel crossing vs RK4
    const double r0 = std::hypot(c.x0.x(), c.x0.y());
    const double rMaxCircle = std::hypot(h.cx, h.cy) + h.R;  // farthest the circle gets from the axis
    const double Rc = r0 + 0.5 * (rMaxCircle - r0);           // reachable by construction
    auto cyl = k4rt::crossCylinder(h, Rc, 0.0, +1);
    check(cyl.ok, "forward cylinder crossing found", cyl.ok ? 0 : 1, 0.5);
    if (cyl.ok) {
      check(std::abs(std::hypot(cyl.pos.x(), cyl.pos.y()) - Rc) < 1e-9, "crossing lies on the cylinder",
            std::abs(std::hypot(cyl.pos.x(), cyl.pos.y()) - Rc), 1e-9);
      const Vec ref = rk4Until(st, c.q, c.B, +1, [&](const Vec& x) { return std::hypot(x.x(), x.y()) - Rc; });
      check((cyl.pos - ref).norm() < tol, "forward cylinder vs RK4", (cyl.pos - ref).norm(), tol);
    }

    // 3. disk crossing vs RK4
    const double zd = c.x0.z() + (c.p0.z() > 0 ? 20.0 : -20.0);
    auto dk = k4rt::crossPlaneZ(h, zd, 0.0, +1);
    check(dk.ok, "forward disk crossing found", dk.ok ? 0 : 1, 0.5);
    if (dk.ok) {
      const Vec ref = rk4Until(st, c.q, c.B, +1, [&](const Vec& x) { return x.z() - zd; });
      check((dk.pos - ref).norm() < tol, "disk vs RK4", (dk.pos - ref).norm(), tol);
    }

    // 4. backward crossing of a smaller radius (needs a displaced start)
    if (r0 > 1.0) {
      const double Rb = 0.5 * r0;
      auto bk = k4rt::crossCylinder(h, Rb, 0.0, -1);
      if (bk.ok) {
        const Vec ref = rk4Until(st, c.q, c.B, -1, [&](const Vec& x) { return std::hypot(x.x(), x.y()) - Rb; });
        check((bk.pos - ref).norm() < tol, "backward cylinder vs RK4", (bk.pos - ref).norm(), tol);
        check(bk.s < 0, "backward crossing has s < 0", bk.s, 0);
      } else {
        std::printf("  (backward cylinder r=%.2f not reached — skipped)\n", Rb);
      }
    }

    // 5. recover the helix from points sampled on it (rot, tanL, zRef, unwrapping)
    std::vector<Vec> pts;
    for (int i = 0; i < 7; ++i) pts.push_back(h.at(0.45 * h.R * i));  // > pi of turning in total
    Helix g;
    std::vector<double> s;
    const bool okFit = k4rt::helixFromHits(h.cx, h.cy, h.R, pts, g, &s);
    check(okFit && g.rot == h.rot, "helixFromHits: rotation sense", okFit ? std::abs(g.rot - h.rot) : 9, 0.5);
    check(std::abs(g.tanL - h.tanL) < 1e-9, "helixFromHits: tanLambda", std::abs(g.tanL - h.tanL), 1e-9);
    check(std::abs(g.zRef - h.zRef) < 1e-9, "helixFromHits: zRef", std::abs(g.zRef - h.zRef), 1e-9);
    check(std::abs(s.back() - 0.45 * h.R * 6) < 1e-9, "helixFromHits: unwrapped arc length > pi R",
          std::abs(s.back() - 0.45 * h.R * 6), 1e-9);

    // 6. exact intersection with a tilted sensor plane near the cylinder crossing
    if (cyl.ok) {
      const Vec n = Vec(cyl.pos.x(), cyl.pos.y(), 0).normalized() + Vec(0.2, -0.1, 0.15);
      const Vec nn = n.normalized();
      const Vec o = cyl.pos + 0.3 * Vec(-nn.y(), nn.x(), 0);  // plane through a nearby point
      auto pl = k4rt::refineToPlane(h, o, nn, cyl.s);
      check(pl.ok, "refineToPlane converged", pl.ok ? 0 : 1, 0.5);
      if (pl.ok) {
        check(std::abs(nn.dot(pl.pos - o)) < 1e-9, "point lies on the plane", std::abs(nn.dot(pl.pos - o)), 1e-9);
        const Vec ref = rk4Until(st, c.q, c.B, +1, [&](const Vec& x) { return nn.dot(x - o); });
        check((pl.pos - ref).norm() < tol, "plane intersection vs RK4", (pl.pos - ref).norm(), tol);
      }
    }
  }
  std::printf("\n%s: %d failure(s)\n", nFail ? "FAILED" : "PASSED", nFail);
  return nFail ? 1 : 0;
}
