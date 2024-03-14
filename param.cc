#include <cmath>
#include <fstream>
#include <limits>
#include <map>

#include "param.hh"

using namespace Geometry;

size_t binomial(size_t n, size_t k) {
  if (k > n)
    return 0;
  size_t result = 1;
  for (size_t d = 1; d <= k; ++d, --n)
    result = result * n / d;
  return result;
}

Point3D rationalEval(const BSCurve &curve, double u, size_t d, VectorVector &der) {
  der.clear(); der.reserve(d + 1);
  VectorVector der3d;
  curve.eval(u, d, der3d);
  for (size_t k = 0; k <= d; ++k) {
    Vector3D v = der3d[k];
    for (size_t i = 1; i <= k; ++i)
      v = v - der[k-i] * der3d[i][2] * binomial(k, i);
    der.push_back(v / der3d[0][2]);
  }
  return der[0];
}

// Newton iteration to find the perpendicular distance from a point
// to the curve (in parametric interval [0.5, 1])
double newton(const BSCurve &curve, const Point3D &p) {
  size_t max_iteration = 10, res = 20;
  double distance_tol = 1e-8, cosine_tol = 1e-8;

  VectorVector der;
  auto lo = 0.5;
  auto hi = 1.0;

  // Find initial value
  double u = 0, min = std::numeric_limits<double>::max();
  for (size_t i = res / 2; i <= res; ++i) {
    auto t = (double)i / res;
    auto q = curve.eval(t); q /= q[2];
    auto d = (p - q).norm();
    if (d < min) {
      u = t;
      min = d;
    }
  }

  Vector3D deviation;
  double distance;
  for (size_t iteration = 0; iteration < max_iteration; ++iteration) {
    deviation = rationalEval(curve, u, 2, der) - p;
    distance = deviation.norm();
    if (distance < distance_tol)
      break;

    double scaled_error = der[1] * deviation;
    double cosine_err = std::abs(scaled_error) / (der[1].norm() * distance);
    if (cosine_err < cosine_tol)
      break;

    double old = u;
    u -= scaled_error / (der[2] * deviation + der[1] * der[1]);
    u = std::min(std::max(u, lo), hi);

    if ((der[1] * (u - old)).norm() < distance_tol)
      break;
  }

  return distance;
}

// Find the intersection in the [0.5, 1] parametric interval
Point2D intersect(const BSCurve &a, const BSCurve &b) {
  size_t iterations = 20;
  double tolerance = 1e-8;
  Point2D result;
  double lo = 0.5, hi = 1.0, phi = (std::sqrt(5) + 1) / 2;
  double d = (phi - 1) * (hi - lo);
  auto u1 = lo + d, u2 = hi - d;
  auto p1 = a.eval(u1); p1 /= p1[2];
  auto d1 = newton(b, p1);
  auto p2 = a.eval(u2); p2 /= p2[2];
  auto d2 = newton(b, p2);
  for (size_t i = 0; i < iterations; ++i) {
    double u;
    if (d1 < d2) {
      u = u1;
      result = { p1[0], p1[1] };
      lo = u2;
      double tmp = u1;
      u1 = u2 + (phi - 1) * (hi - u2);
      u2 = tmp;
      p2 = p1; d2 = d1;
      p1 = a.eval(u1); p1 /= p1[2];
      d1 = newton(b, p1);
    } else {
      u = u2;
      result = { p2[0], p2[1] };
      hi = u1;
      double tmp = u2;
      u2 = u1 - (phi - 1) * (u1 - lo);
      u1 = tmp;
      p1 = p2; d1 = d2;
      p2 = a.eval(u2); p2 /= p2[2];
      d2 = newton(b, p2);
    }
    if (u != 0 && (2 - phi) * std::abs((hi - lo) / u) < tolerance)
      break;
  }
  return result;
}

struct Compare {
  bool operator()(const std::pair<size_t, Point2D> &a,
                  const std::pair<size_t, Point2D> &b) const {
    return a.first < b.first ||
      (a.first == b.first &&
       (a.second[0] < b.second[0] ||
        (a.second[0] == b.second[0] && a.second[1] < b.second[1])));
  }
};

std::map<std::pair<size_t, Point2D>, Point2D, Compare> param_cache;

bool loadCache(std::string filename) {
  param_cache.clear();
  std::ifstream f(filename, std::ios::binary);
  if (!f.good())
    return false;
  f.exceptions(std::ios::failbit | std::ios::badbit);
  size_t n;
  f.read(reinterpret_cast<char *>(&n), sizeof(n));
  for (size_t i = 0; i < n; ++i) {
    std::pair<size_t, Point2D> k;
    Point2D v;
    f.read(reinterpret_cast<char *>(&k), sizeof(k));
    f.read(reinterpret_cast<char *>(&v), sizeof(v));
    param_cache[k] = v;
  }
  return true;
}

bool saveCache(std::string filename) {
  std::ofstream f(filename, std::ios::binary);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  size_t n = param_cache.size();
  f.write(reinterpret_cast<char *>(&n), sizeof(n));
  for (const auto &[k, v] : param_cache) {
    f.write(reinterpret_cast<const char *>(&k), sizeof(k));
    f.write(reinterpret_cast<const char *>(&v), sizeof(v));
  }
  return true;
}

Point2D param(size_t n, const Point2D &uv) {
  if (param_cache.contains({n, uv}))
    return param_cache[{n, uv}];

  if (false) {
    auto c = std::cos(2 * M_PI / n), s = std::sin(2 * M_PI / n);
    Point2DVector corners = {
      { 0, 0 }, { 0.5 + c / 2, -s / 2 }, { 1, 0 }, { 0.5 + c / 2, s / 2 }
    };
    auto p1 = corners[0] + (corners[1] - corners[0]) * uv[0];
    auto p2 = corners[3] + (corners[2] - corners[3]) * uv[0];
    auto result = p1 + (p2 - p1) * uv[1];
    param_cache[{n, uv}] = result;
    return result;
  }

  auto u = uv[0], v = uv[1];
  if (u == 0.0 && v == 0.0)
    return { 0, 0 };

  auto c = std::cos(2 * M_PI / n), s = std::sin(2 * M_PI / n);
  auto c2 = std::cos(4 * M_PI / n), s2 = std::sin(4 * M_PI / n);
  PointVector pv = {
    { (c2 + c) / 2, -(s2 + s) / 2, 1 },
    { c, -s, 1 },
    { 0.5 + c / 2, -s / 2, 1 },
    { 1, 0, 1 },
    { 0.5 + c / 2, s / 2, 1 },
    { c, s, 1 },
    { (c2 + c) / 2, (s2 + s) / 2, 1 }
  };
  BSCurve ucurve({
      pv[0] * (1 - u) + pv[1] * u,
      pv[2] * exp(u * u - u),
      pv[4] * (1 - u) + pv[3] * u
    });
  ucurve.controlPoints()[1][2] = 1 / u;
  BSCurve vcurve({
      pv[6] * (1 - v) + pv[5] * v,
      pv[4] * exp(v * v - v),
      pv[2] * (1 - v) + pv[3] * v
    });
  vcurve.controlPoints()[1][2] = 1 / v;

  Point2D result;
  if (u == 0)
    result = intersect(BSCurve({{-pv[4][0], -pv[4][1], 1}, {0, 0, 1}, pv[4]}), vcurve);
  else if (v == 0)
    result = intersect(BSCurve({{-pv[2][0], -pv[2][1], 1}, {0, 0, 1}, pv[2]}), ucurve);
  else {
    result = intersect(ucurve, vcurve);
  }
  param_cache[{n, uv}] = result;
  return result;
}
