#include <cmath>
#include <fstream>
#include <functional>
#include <map>

#include "param.hh"

using namespace Geometry;

double golden(const std::function<double(double)> &f, double xl, double xh,
              size_t iterations = 20, double tolerance = 1e-8) {
  double x = 0;
  static constexpr double phi = (1 + std::sqrt(5.0)) / 2;
  auto d = (phi - 1) * (xh - xl);
  auto x1 = xl + d, x2 = xh - d;
  auto f1 = f(x1), f2 = f(x2);
  for (size_t i = 0; i < iterations; ++i) {
    if (f1 < f2) {
      x = x1;
      xl = x2;
      auto tmp = x1;
      x1 = x2 + (phi - 1) * (xh - x2);
      x2 = tmp;
      f2 = f1;
      f1 = f(x1);
    } else {
      x = x2;
      xh = x1;
      auto tmp = x2;
      x2 = x1 - (phi - 1) * (x1 - xl);
      x1 = tmp;
      f1 = f2;
      f2 = f(x2);
    }
    if (x != 0 && (2 - phi) * std::abs((xh - xl) / x) < tolerance)
      break;
  }
  return x;
}

// Find the intersection in the [0.5, 1] parametric interval
Point2D intersect(const BSCurve &a, const BSCurve &b) {
  auto f = [&](double u) {
    auto p = a.eval(u); p /= p[2];
    auto g = [&](double v) {
      auto q = b.eval(v); q /= q[2];
      return (p - q).norm();
    };
    auto v = golden(g, 0.5, 1.0);
    return g(v);
  };
  auto u = golden(f, 0.5, 1.0);
  auto p = a.eval(u);
  return { p[0] / p[2], p[1] / p[2] };
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
