#include <cmath>
#include <fstream>
#include <map>

#include "param.hh"

using namespace Geometry;

// Find the intersection in the [0.5, 1] parametric interval
Point2D intersect(const BSCurve &a, const BSCurve &b) {
  size_t res = 1000;
  Point2D result(0.5, 0.5);
  auto pa = a.eval(0.5); pa /= pa[2];
  auto pb = b.eval(0.5); pb /= pb[2];
  auto dmin = (pa - pb).norm();
  for (size_t i = 50; i <= res; ++i) {
    auto u = i / (double)res;
    auto pa = a.eval(u); pa /= pa[2];
    for (size_t j = 50; j <= res; ++j) {
      auto v = j / (double)res;
      auto pb = b.eval(v); pb /= pb[2];
      auto d = (pa - pb).norm();
      if (d < dmin) {
        dmin = d;
        auto p = (pa + pb) / 2;
        result = { p[0], p[1] };
      }
    }
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
