#include <cmath>

#include <geometry.hh>

using namespace Geometry;

#define PS_OUTPUT

#ifdef PS_OUTPUT

void init() {
  std::cout << "%!PS" << std::endl;
}
void writeSegments(const PointVector &pv) {
  std::cout << "newpath" << std::endl;
  size_t n = pv.size();
  for (size_t k = 0; k < n; ++k) {
    const auto &p = pv[k];
    auto x = p[0] * 200 + 250;
    auto y = 400 - p[1] * 200;
    std::cout << x << ' ' << y << (k ? " lineto" : " moveto") << std::endl;
  }
  std::cout << "stroke" << std::endl;
}
void finalize() {
  std::cout << "showpage" << std::endl;
}

#else // !PS_OUTPUT

void init() { }
void writeSegments(const PointVector &pv) {
  static size_t index = 1;
  size_t n = pv.size();
  for (const auto &p : pv)
    std::cout << "v " << p << std::endl;
  std::cout << 'l';
  for (size_t k = 0; k < n; ++k)
    std::cout << ' ' << index + k;
  std::cout << std::endl;
  index += n;
}
void finalize() { }

#endif

int main() {
  bool u_dir = true;
  size_t n = 5, res = 5;

  init();

  if (true) {                   // domain polygon
    PointVector pv;
    for (size_t k = 0; k <= n; ++k)
      pv.emplace_back(std::cos(2 * k * M_PI / n), std::sin(2 * k * M_PI / n), 1);
    writeSegments(pv);
  }

  if (true) {                   // bilinear mapping
    auto c = std::cos(2 * M_PI / n), s = std::sin(2 * M_PI / n);
    PointVector corners = {
      { 0, 0, 1 }, { 0.5 + c / 2, -s / 2, 1 }, { 1, 0, 1 }, { 0.5 + c / 2, s / 2, 1 }
    };
    for (size_t i = 0; i <= res; ++i) {
      PointVector pv;
      auto u = i / (double)res;
      for (size_t j = 0; j <= 100; ++j) {
        auto v = j / 100.0;
        if (!u_dir)
          std::swap(u, v);
        auto p1 = corners[0] + (corners[1] - corners[0]) * u;
        auto p2 = corners[3] + (corners[2] - corners[3]) * u;
        pv.push_back(p1 + (p2 - p1) * v);
        if (!u_dir)
          std::swap(u, v);
      }
      writeSegments(pv);
    }
    finalize();
    return 0;
  }

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
  PointVector points;
  points.push_back(u_dir ? pv[0] : pv[6]);
  points.emplace_back(0, 0, 1);
  points.push_back(u_dir ? pv[4] : pv[2]);
  writeSegments(points);
  for (size_t i = 1; i <= res; ++i) {
    auto u = i / (double)res;
    BSCurve curve;
    if (u_dir)
      curve = BSCurve({
          pv[0] * (1 - u) + pv[1] * u,
          pv[2] * exp(u * u - u),
          pv[4] * (1 - u) + pv[3] * u
        });
    else
      curve = BSCurve({
          pv[6] * (1 - u) + pv[5] * u,
          pv[4] * exp(u * u - u),
          pv[2] * (1 - u) + pv[3] * u
        });
    curve.controlPoints()[1][2] = 1 / u;
    points.clear();
    for (size_t j = 0; j <= 100; ++j) {
      auto v = j / 100.0;
      auto p = curve.eval(v);
      points.emplace_back(p[0] / p[2], p[1] / p[2], 1);
    }
    writeSegments(points);
  }

  finalize();
}
