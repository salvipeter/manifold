#include <cmath>

#include <geometry.hh>

using namespace Geometry;

int main() {
  bool u_dir = true;
  size_t n = 5, res = 10;

  size_t index = 1;

  if (false) {                   // domain polygon
    for (size_t k = 0; k <= n; ++k) {
      Point3D p(std::cos(2 * k * M_PI / n), std::sin(2 * k * M_PI / n), 1);
      std::cout << "v " << p << std::endl;
    }
    std::cout << 'l';
    for (size_t k = 0; k <= n; ++k)
      std::cout << ' ' << index + k;
    std::cout << std::endl;
    index += n + 1;
  }

  if (false) {                   // bilinear mapping
    auto c = std::cos(2 * M_PI / n), s = std::sin(2 * M_PI / n);
    PointVector corners = {
      { 0, 0, 1 }, { 0.5 + c / 2, -s / 2, 1 }, { 1, 0, 1 }, { 0.5 + c / 2, s / 2, 1 }
    };
    for (size_t i = 0; i <= res; ++i) {
      auto u = i / (double)res;
      for (size_t j = 0; j <= 100; ++j) {
        auto v = j / 100.0;
        if (!u_dir)
          std::swap(u, v);
        auto p1 = corners[0] + (corners[1] - corners[0]) * u;
        auto p2 = corners[3] + (corners[2] - corners[3]) * u;
        auto p = p1 + (p2 - p1) * v;
        if (!u_dir)
          std::swap(u, v);
        std::cout << "v " << p << std::endl;
      }
      std::cout << 'l';
      for (size_t j = 0; j <= 100; ++j)
        std::cout << ' ' << index + j;
      std::cout << std::endl;
      index += 101;
    }
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
  if (u_dir)
    std::cout << "v " << pv[0] << std::endl;
  else
    std::cout << "v " << pv[6] << std::endl;
  std::cout << "v 0 0 1" << std::endl;
  if (u_dir)
    std::cout << "v " << pv[4] << std::endl;
  else
    std::cout << "v " << pv[2] << std::endl;
  std::cout << "l 1 2 3" << std::endl;
  index += 3;
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
    for (size_t j = 0; j <= 100; ++j) {
      auto v = j / 100.0;
      auto p = curve.eval(v);
      std::cout << "v " << p[0] / p[2] << ' ' << p[1] / p[2] << " 1" << std::endl;
    }
    std::cout << 'l';
    for (size_t j = 0; j <= 100; ++j)
      std::cout << ' ' << index + j;
    std::cout << std::endl;
    index += 101;
  }
}
