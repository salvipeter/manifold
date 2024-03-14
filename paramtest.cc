#include <cmath>

#include <geometry.hh>

using namespace Geometry;

int main() {
  size_t n = 5;
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
  std::cout << "v " << pv[0] << std::endl;
  std::cout << "v 0 0 1" << std::endl;
  std::cout << "v " << pv[4] << std::endl;
  std::cout << "l 1 2 3" << std::endl;
  size_t index = 4;
  for (size_t i = 1; i <= 50; ++i) {
    auto u = i / 50.0;
    BSCurve curve({
        pv[0] * (1 - u) + pv[1] * u,
        pv[2] * exp(u * u - u),
        pv[4] * (1 - u) + pv[3] * u
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
