#include <cmath>
#include <functional>

#include "blend.hh"

static double trapezoid(const std::function<double(double)> &f, double a, double b,
                        size_t n, double s) {
  if (n == 1)
    return (f(a) + f(b)) / 2.0 * (b - a);
  double k = std::pow(2, n - 2);
  double h = (b - a) / k;
  double sum = 0;
  for (double x = a + h / 2.0; x <= b; x += h)
    sum += f(x);
  return (s + h * sum) / 2.0;
}

static double simpson(const std::function<double(double)> &f, double a, double b,
                      size_t iterations = 100, double epsilon = 1.0e-7) {
  double s = 0.0, s_prev = std::numeric_limits<double>::lowest(), st_prev = s_prev;
  for (size_t i = 1; i <= iterations; ++i) {
    double st = trapezoid(f, a, b, i, st_prev);
    s = (4.0 * st - st_prev) / 3.0;
    if (i > 5 && (std::abs(s - s_prev) < epsilon * std::abs(s_prev) || (s == 0 && s_prev == 0)))
      break;
    s_prev = s;
    st_prev = st;
  }
  return s;
}

static double erbsBlend(double t) {
  if (t < 1e-5)
    return 0;
  size_t iterations = 20;
  static constexpr double Sd = 1.6571376796460222;
  auto phi = [](double s) { return std::exp(-std::pow(s - 0.5, 2) / (s * (1 - s))); };
  return Sd * simpson(phi, 0, t, iterations);
}

// f(0) = 1, f(1) = 0, and f^k(0) = f^k(1) = 0 for some k > 0
double blendFunction(double x, BlendType type) {
  switch(type) {
  case BlendType::Linear:
    return 1 - x;
  case BlendType::G1:
    return std::pow(1 - x, 3) + 3 * std::pow(1 - x, 2) * x;
  case BlendType::G2:
    return std::pow(1 - x, 5) + 5 * std::pow(1 - x, 4) * x + 10 * std::pow(1 - x, 3) * x * x;
  case BlendType::Bump:
    return std::exp(-1 / (1 - x)) / (std::exp(-1 / x) + std::exp(-1 / (1 - x)));
  case BlendType::ERBS:
    return 1 - erbsBlend(x);
  }
  throw "Illegal blend type";
}
