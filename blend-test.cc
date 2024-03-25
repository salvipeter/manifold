#include <fstream>

#include "blend.hh"

int main(int argc, char **argv) {
  size_t res = 1000;
  std::ofstream f("/tmp/blend.txt");
  for (size_t i = 0; i <= res; ++i) {
    double x = (double)i / res;
    f << x
      << ' ' << blendFunction(x, BlendType::G1)
      << ' ' << blendFunction(x, BlendType::G2)
      << ' ' << blendFunction(x, BlendType::Bump)
      << ' ' << blendFunction(x, BlendType::ERBS)
      << std::endl;
  }
}
