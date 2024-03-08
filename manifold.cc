#include <fstream>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <transfinite/domain.hh>
#include <transfinite/surface-c0coons.hh>
#include <transfinite/surface-generalized-bezier.hh>
#include <transfinite/surface-midpoint.hh>

using namespace Geometry;

struct Traits : public OpenMesh::DefaultTraits {
  using Point  = OpenMesh::Vec3d;
};
using Mesh = OpenMesh::PolyMesh_ArrayKernelT<Traits>;

class NPatch {
public:
  virtual ~NPatch() {}
  virtual Point3D eval(const Point2D &) const = 0;
  virtual size_t size() const = 0;
};

class MP : public NPatch {
public:
  MP(const Transfinite::SurfaceMidpoint &mp) : patch(mp) {}
  ~MP() {}
  Point3D eval(const Point2D &p) const override { return patch.eval(p); }
  size_t size() const override { return patch.domain()->size(); }
private:
  Transfinite::SurfaceMidpoint patch;
};

// Quadratic Bezier patch in [-1,1]^2
class QB : public NPatch {
public:
  QB(const PointMatrix &curves, const Point3D &center) {
    const auto &P00 = curves[2][2];
    const auto &P10 = curves[3][1];
    const auto &P20 = curves[0][0];
    const auto &P01 = curves[2][1];
    const auto &P21 = curves[0][1];
    const auto &P02 = curves[1][2];
    const auto &P12 = curves[1][1];
    const auto &P22 = curves[0][2];
    patch = BSSurface(2, 2, {
        P00, P01, P02,
        P10, (center*16-(P00+P02+P20+P22)-(P10+P01+P12+P21)*2)/4, P12,
        P20, P21, P22
      });
  }
  ~QB() {}
  Point3D eval(const Point2D &p) const override {
    return patch.eval(0.5 + p[0] / 2, 0.5 + p[1] / 2);
  }
  size_t size() const override { return 4; }
private:
  BSSurface patch;
};

// Generates the N-patch based on the quads around
// the vertex at the origin of the given half-edge.
std::unique_ptr<NPatch> generateNPatch(const Mesh &mesh, OpenMesh::SmartHalfedgeHandle he) {
  Point3D center(mesh.point(he.from()).data());
  PointMatrix curves;
  auto start = he.next(), it = start;
  do {
    Point3D p0(mesh.point(it.to()).data());
    it = it.next();
    Point3D p1(mesh.point(it.to()).data());
    it = it.next().opp().next();
    Point3D p2(mesh.point(it.to()).data());
    curves.push_back({ p0, p1 * 2 - (p0 + p2) / 2, p2 });
  } while (it != start);
  curves.insert(curves.begin(), curves.back());
  curves.pop_back();

  if (curves.size() == 4)
    return std::make_unique<QB>(curves, center);

  Transfinite::SurfaceMidpoint surface;
  std::vector<std::shared_ptr<Transfinite::Curve>> boundaries;
  std::transform(curves.begin(), curves.end(), std::back_inserter(boundaries),
                 [](const PointVector &p) {
                   return std::make_shared<Transfinite::BSplineCurve>(p);
                 });
  surface.setCurves(boundaries);

  // GB patch specific part
  // double a = 1.0/3.0, b = 2.0/3.0;
  // surface.initNetwork(n, 3);

  // center = Point3D(0, 0, 0);
  // for (size_t i = 0; i < n; ++i)
  //   center += curves[i][1];
  // center /= n;

  // PointMatrix opposite;
  // size_t n = curves.size();
  // for (size_t i = 0; i < n; ++i) {
  //   size_t im = (i + n - 1) % n, ip = (i + 1) % n;
  //   auto q0 = curves[im][1];
  //   auto q1 = center;
  //   auto q2 = curves[ip][1];
  //   opposite.push_back({ q0, q1, q2 });
  // }

  // surface.setCentralControlPoint(center);
  // for (size_t i = 0; i < n; ++i) {
  //   const auto &p = curves[i];
  //   const auto &q = opposite[i];
  //   surface.setControlPoint(i, 0, 0, p[0]);
  //   surface.setControlPoint(i, 1, 0, p[0] * a + p[1] * b);
  //   surface.setControlPoint(i, 2, 0, p[2] * a + p[1] * b);
  //   surface.setControlPoint(i, 3, 0, p[2]);
  //   surface.setControlPoint(i, 0, 1, p[0] * a + q[0] * b);
  //   surface.setControlPoint(i, 1, 1, p[0] * a * a + (p[1] + q[0]) * a * b + q[1] * b * b);
  //   surface.setControlPoint(i, 2, 1, p[2] * a * a + (p[1] + q[2]) * a * b + q[1] * b * b);
  //   surface.setControlPoint(i, 3, 1, p[2] * a + q[2] * b);
  // }

  surface.setupLoop();
  surface.update();
  surface.setMidpoint(center);
  return std::make_unique<MP>(surface);
}

double trapezoid(const std::function<double(double)> &f, double a, double b,
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

double simpson(const std::function<double(double)> &f, double a, double b,
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

double erbsBlend(double t) {
  if (t < 1e-5)
    return 0;
  size_t iterations = 20;
  constexpr double Sd = 1.6571376796460222;
  auto phi = [](double s) { return std::exp(-std::pow(s - 0.5, 2) / (s * (1 - s))); };
  return Sd * simpson(phi, 0, t, iterations);
}

// f(0) = 1, f(1) = 0, and f^k(0) = f^k(1) = 0 for some k > 0
double blendFunction(double x) {
  // Linear
  // return 1 - x;
  // G1 Hermite
  // return std::pow(1 - x, 3) + 3 * std::pow(1 - x, 2) * x;
  // G2 Hermite
  // return std::pow(1 - x, 5) + 5 * std::pow(1 - x, 4) * x + 10 * std::pow(1 - x, 3) * x * x;
  // Bump
  // return std::exp(-1 / (1 - x)) / (std::exp(-1 / x) + std::exp(-1 / (1 - x)));
  // ERBS
  return 1 - erbsBlend(x);
}

double blend(const Point2D &uv) {
  return blendFunction(uv[0]) * blendFunction(uv[1]);
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input-mesh> [resolution]" << std::endl;
    return 1;
  }

  Mesh cage;
  if (!OpenMesh::IO::read_mesh(cage, argv[1]))
    return 2;

  size_t resolution = 50;
  if (argc == 3)
    resolution = std::atoi(argv[2]);

  TriMesh mesh;

  for (auto face : cage.faces()) {
    std::vector<std::unique_ptr<NPatch>> patches;
    for (auto he : face.halfedges())
      patches.push_back(generateNPatch(cage, he));

    PointVector points;
    for (size_t j = 0; j < resolution; ++j) {
      double u = (double)j / (resolution - 1);
      for (size_t k = 0; k < resolution; ++k) {
        double v = (double)k / (resolution - 1);
        Point3D p(0, 0, 0);
        for (size_t i = 0; i < 4; ++i) {
          const auto &surface = patches[i];
          Point2D uv;
          switch (i) {
          case 0: uv[0] = u; uv[1] = v; break;
          case 1: uv[0] = v; uv[1] = 1 - u; break;
          case 2: uv[0] = 1 - u; uv[1] = 1 - v; break;
          case 3: uv[0] = 1 - v; uv[1] = u; break;
          }
          auto b = blend(uv);
          auto n = surface->size();
          if (n != 4) {
            auto c = std::cos(2 * M_PI / n), s = std::sin(2 * M_PI / n);
            Point2DVector corners = {
              { 0, 0 }, { 0.5 + c / 2, -s / 2 }, { 1, 0 }, { 0.5 + c / 2, s / 2 }
            };
            auto p1 = corners[0] + (corners[1] - corners[0]) * uv[0];
            auto p2 = corners[3] + (corners[2] - corners[3]) * uv[0];
            uv = p1 + (p2 - p1) * uv[1];
          }
          p += surface->eval(uv) * b;
        }
        points.push_back(p);
      }
    }

    TriMesh localmesh;
    localmesh.setPoints(points);
    for (size_t j = 1; j < resolution; ++j)
      for (size_t k = 1; k < resolution; ++k) {
        size_t index = j * resolution + k;
        localmesh.addTriangle(index - resolution - 1, index, index - resolution);
        localmesh.addTriangle(index - 1, index, index - resolution - 1);
      }
    mesh.append(localmesh);
  }
  mesh.writeSTL("/tmp/surface.stl");
}
