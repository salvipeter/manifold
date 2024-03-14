#include <fstream>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/CatmullClarkT.hh>

#include <qgb.hh>

#include "blend.hh"
#include "param.hh"

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

class QGBWrap : public NPatch {
public:
  QGBWrap(const QGB &qgb) : patch(qgb) {}
  ~QGBWrap() {}
  Point3D eval(const Point2D &p) const override { return patch.eval(p); }
  size_t size() const override { return patch.size(); }
private:
  QGB patch;
};

class CBT : public NPatch {
public:
  CBT(const std::vector<QGB::Boundary> &curves, const Point3D &center) {
    const auto &c = curves;
    cpts = PointVector{
      c[0][0], (c[0][0] + c[0][1] * 2) / 3, (c[0][2] + c[0][1] * 2) / 3, c[0][2],
      (c[0][0] + c[2][1] * 2) / 3, { 0, 0, 0 }, (c[0][2] + c[1][1] * 2) / 3,
      (c[2][0] + c[2][1] * 2) / 3, (c[1][2] + c[1][1] * 2) / 3,
      c[1][2]
    };
    cpts[5] = (center * 27 - (cpts[0] + cpts[3] + cpts[9]) -
               (cpts[1] + cpts[2] + cpts[4] + cpts[6] + cpts[7] + cpts[8]) * 3) / 6;
  }
  ~CBT() {}
  Point3D eval(const Point2D &p) const override {
    double l1 = (1 + p[0] * 2) / 3;
    double l2 = (1 - p[0] + p[1] * std::sqrt(3)) / 3;
    double l3 = 1 - l1 - l2;
    return
      cpts[0] * l1 * l1 * l1 +
      cpts[1] * 3 * l1 * l1 * l2 +
      cpts[2] * 3 * l1 * l2 * l2 +
      cpts[3] * l2 * l2 * l2 +
      cpts[4] * 3 * l1 * l1 * l3 +
      cpts[5] * 6 * l1 * l2 * l3 +
      cpts[6] * 3 * l2 * l2 * l3 +
      cpts[7] * 3 * l1 * l3 * l3 +
      cpts[8] * 3 * l2 * l3 * l3 +
      cpts[9] * l3 * l3 * l3;
  }
  size_t size() const override { return 3; }
private:
  PointVector cpts;
};

// Generates the N-patch based on the quads around
// the vertex at the origin of the given half-edge.
std::unique_ptr<NPatch> generateNPatch(const Mesh &mesh, OpenMesh::SmartHalfedgeHandle he) {
  Point3D center(mesh.point(he.from()).data());
  std::vector<QGB::Boundary> curves;
  auto start = he.next(), it = start;
  do {
    Point3D p0(mesh.point(it.to()).data());
    it = it.next();
    Point3D p1(mesh.point(it.to()).data());
    it = it.next().opp().next();
    Point3D p2(mesh.point(it.to()).data());
    curves.push_back({ p0, p1 * 2 - (p0 + p2) / 2, p2 });
  } while (it != start);

  auto n = curves.size();
  if (n == 3)
    return std::make_unique<CBT>(curves, center);

  QGB surface(n);
  for (size_t i = 0; i < n; ++i)
    surface.setBoundary(i, curves[(i+n-1)%n]);
  surface.setMidpoint(center);
  return std::make_unique<QGBWrap>(surface);
}

double blend(const Point2D &uv) {
  auto type = BlendType::G2;
  return blendFunction(uv[0], type) * blendFunction(uv[1], type);
}

int main(int argc, char **argv) {
  if (argc != 2 && argc != 3) {
    std::cerr << "Usage: " << argv[0] << " <input-mesh> [resolution]" << std::endl;
    return 1;
  }

  Mesh cage;
  if (!OpenMesh::IO::read_mesh(cage, argv[1]))
    return 2;

  // Check if it is a quad mesh - if not, perform 1 Catmull-Clark subdivision step
  bool only_quads = true;
  for (auto face : cage.faces()) {
    size_t nf = 0;
    for (auto it = cage.cfv_iter(face); it.is_valid(); ++it)
      nf++;
    if (nf != 4) {
      only_quads = false;
      break;
    }
  }
  if (!only_quads) {
    OpenMesh::Subdivider::Uniform::CatmullClarkT<Mesh> cc;
    cc.attach(cage);
    cc(1);
    cc.detach();
  }

  size_t resolution = 50;
  if (argc == 3)
    resolution = std::atoi(argv[2]);
  auto cachefile = std::string("cache") + std::to_string(resolution) + ".dat";
  loadCache(cachefile);

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
          if (n != 4)
            uv = param(n, uv);
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

  saveCache(cachefile);
}
