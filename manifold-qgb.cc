#include <fstream>

#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>

#include <qgb.hh>

#include "blend.hh"
#include "param.hh"

using namespace Geometry;

struct Traits : public OpenMesh::DefaultTraits {
  using Point  = OpenMesh::Vec3d;
};
using Mesh = OpenMesh::PolyMesh_ArrayKernelT<Traits>;

// Generates the N-patch based on the quads around
// the vertex at the origin of the given half-edge.
std::unique_ptr<QGB> generateNPatch(const Mesh &mesh, OpenMesh::SmartHalfedgeHandle he) {
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
  auto surface = std::make_unique<QGB>(n);
  for (size_t i = 0; i < n; ++i)
    surface->setBoundary(i, curves[(i+n-1)%n]);
  surface->setMidpoint(center);
  return surface;
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

  size_t resolution = 50;
  if (argc == 3)
    resolution = std::atoi(argv[2]);
  auto cachefile = std::string("cache") + std::to_string(resolution) + ".dat";
  loadCache(cachefile);

  TriMesh mesh;

  for (auto face : cage.faces()) {
    std::vector<std::unique_ptr<QGB>> patches;
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
