#pragma once
// Scene pieces: the manifest written by the scene preparation tools and the scene <-> library
// mapping.
//
// Piece vertices are exported in WORLD space (scene units, meters, Z up), with all
// object and armature transforms already applied; the manifest carries no vertex
// transform. In the library (mm) a piece lives at
//     p_lib = sigma * (p_scene - c_p),
// with one global scale sigma [mm per scene unit] for all pieces and the piece's own
// area centroid c_p (re-centring keeps CurveSet::extent, and so every allowance,
// at the piece's size). The inverse is p_scene = p_lib / sigma + c_p.
#include <string>
#include <vector>

#include "model/coordinates.h"
#include "model/mesh.h"

namespace knit {

struct PieceEntry {
    std::string file, object, material;
    int island = 0, rings = 0, segments = 0;
    double area = 0;
    Vec3 centroid, poleA, poleB;  // scene units (as exported)
};

struct ExtrasEntry {
    std::string file, object;
    int islands = 0;
};

struct Manifest {
    std::string dir;  // directory of manifest.txt (piece files are relative to it)
    std::vector<PieceEntry> pieces;
    std::vector<ExtrasEntry> extras;
};

// Reads v / vt / f (v, v/vt, v/vt/vn, v//vn), polygons as triangle fans, and computes the normals.
// Returns false with a message on error.
bool readObj(const std::string& path, TriMesh& m, std::string* err = nullptr);

bool readManifest(const std::string& path, Manifest& m, std::string* err = nullptr);
std::string piecePath(const Manifest& m, const std::string& file);

// Gauge: sigma = stitches * w / L, with L the longest chart ring (scene units) of the
// gauge piece. Returns false if the piece is not in the manifest or cannot be read.
struct Gauge {
    std::string file;
    double stitches = 0, w = 0;
    double ringLength = 0, aAt = 0;  // scene units
    double sigma = 0;                // mm per scene unit
};
struct GaugeRequest {
    std::string file;  // the gauge piece
    double stitches, w;
};
bool computeGauge(const Manifest& m, const GaugeRequest& request, Gauge& g, std::string* err = nullptr);

// A piece in library units (mm, re-centred at its area centroid).
struct LoadedPiece {
    PieceEntry entry;
    SceneMapping map;
    TriMesh mesh;
};
bool loadPiece(const Manifest& m, int index, double sigma, LoadedPiece& out, std::string* err = nullptr);

} // namespace knit
