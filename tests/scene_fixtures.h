#pragma once
// Meshes, GLB files and tiled knitting parameters the scene and export test programs build by hand
// (knit_tests_s6, knit_tests_export).
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "knitting/knit/mesh_knit.h"
#include "model/mesh.h"
#include "model/vec.h"

namespace knit {

// Writes a .gltf with an external .bin: one triangle (0,0,0), (1,0,0), (0,1,0) with uvs,
// instanced by a child node under a matrix parent and by a mirrored root node.
inline std::string writeTinyGltf(const std::string& dir) {
    namespace fs = std::filesystem;
    fs::create_directories(dir);
    const float pos[9] = {0, 0, 0, 1, 0, 0, 0, 1, 0};
    const float uv[6] = {0, 0, 1, 0, 0, 1};
    const uint16_t idx[3] = {0, 1, 2};
    {
        std::ofstream b(dir + "/tiny.bin", std::ios::binary);
        b.write(reinterpret_cast<const char*>(pos), sizeof pos);
        b.write(reinterpret_cast<const char*>(uv), sizeof uv);
        b.write(reinterpret_cast<const char*>(idx), sizeof idx);
        const char pad[2] = {0, 0};
        b.write(pad, 2);
    }
    // Parent matrix (column-major): rotation by +90 deg about Y, then translation (1, 2, 3).
    // Child TRS: translation (0, 1, 0), scale 2. Mirror: scale (-1, 1, 1), translation (5, 0, 0).
    std::ofstream f(dir + "/tiny.gltf");
    f << R"({"asset":{"version":"2.0"},"scene":0,"scenes":[{"nodes":[0,2]}],
"nodes":[{"name":"parent","matrix":[0,0,-1,0, 0,1,0,0, 1,0,0,0, 1,2,3,1],"children":[1]},
 {"name":"child","translation":[0,1,0],"scale":[2,2,2],"mesh":0},
 {"name":"mirror","translation":[5,0,0],"scale":[-1,1,1],"mesh":0}],
"meshes":[{"name":"tri","primitives":[{"attributes":{"POSITION":0,"TEXCOORD_0":1},"indices":2,"material":0}]}],
"materials":[{"name":"wool","pbrMetallicRoughness":{"baseColorFactor":[0.5,0.25,0.125,1]}}],
"buffers":[{"uri":"tiny.bin","byteLength":68}],
"bufferViews":[{"buffer":0,"byteOffset":0,"byteLength":36},{"buffer":0,"byteOffset":36,"byteLength":24},{"buffer":0,"byteOffset":60,"byteLength":6}],
"accessors":[{"bufferView":0,"componentType":5126,"count":3,"type":"VEC3","min":[0,0,0],"max":[1,1,0]},
 {"bufferView":1,"componentType":5126,"count":3,"type":"VEC2"},
 {"bufferView":2,"componentType":5123,"count":3,"type":"SCALAR"}]})";
    return dir + "/tiny.gltf";
}

// A triangulated grid over [0, nx c] x [0, ny c] (two triangles per cell, CCW seen from
// the side the position map's normal points to), positions and corner uvs from maps of
// the grid point; `keep(x, y)` drops cells (holes).
inline TriMesh gridMesh(int nx, int ny, double c, const std::function<Vec3(double, double)>& pos,
                 const std::function<Vec2(double, double)>& uv, const std::function<bool(double, double)>& keep = nullptr) {
    TriMesh m;
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) m.pos.push_back(pos(i * c, j * c));
    auto id = [&](int i, int j) { return j * (nx + 1) + i; };
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
            const double cx = (i + 0.5) * c, cy = (j + 0.5) * c;
            if (keep && !keep(cx, cy)) continue;
            const int a = id(i, j), b = id(i + 1, j), d = id(i + 1, j + 1), e = id(i, j + 1);
            auto corner = [&](int v) { return uv((v % (nx + 1)) * c, (v / (nx + 1)) * c); };
            m.tri.push_back({a, b, d});
            m.uv.push_back({corner(a), corner(b), corner(d)});
            m.tri.push_back({a, d, e});
            m.uv.push_back({corner(a), corner(d), corner(e)});
        }
    // Drop unused vertices (holes) so the mesh stays compact.
    std::vector<int> used(m.pos.size(), -1);
    TriMesh out;
    for (auto& t : m.tri)
        for (int& v : t) {
            if (used[v] < 0) {
                used[v] = static_cast<int>(out.pos.size());
                out.pos.push_back(m.pos[v]);
            }
            v = used[v];
        }
    out.tri = m.tri;
    out.uv = m.uv;
    out.computeNormals();
    return out;
}

inline PatternParams expandedD1() {
    PatternParams pp;  // the expanded preset: 8 x 6.5 mm, D 1.0
    return pp;
}

inline MeshKnitParams tiledParams(const PatternParams& pp) {
    MeshKnitParams p;
    p.tile.pp = pp;
    p.tile.g.w = pp.waleSpacing;
    p.tile.g.h = pp.courseSpacing;
    return p;
}

} // namespace knit
