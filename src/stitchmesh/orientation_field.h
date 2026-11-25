#pragma once
// Method B's orientation field on the reference surface
// (Wu et al. 2018 section 4: RIM with a 2-RoSy symmetry): solved by RIM's smoothing, or fixed by an axis guide,
// without RIM's position field or extraction (checkpoint 1 analyses the field; no coarse mesh is made here).
#include <cstdint>
#include <string>
#include <vector>

#include "model/mesh.h"

namespace knit {

struct OrientationFieldRequest {
    double spacing = 11;          // RIM's target edge [mm]: the coarse edge, 2 x the stitch width
    int rosy = 2;                 // 2 (Wu 2018 section 4) or 4 (RIM's own)
    bool extrinsic = false;       // Wu 2018's default: the intrinsic field
    int sweepsPerLevel = 20;      // dbf5f41's setting (RIM itself runs 100)
    uint64_t seed = 1;            // RIM's random start (0: its own default)
    std::vector<Vec3> guide;      // per vertex, a direction that fixes the field (empty: solved)
};

struct OrientationField {
    std::vector<Vec3> q;              // per input vertex: the field's representative direction, unit, tangent
    std::vector<Vec3> singularities;  // RIM's orientation singularities, as points on the surface
    int finestVertices = 0;           // RIM's finest level; equals the input's vertex count when RIM subdivided nothing
    std::vector<int> levelSizes;
    double ms = 0;
    std::string log;                  // RIM's console output
};

// False when RIM refuses the input or subdivides it (then q would not be per input vertex).
bool computeOrientationField(const TriMesh& m, const OrientationFieldRequest& r, OrientationField& out, std::string* err);

// The guide "rows around the axis": at every vertex the direction around the axis line through the origin, in the
// vertex's tangent plane (dbf5f41 remesh_stage `globalAxisGuide("ring-z")` for axis z).
std::vector<Vec3> ringAxisGuide(const TriMesh& m, const Vec3& axis);

// The angle [degrees] between two tangent directions taken modulo `symmetryDegrees` (180 for a line field, 90 for
// the lattice that a 2-RoSy field and its normal-rotated partner span): in [0, symmetryDegrees / 2].
double directionDeviation(const Vec3& a, const Vec3& b, const Vec3& normal, double symmetryDegrees);

} // namespace knit
