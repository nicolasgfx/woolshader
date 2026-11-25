#pragma once
// Method A's row-time field: the geodesic
// distance on a triangle mesh from a set of source vertices (a cast-on boundary), by fast marching (Kimmel and
// Sethian 1998). The rings of a component are this field's level sets at the row spacing (row_rings.h).
#include <string>
#include <vector>

#include "model/mesh.h"

namespace knit {

struct TimeFieldReport {
    int accepted = 0;            // vertices that received a distance
    int triangleUpdates = 0;     // accepted values that came from a planar wavefront through a triangle
    int edgeUpdates = 0;         // ... from an edge (no causal triangle update: obtuse corner or a one-sided front)
    double ms = 0;
};

// t[v] = 0 on the sources. A vertex's value comes from the planar wavefront through two accepted corners of a triangle
// when its characteristic through the vertex crosses their edge (the causal case), otherwise from the cheaper of the
// two edges. Unreachable vertices keep +infinity. False only when there is no source or a source index is invalid.
// seedValues (optional, one per source): the sources' starting values instead of 0, for a source point inside a triangle
// (its corners seeded with their distances from the point).
bool geodesicDistance(const TriMesh& m, const std::vector<int>& sources, std::vector<double>& t, TimeFieldReport& report,
                      std::string* err, const std::vector<double>* seedValues = nullptr);

} // namespace knit
