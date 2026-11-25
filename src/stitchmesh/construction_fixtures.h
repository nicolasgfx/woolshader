#pragma once
// The minimal construction fixtures of checkpoint 2 and their negative controls (a
// construction is admitted on a model only after its minimal hand-built neighbourhood passes the independent
// check and its negative controls fail). Library code so the lab can write and draw the same cases the tests assert.
//
// Ring tubes: ring k has counts[k] vertices on a cylinder at height k * rowHeight, vertex 0 on the cut (phi = 0); strip k's
// faces come from assembleStrip with the given shaping slots and the lookup; ring 0 is the cast-on, the last ring the
// bind-off; the shift-path is the chain of vertex-0 wale edges. The special-quad neighbourhood is the recovery branch's
// hand-authored SpecialQuads coarse fixture (dbf5f41 stitch_fixtures) subdivided by the ported buildStitchMesh.
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"

namespace knit {

struct FixtureCase {
    std::string name;
    std::string entry;          // the catalog entry it admits or attacks
    std::string purpose;        // what it shows, or the defect a negative control carries
    bool negative = false;
    TypedStitchMesh mesh;
    TriMesh reference;          // the surface the check measures geometry against
    std::vector<DeclaredBoundary> declared;
    bool runShift = false;
    ShiftReading reading = ShiftReading::PathSideDown;
    Status expectConnectivity = Status::Pass;
    Status expectShift = Status::Unresolved;   // Unresolved: not run
};

// One ring tube. shaping[k]: the shaping slots of strip k (sizes must match the count changes).
bool ringTubeFixture(const std::vector<int>& counts, const std::vector<std::vector<int>>& shaping, double radius, double rowHeight,
                     TypedStitchMesh& mesh, TriMesh& reference, std::vector<DeclaredBoundary>& declared, std::string* err, bool boundaryShaping = false);

// A cylinder of radius r about +z from z0 to z1 (192 segments, rings `step` apart), counter-clockwise from outside: the
// reference of the ring-tube fixtures (relaxation uses one that ends exactly at the fixture's border rings).
TriMesh cylinderReference(double r, double z0, double z1, double step);

// A surface of revolution about +z through rings of radius rho[k] at height z[k], piecewise linear between them, `around`
// segments, `sub` bands per row; it ends exactly at the first and last ring (the shaping and yarn fixtures' reference).
TriMesh revolutionThroughRings(const std::vector<double>& rho, const std::vector<double>& z, int around, int sub);

// The boundary loops of a triangle mesh, each as vertex ids in walking order (the triangle lies on the left).
std::vector<std::vector<int>> triangleBoundaryLoops(const TriMesh& m);

// Every fixture and negative control, in report order.
std::vector<FixtureCase> constructionFixtureCases();

// The imported arm's constructions, kept apart from the CP2 set:
// A pole disc: the pole at the origin of the plane z = 0, ring k (k = 1 .. K) a circle of radius k rowHeight with counts[k - 1]
// vertices, clockwise seen from +z (so the faces are counter-clockwise from outside, +z), vertex 0 on the +x axis. Strip 0 is
// the gathered pole start (counts[0] triangles (pole, ring 1 j + 1, ring 1 j), roles W, T, W); strip k >= 1 comes from
// assembleStrip with shaping[k - 1]; the last ring is the bind-off; the shift-path runs from ring 1 to ring K.
bool poleDiscFixture(const std::vector<int>& counts, const std::vector<std::vector<int>>& shaping, double rowHeight, TypedStitchMesh& mesh, TriMesh& reference,
                     std::vector<DeclaredBoundary>& declared, std::string* err);
// The pole-start, adjacent-shaping and boundary-shaping fixtures and their negative controls.
std::vector<FixtureCase> armConstructionFixtureCases();

} // namespace knit
