#pragma once
// The hand-authored coarse inputs of stage 4: small labelled, directed coarse
// meshes whose typed stitch mesh, ports and strands are asserted by the tests. They
// are library code so the lab can write and draw them too; nothing here is extracted from a surface.
#include <string>
#include <vector>

#include "model/mesh.h"
#include "stitchmesh/polygon_mesh.h"

namespace knit {

enum class StitchFixture {
    ClosedRows,            // 1: a cylinder of 3 around, 2 rows (9 / 15 / 6): four closed stitch rows
    SpecialQuads,          // 2: 4 around, the middle row short with two row-end triangles: two special quads
    Increase,              // 3: the middle row grows from 3 to 4 through a triangle with two wale edges: two kp pentagons
    Decrease,              // 4: the mirror: two d12k pentagons
    MismatchBetweenRows,   // 5a: fixture 1 with the second row growing down (a ring of both-top edges: joining)
    MismatchSplitting,     // 5c: fixture 1 with the FIRST row growing down (a ring of both-bottom edges: splitting)
    MismatchFlatJoining,   // 5d: two rows of three rectangles on a sheet, the upper growing down (an open chain of both-top edges)
    MismatchFlatSplitting, // 5e: the same with the lower row growing down (an open chain of both-bottom edges)
    MismatchWithinRow,     // 5b: a hairpin on a flat sheet, the fold edges top for both arms
    IncreaseFlat,          // 3 flat: rows of rectangles on a sheet, the middle row one stitch wider by a full-width
                           //         triangle in the middle: the kp pair with no shear anywhere (a diagnostic, open rows)
    DecreaseFlat,          // 4 flat: the mirror, one stitch narrower
    // From the frozen sphere (recovery M7, stitchmesh_lab --cutout-vertex 24 --cutout-rings 1 on m4/pipeline-sphere): the four
    // coarse faces round a pinch vertex where a both-bottom edge's two receivers, one through a special quad, stand side by
    // side in one row (a cast-on pulled by an increase), next to a facing split cast-on; 5.5 mm gauge, its own facets as surface
    CastOnSideBySide,
    // From the frozen sphere (--cutout-vertex 28 --cutout-rings 1): five coarse faces round a vertex where two both-bottom
    // edges face a plain knit on one side and a kp pentagon on the other (the M5 contract: the kp pair is drawn through its
    // half of the split cast-on's wrap); a third both-bottom edge's receiver lies beyond the cut (unlinked, by the cut)
    CastOnKpReceivers,
    // From the frozen sphere (--cutout-vertex 11 --cutout-rings 1): the pole where three special quads pass one cast-on's
    // loop to two kp pentagons standing side by side in ONE row: the loop pulled by four neighbouring stitches (k, p, k, p)
    CastOnPoleKp
};

struct CoarseFixture {
    std::string name;
    PolygonMesh mesh;
    TriMesh surface;   // the smooth reference surface the coarse mesh approximates (the cylinder of radius 10, or the
                       // sheet's plane), for the relaxation and the yarn's feet; the coarse mesh's own facets are not it
    std::vector<std::vector<char>> top;   // DirectedMesh::top: per face, per corner edge: 1 top, 0 bottom, -1 wale
};

// The document's fixtures (around = 0: 3 or 4 coarse faces round, 10 mm coarse rows, radius 10), or the same patterns
// with `around` coarse faces round the cylinder and 5 mm coarse rows: the fine versions the yarn is placed on, where a
// stitch bends by 360 / (2 around) degrees instead of 60 or 45 (the document's counts do not hold for them).
// `radius`: the ring fixtures on a cylinder of that radius instead of 10 (x and y scaled, the topology and the rows' height
// unchanged): with around = 48 and radius 40 the stitches are as wide as around = 12 at radius 10, and a connector row
// extruded inwards (two stitch heights and a bit, StitchYarnOptions) shrinks by 13 % instead of by half.
CoarseFixture stitchFixture(StitchFixture which, int around = 0, double radius = 10.0);
const char* stitchFixtureName(StitchFixture which);   // closed-rows, special-quads, kp, d12k, mismatch-between-rows, mismatch-within-row, kp-flat, d12k-flat
bool parseStitchFixture(const std::string& name, StitchFixture& out);

} // namespace knit
