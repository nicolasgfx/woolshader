#pragma once
// Our code called from the patched RIM sources (third_party/rim/rim.cmake), written from the text of Wu, Gao,
// Ferguson, Panozzo, Yuksel 2018 ("Stitch Meshing", section 4). Included INSIDE namespace rim (the wrapped RIM
// translation units, rim_remesh.cpp), after positions.h.

// The orientation symmetry of smoothOrientationsTri: 4 (RIM's original functions) or 2 (Wu 2018 section 4: "a 2-RoSy
// symmetry instead of a 4-RoSy symmetry"). Set by rim_remesh.cpp per remesh (one remesh at a time).
extern int knit_rim_rosy;

// 2-RoSy compatibility: the neighbour's direction (rotated into this tangent plane, intrinsic) or itself
// (extrinsic), times the sign of its dot product with ours: q1 * sign(q0 . q1). 4-RoSy: RIM's applyRotation.
inline Vector3f knitRotation(const Vector3f& q0, const Vector3f& n0, const Vector3f& q1, const Vector3f& n1) {
    if (knit_rim_rosy != 2) return applyRotation(q0, n0, q1, n1);
    const Vector3f r = rotateVectorIntoPlane(q1, n1, n0);
    return r * std::copysign((Float) 1, q0.dot(r));
}
// A boundary vertex of a 2-RoSy field: of its boundary's two directions (t, its fixed boundary direction, and n x t,
// across), the one nearer the smoothed average q of its neighbours, with q's sign; t itself when q is zero. RIM's
// 4-RoSy boundary constraint ("aligned with the boundary") written for 2-RoSy (rim.cmake, the boundary of a 2-RoSy field).
inline Vector3f knitBoundarySnap(const Vector3f& t, const Vector3f& n, const Vector3f& q) {
    if (q == Vector3f::Zero()) return t;
    const Vector3f across = n.cross(t);
    const Float a = q.dot(t), b = q.dot(across);
    return std::abs(a) >= std::abs(b) ? t * std::copysign((Float) 1, a) : across * std::copysign((Float) 1, b);
}

inline Vector3f knitRotationExtrinsic(const Vector3f& q0, const Vector3f& n0, const Vector3f& q1, const Vector3f& n1) {
    if (knit_rim_rosy != 2) return applyRotationExtrinsic(q0, n0, q1, n1);
    return q1 * std::copysign((Float) 1, q0.dot(q1));
}

// Wu 2018 section 4 post-processing inside the extraction (rim_remesh.cpp; a no-op unless enabled for the remesh):
// knit_wu_loop, every round of meshExtraction2D's collapse / split loop, before its collapse queue is
// rebuilt: (2) polygons with more than 5 sides are split along the diagonal most aligned with the field,
// recursively; (1) in each triangle the edge whose opposite angle is closest to 90 degrees (not a sharp
// feature: the two faces' normals have a negative dot product) is marked a diagonal (RIM's D tag: the
// collapse merges its two faces when that gives a simple polygon). The loop's own convergence rule ("no
// change of the face count") ends the interleaving. The step after the extraction (large polygons split once
// more, pentagons into three triangles) is normalizePolygons in src/stitchmesh/polygon_mesh.cpp.
class MultiResolutionHierarchy;
void knit_wu_loop(MultiResolutionHierarchy& h);
