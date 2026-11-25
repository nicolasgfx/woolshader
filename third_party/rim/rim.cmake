# Stitch Meshing v2, stage 1: build patch for RIM (Gao, Jakob, Tarini, Panozzo 2017, "Robust Hex-Dominant Mesh
# Generation using Field-Guided Polyhedral Agglomeration"; gaoxifeng/robust_hex_dominant_meshing at 83dd574,
# MIT), run as the fetch's PATCH_COMMAND in the source directory (cmake/dependencies.cmake, KNIT_RIM_PATCH 5). Surface path only:
#   - TetGen removed (the include and the body of the volume path's tet_meshing);
#   - an Eigen 3.4 type fix (Float is double: the per-face centres are Vector4f = Vector4d, not Eigen's float);
#   - system("PAUSE") on three inconsistency paths becomes an exception (headless: the remesh fails loudly);
#   - the seed: build()'s random start (the frames and lattice origins) reads knit_rim_seed (third_party/rim/rim_remesh.cpp);
#     0 keeps pcg32's default state, so seed 0 is RIM's own start;
#   - the 2-RoSy boundary (2026-11-14): boundary vertices smoothed and snapped along or across the boundary (below);
#   - two hook calls into our code (third_party/rim/knit_rim_hooks.h, written from Wu et al. 2018 section 4, not from
#     any other implementation): the 2-RoSy orientation compatibility in smoothOrientationsTri (the
#     original 4-RoSy functions unless knit_rim_rosy == 2), and the post-processing inside
#     meshExtraction2D's collapse / split loop (the step after the extraction is our own normalizePolygons,
#     src/stitchmesh/polygon_mesh.cpp, not a hook);
#   - the rectangular lattice (aspect, below; inactive at 1);
#   - three crash fixes: simple_polygon (an empty boundary; a discarded vertex on the boundary) and remove_doublets2D
#     (a shared-edge chain without two end vertices), all below.
# No other change. Line endings are normalized to LF first (the sources are CRLF). Every replacement must
# apply exactly as written, else the configure step fails.
cmake_minimum_required(VERSION 3.20)
set(SRC "${CMAKE_CURRENT_SOURCE_DIR}/src")

# One replacement: `count` occurrences of `old` (exactly) become `new`. Arguments, not lists: the strings
# contain semicolons.
function(knit_replace file count old new)
    file(READ "${SRC}/${file}" S)
    string(REPLACE "\r\n" "\n" S "${S}")
    set(rest "${S}")
    set(found 0)
    string(LENGTH "${old}" len)
    while(TRUE)
        string(FIND "${rest}" "${old}" at)
        if(at EQUAL -1)
            break()
        endif()
        math(EXPR found "${found} + 1")
        math(EXPR at "${at} + ${len}")
        string(SUBSTRING "${rest}" ${at} -1 rest)
    endwhile()
    if(NOT found EQUAL count)
        message(FATAL_ERROR "rim patch: ${file}: expected ${count} of '${old}', found ${found}")
    endif()
    string(REPLACE "${old}" "${new}" S "${S}")
    file(WRITE "${SRC}/${file}" "${S}")
endfunction()

# Crash fix (2026-10-29, the stitch mesh of the Stanford bunny): simple_polygon reads pes[0], the merged
# region's boundary, also when every edge of the region is interior (no boundary: pes empty, a read at
# address 0). An empty boundary is not a simple polygon: false, which every caller already handles (that
# collapse is skipped). Its own probe, so a source tree patched before this fix gets it too.
file(READ "${SRC}/extract.cpp" probe)
string(FIND "${probe}" "knit: no boundary edge" fixed)
if(fixed EQUAL -1)
    knit_replace(extract.cpp 1 "\t//test nvs for each v\n"
        "\tif (pes.empty()) return false;  // knit: no boundary edge (every edge interior): not a polygon\n\t//test nvs for each v\n")
endif()

# Crash fix (2026-11-13, the second intermittent access violation, in tagging_collapseTri; the same reproducer file):
# simple_polygon maps the boundary edges' vertices through V_map, where a discarded vertex is INVALID_V = one past the
# last vertex (meshExtraction2D, = V_map.size()); such an id then indexes Reverse_V_nts one past its end in the collapse
# loop (extract.cpp, std::replace on Reverse_V_nts[vid]), which reads whatever follows on the heap: a fault or silent
# corruption (measured on the car, seed 1, extrinsic attempt: vid 950 of 950, in 5 of 6 runs a crash). A polygon whose
# boundary touches a discarded vertex is not a simple polygon: false, as the other failure paths of the function.
# (mV_flag is not the bound: it holds 3 x the vertex count, Eigen's size() of a 3 x N matrix.) Its own probe, as above.
file(READ "${SRC}/extract.cpp" probe)
string(FIND "${probe}" "knit: a discarded vertex" fixed)
if(fixed EQUAL -1)
    knit_replace(extract.cpp 1 "\t//test nvs for each v\n"
        "\tfor (uint32_t i = 0; i < pes.size(); ++i) {  // knit: a discarded vertex (V_map one past the end) on the boundary: not a polygon\n\t\tconst uint32_t a = V_map[std::get<0>(mEs[pes[i]])], b = V_map[std::get<1>(mEs[pes[i]])];\n\t\tif (a >= V_map.size() || b >= V_map.size()) return false;\n\t}\n\t//test nvs for each v\n")
endif()

# Crash fix (2026-11-13, intermittent access violation in the extraction, reproduced on a car ring with seed 1):
# remove_doublets2D collects the two end vertices of the edge chain two faces
# share and refuses more than two, but not fewer; when the shared edges close a loop the list is empty and vs_vd[1]
# reads stale heap memory, which indexes Vs_nes far out of bounds. Fewer than two ends is skipped like more than two
# (the doublet stays; the extraction's result is validated by the caller). Its own probe, as above.
file(READ "${SRC}/extract.cpp" probe)
string(FIND "${probe}" "vs_vd.size() != 2" fixed)
if(fixed EQUAL -1)
    knit_replace(extract.cpp 1 "if (vs_vd.size() > 2) { continue;" "if (vs_vd.size() != 2) { continue;")
endif()

# The boundary of a 2-RoSy field (2026-11-14, Stitch Meshing recovery M2): RIM fixes every boundary vertex to its
# boundary edge's direction, which is RIM's 4-RoSy constraint "aligned with the boundary" (along it and across it are the
# same 4-RoSy direction). Kept as a fixed 2-RoSy direction it says "along" on every side, so a rectangle's adjacent sides
# disagree by 90 degrees and the field needs a singularity at each corner (the solved flat: 4 orientation singular
# points, 42 triangles; m2/flat-rosy2.png). With 2-RoSy a boundary vertex is smoothed like an interior one and then
# snapped to the nearer of its boundary's two directions (knitBoundarySnap): the same constraint, written for 2-RoSy.
# Its own probe, as above.
file(READ "${SRC}/orientations.cpp" probe)
string(FIND "${probe}" "knitBoundarySnap" fixed)
if(fixed EQUAL -1)
    knit_replace(orientations.cpp 1
        "\t\t\t\tif (nV_boundary_flag[l][i]) {\n\t\t\t\t\tQ_new.col(i) = Q.col(i);\n\t\t\t\t\tcontinue;\n\t\t\t\t}\n"
        "\t\t\t\tconst bool knitBoundary2 = nV_boundary_flag[l][i] && knit_rim_rosy == 2;  // knit: smoothed, then snapped\n\t\t\t\tif (nV_boundary_flag[l][i] && !knitBoundary2) {\n\t\t\t\t\tQ_new.col(i) = Q.col(i);\n\t\t\t\t\tcontinue;\n\t\t\t\t}\n")
    knit_replace(orientations.cpp 1
        "                if (q_i != Vector3f::Zero())\n                    Q_new.col(i) = (q_i - n_i.dot(q_i) * n_i).normalized();\n"
        "                if (knitBoundary2)\n                    Q_new.col(i) = knitBoundarySnap(Q.col(i), n_i, q_i);  // knit: along or across the boundary\n                else if (q_i != Vector3f::Zero())\n                    Q_new.col(i) = (q_i - n_i.dot(q_i) * n_i).normalized();\n")
endif()

file(READ "${SRC}/hierarchy.h" probe)
string(FIND "${probe}" "knit build" done)
if(NOT done EQUAL -1)
    return()
endif()

knit_replace(hierarchy.h 1 "#include \"tetgen.h\""
    "// tetgen.h removed (knit build: surface path only)\nextern uint64_t knit_rim_seed;  // knit: third_party/rim/rim_remesh.cpp")

# The volume path's tet_meshing (TetGen) is replaced by a stub.
file(READ "${SRC}/hierarchy.cpp" S)
string(REPLACE "\r\n" "\n" S "${S}")
string(FIND "${S}" "bool MultiResolutionHierarchy::tet_meshing()" a)
string(FIND "${S}" "void MultiResolutionHierarchy::build()" b)
if(a EQUAL -1 OR b EQUAL -1 OR NOT a LESS b)
    message(FATAL_ERROR "rim patch: hierarchy.cpp: tet_meshing / build not found in order")
endif()
string(SUBSTRING "${S}" 0 ${a} head)
string(SUBSTRING "${S}" ${b} -1 tail)
file(WRITE "${SRC}/hierarchy.cpp"
     "${head}bool MultiResolutionHierarchy::tet_meshing()\n{\n\treturn false;  // knit build: no TetGen (volume path unused)\n}\n${tail}")

# The seed: the one pcg32 object of build() (the start frames and lattice origins). The sweeps' shuffles are seeded
# by RIM's iteration counters and the row index, unchanged. 0 keeps pcg32's default state and stream.
knit_replace(hierarchy.cpp 1 "\tpcg32 rng;\n"
    "\tpcg32 rng(knit_rim_seed ? knit_rim_seed : PCG32_DEFAULT_STATE, PCG32_DEFAULT_STREAM);  // knit: the request's seed\n")

knit_replace(orientations.cpp 1 "#include \"positions.h\"\n" "#include \"positions.h\"\n#include \"knit_rim_hooks.h\"\n")
knit_replace(orientations.cpp 1
    "q_j = applyRotationExtrinsic((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);"
    "q_j = knitRotationExtrinsic((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);")
knit_replace(orientations.cpp 1
    "q_j = applyRotation((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);"
    "q_j = knitRotation((q_i == Vector3f::Zero()) ? Q.col(i) : q_i, n_i, q_j, n_j);")

knit_replace(extract.cpp 1 "#include \"positions.h\"\n" "#include \"positions.h\"\n#include \"knit_rim_hooks.h\"\n")
knit_replace(extract.cpp 1 "Eigen::Vector4f::Zero()" "Vector4f::Zero()")
knit_replace(extract.cpp 3 "system(\"PAUSE\");" "throw std::runtime_error(\"rim extraction: inconsistent edge data\");")
knit_replace(extract.cpp 1
    "\t\tif (splitting) {\n\t\t\twhile (split_face2D(true));\n\t\t}\n"
    "\t\tif (splitting) {\n\t\t\twhile (split_face2D(true));\n\t\t}\n\t\tknit_wu_loop(*this);  // knit: Wu 2018 section 4, post-processing interleaved\n")
# The rectangular lattice (M5; inactive while the aspect is 1), in RIM's 2D position functions. knit_rim_aspect (rim_driver.cpp) = t spacing / q spacing, t = n x q: the t coordinate is divided by
# it before rounding and multiplied back after (aspect 1: exact, the original arithmetic). While it is not 1 the
# neighbour matching is 2-way (q never swapped with t: the spacings differ): applyRotation, applyRotationKeep and
# findRotation return the neighbour's q (rotated into the plane) up to sign.
knit_replace(orientations.h 1 "#include \"quat.h\"\n" "#include \"quat.h\"\n\nextern Float knit_rim_aspect;  // knit: third_party/rim/rim.cmake\n")
knit_replace(orientations.h 1
    "inline int findRotation(const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) {\n"
    "inline int findRotation(const Vector3f &q0, const Vector3f &n0, const Vector3f &_q1, const Vector3f &n1) {\n    if (knit_rim_aspect != 1) return rotateVectorIntoPlane(_q1, n1, n0).dot(q0) > 0 ? 0 : 2;  // knit: 2-way\n")
knit_replace(orientations.h 1
    "inline Vector3f applyRotation(const Vector3f &q0, const Vector3f &n0,\n                              const Vector3f &_q1, const Vector3f &n1) {\n"
    "inline Vector3f applyRotation(const Vector3f &q0, const Vector3f &n0,\n                              const Vector3f &_q1, const Vector3f &n1) {\n    if (knit_rim_aspect != 1) { const Vector3f r = rotateVectorIntoPlane(_q1, n1, n0); return r * std::copysign((Float) 1, r.dot(q0)); }  // knit: 2-way\n")
knit_replace(orientations.h 1
    "inline Vector3f applyRotationKeep(const Vector3f &q0, const Vector3f &n0,\n                                  const Vector3f &_q1, const Vector3f &n1) {\n"
    "inline Vector3f applyRotationKeep(const Vector3f &q0, const Vector3f &n0,\n                                  const Vector3f &_q1, const Vector3f &n1) {\n    if (knit_rim_aspect != 1) return _q1 * std::copysign((Float) 1, rotateVectorIntoPlane(_q1, n1, n0).dot(q0));  // knit: 2-way\n")
knit_replace(positions.h 1
    "    Vector2f rel = round(Vector2f(M.transpose() * ((o0 - o1) * invScale)));\n    return std::make_pair(o1 + (M * rel) * scale, rel.cast<int>());"
    "    Vector2f rel = Vector2f(M.transpose() * ((o0 - o1) * invScale));\n    rel[1] /= knit_rim_aspect;  // knit: rectangular lattice\n    rel = round(rel);\n    const Vector2i ri = rel.cast<int>();\n    rel[1] *= knit_rim_aspect;\n    return std::make_pair(o1 + (M * rel) * scale, ri);")
knit_replace(positions.h 2
    "\tVector2f vv = Vector2f(M.transpose() * ((o0 - o1) * invScale));\n"
    "\tVector2f vv = Vector2f(M.transpose() * ((o0 - o1) * invScale));\n\tvv[1] /= knit_rim_aspect;  // knit: rectangular lattice\n")
knit_replace(positions.h 1
    "    Vector2f rel = M.transpose() * ((ref - o) * invScale);\n    return o + (M * round(rel)) * scale;"
    "    Vector2f rel = M.transpose() * ((ref - o) * invScale);\n    rel[1] /= knit_rim_aspect;  // knit: rectangular lattice\n    rel = round(rel);\n    rel[1] *= knit_rim_aspect;\n    return o + (M * rel) * scale;")
message(STATUS "rim patch applied to ${SRC}")
