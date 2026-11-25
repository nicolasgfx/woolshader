// A built stitch mesh's display gate and its findings in words (stitch_mesh_data.h).
#include "model/stitch_mesh_data.h"

#include <cmath>
#include <cstdio>

namespace knit {

namespace {

const char* const kSeparator = " \xC2\xB7 ";   // a middle dot between the headline's parts

std::string countOf(int n, const char* one, const char* many) { return std::to_string(n) + " " + (n == 1 ? one : many); }

// A finding's short form: the line up to its parenthesis.
std::string shortForm(const std::string& line) {
    const size_t paren = line.find(" (");
    return paren == std::string::npos ? line : line.substr(0, paren);
}

bool finite(const Vec3& p) { return std::isfinite(p.x) && std::isfinite(p.y) && std::isfinite(p.z); }

} // namespace

StitchMeshUse stitchMeshUse(const StitchMeshData& m) {
    StitchMeshUse use;
    const int faces = m.faces();
    if (faces == 0) {
        use.reason = "the build made no stitch faces";
        return use;
    }
    if (m.faceStart.front() != 0 || m.faceStart.back() != static_cast<int>(m.corner.size())) {
        use.reason = "the face table does not match its corner list";
        return use;
    }
    const int vertices = static_cast<int>(m.pos.size());
    for (int v = 0; v < vertices; ++v)
        if (!finite(m.pos[v])) {
            use.reason = "vertex " + std::to_string(v) + " has a non-finite coordinate";
            return use;
        }
    for (int f = 0; f < faces; ++f) {
        if (m.sides(f) < 3) {
            use.reason = "face " + std::to_string(f) + " has " + std::to_string(m.sides(f)) + " corners";
            return use;
        }
        for (int k = m.faceStart[f]; k < m.faceStart[f + 1]; ++k)
            if (m.corner[k] < 0 || m.corner[k] >= vertices) {
                use.reason = "face " + std::to_string(f) + " names vertex " + std::to_string(m.corner[k]) + " of " + std::to_string(vertices);
                return use;
            }
    }
    use.drawable = true;
    return use;
}

std::vector<std::string> stitchMeshFindings(const StitchMeshQuality& q) {
    std::vector<std::string> lines;
    char b[240];
    // structure: faces that are wrong as polygons, whatever their fit
    if (q.inverted > 0) lines.push_back(countOf(q.inverted, "face is", "faces are") + " wound against the piece (inverted)");
    if (q.degenerate > 0) lines.push_back(countOf(q.degenerate, "face has", "faces have") + " almost no area (degenerate)");
    if (q.untriangulable > 0) lines.push_back(countOf(q.untriangulable, "face cannot", "faces cannot") + " be triangulated in its own plane");
    if (q.intersectingPairs > 0) lines.push_back(countOf(q.intersectingPairs, "pair", "pairs") + " of faces intersect");
    if (q.unsupported > 0) lines.push_back(countOf(q.unsupported, "cell has", "cells have") + " no knitting construction yet");
    if (!q.connectivityPass) lines.push_back("Connectivity does not pass (" + (q.connectivityDetail.empty() ? std::string("unresolved") : q.connectivityDetail) + ")");
    // fit to the piece
    if (q.facesOverFit > 0) {
        std::snprintf(b, sizeof b, "%s the %.2f mm fit target (worst %.2f mm)", countOf(q.facesOverFit, "face exceeds", "faces exceed").c_str(), q.fitTargetMm, q.maxDeviationMm);
        lines.push_back(b);
    }
    // coverage and boundary
    if (!q.uncoveredMeasured) lines.push_back("Coverage cannot be measured (the whole piece lies within half a row of its border)");
    else if (q.uncoveredFraction > q.uncoveredLimit) {
        std::snprintf(b, sizeof b, "%.1f %% of the piece is not covered (limit %.1f %%)", 100 * q.uncoveredFraction, 100 * q.uncoveredLimit);
        lines.push_back(b);
    }
    if (q.coveredTwice > 0) lines.push_back(countOf(q.coveredTwice, "place is", "places are") + " covered twice (a fold or an overlap)");
    if (q.areaRatio < q.areaLow || q.areaRatio > q.areaHigh) {
        std::snprintf(b, sizeof b, "The mesh area is %.3f of the piece's (expected %.2f to %.2f)", q.areaRatio, q.areaLow, q.areaHigh);
        lines.push_back(b);
    }
    if (q.boundaryMm > q.boundaryLimitMm) {
        std::snprintf(b, sizeof b, "The border lies up to %.2f mm from the piece's edge (limit %.2f mm)", q.boundaryMm, q.boundaryLimitMm);
        lines.push_back(b);
    }
    return lines;
}

std::string stitchMeshRelaxationNote(const StitchMeshQuality& q) {
    if (!q.relaxed) return "";
    if (!q.relaxRan) return "The relaxation did not run (" + q.relaxNote + ")";
    if (q.relaxStop != "converged") return "The relaxation stopped without converging (" + q.relaxStop + ")";
    return "";
}

std::string stitchMeshHeadline(const StitchMeshUse& use, const StitchMeshQuality& q) {
    if (!use.drawable) return "No stitch mesh: " + use.reason;
    std::string line = "Mesh generated" + std::string(kSeparator) + countOf(q.faces, "face", "faces") + kSeparator + countOf(q.rows, "row", "rows") + kSeparator;
    const std::vector<std::string> findings = stitchMeshFindings(q);
    if (findings.empty()) return line + (q.clean ? "passes its checks" : "does not pass its checks (see the findings)");
    line += shortForm(findings.front());
    if (findings.size() > 1) line += kSeparator + countOf(static_cast<int>(findings.size()) - 1, "more finding", "more findings");
    return line;
}

void stitchMeshBounds(const StitchMeshData& m, Vec3& lo, Vec3& hi) {
    lo = hi = Vec3{};
    for (size_t v = 0; v < m.pos.size(); ++v) {
        lo = v == 0 ? m.pos[v] : componentMin(lo, m.pos[v]);
        hi = v == 0 ? m.pos[v] : componentMax(hi, m.pos[v]);
    }
}

} // namespace knit
