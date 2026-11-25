// The yarn-generation commands of stitchmesh_lab:
//
//   yarn-fixtures   minimal neighbourhoods (regular patch with cast-on and bind-off, isolated increase and decrease, the
//                   boundary-shaping variants, the gathered pole start, dense shaping) relaxed with the CP3 model on a
//                   surface that fits their counts, then yarn embedded and checked; plus negative controls of the checks
//   yarn-embed      one frozen relaxed stitch mesh (a marked run's after/ files and reference.obj, read back, not rebuilt):
//                   yarn embedded and checked, centrelines, controls and tubes written
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "stitchmesh_lab_shared.h"
#include "stitchmesh/construction_fixtures.h"
#include "model/mesh_query.h"
#include "stitchmesh/coverage.h"
#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"
#include "stitchmesh/yarn_embed.h"

namespace fs = std::filesystem;
using namespace knit;

namespace {

using lab::YarnFixture;

constexpr double kPi = 3.14159265358979323846;

struct Opts {
    std::map<std::string, std::string> kv;
    std::string get(const std::string& k, const std::string& d = "") const {
        const auto it = kv.find(k);
        return it == kv.end() ? d : it->second;
    }
    double number(const std::string& k, double d) const {
        const auto it = kv.find(k);
        return it == kv.end() ? d : std::atof(it->second.c_str());
    }
};

Opts parse(int argc, char** argv) {
    Opts o;
    for (int i = 2; i + 1 < argc; i += 2) o.kv[argv[i]] = argv[i + 1];
    return o;
}

void putf(std::ostream& o, const char* fmt, ...) {
    char b[2048];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(b, sizeof b, fmt, ap);
    va_end(ap);
    o << b;
}

void writeTriObj(const std::string& path, const TriMesh& m) {
    std::ofstream f(path, std::ios::binary);
    for (const Vec3& q : m.pos) putf(f, "v %.6f %.6f %.6f\n", q.x, q.y, q.z);
    for (const auto& t : m.tri) putf(f, "f %d %d %d\n", t[0] + 1, t[1] + 1, t[2] + 1);
}

bool readObjTriangles(const std::string& path, TriMesh& ref) {
    ref = TriMesh{};
    std::ifstream f(path);
    std::string line;
    while (std::getline(f, line)) {
        if (line.rfind("v ", 0) == 0) {
            Vec3 p;
            std::sscanf(line.c_str() + 2, "%lf %lf %lf", &p.x, &p.y, &p.z);
            ref.pos.push_back(p);
        } else if (line.rfind("f ", 0) == 0) {
            std::array<int, 3> t{};
            std::sscanf(line.c_str() + 2, "%d %d %d", &t[0], &t[1], &t[2]);
            ref.tri.push_back({t[0] - 1, t[1] - 1, t[2] - 1});
        }
    }
    if (ref.tri.empty()) return false;
    ref.computeNormals();
    return true;
}

YarnSettings settingsFrom(const Opts& a) {
    YarnSettings st;
    st.diameter = a.number("--diameter", st.diameter);
    st.joinPoleRound = a.get("--join-pole", "1") == "1";
    st.bindingThread = a.get("--binding", "1") == "1";
    st.spanTangents = a.get("--tangents", "catmull") == "span";
    st.templateSet = a.get("--templates", "baseline");
    st.placement = a.get("--placement", toString(YarnSettings{}.placement)) == "face" ? KnotPlacement::FaceInterpolation : KnotPlacement::ReferenceProjection;
    return st;
}

// Everything a yarn run writes: centrelines, controls, tubes, the check and its per-face table.
void writeYarn(const std::string& dir, const TypedStitchMesh& s, const YarnEmbedding& e, const YarnCheck& c, const YarnSettings& st, const std::string& title) {
    fs::create_directories(dir);
    {
        std::ofstream f(dir + "/yarn_centrelines.obj", std::ios::binary);
        f << "# " << title << ": yarn centrelines [mm], one object per strand\n";
        int base = 1;
        for (size_t k = 0; k < e.strands.size(); ++k) {
            const YarnStrand& sd = e.strands[k];
            putf(f, "o strand_%zu\n", k);
            for (const Vec3& p : sd.samples) putf(f, "v %.6f %.6f %.6f\n", p.x, p.y, p.z);
            f << "l";
            for (size_t i = 0; i < sd.samples.size(); ++i) putf(f, " %d", base + static_cast<int>(i));
            if (sd.closed) putf(f, " %d", base);
            f << "\n";
            base += static_cast<int>(sd.samples.size());
        }
    }
    {
        std::ofstream f(dir + "/yarn_pieces.obj", std::ios::binary);
        f << "# " << title << ": yarn pieces [mm], one group per kind\n";
        int base = 1;
        for (const char* kind : {"head", "row", "cap", "binding", "pole_join"}) {
            putf(f, "g %s\n", kind);
            for (const YarnPiece& p : e.pieces) {
                if (std::string(toString(p.kind)) != kind) continue;
                for (const Vec3& q : p.samples) putf(f, "v %.6f %.6f %.6f\n", q.x, q.y, q.z);
                f << "l";
                for (size_t i = 0; i < p.samples.size(); ++i) putf(f, " %d", base + static_cast<int>(i));
                f << "\n";
                base += static_cast<int>(p.samples.size());
            }
        }
    }
    {
        std::ofstream f(dir + "/yarn_controls.csv", std::ios::binary);
        f << "strand,order,piece,kind,face,template,x,y,z,tx,ty,tz,depth_d\n";
        for (size_t k = 0; k < e.strands.size(); ++k) {
            int order = 0;
            for (int pi : e.strands[k].pieces) {
                const YarnPiece& p = e.pieces[pi];
                for (size_t i = 0; i < p.controls.size(); ++i) {
                    if (i == 0 && order > 0) continue;   // a shared port once
                    putf(f, "%zu,%d,%d,%s,%d,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.3f\n", k, order++, pi, toString(p.kind), p.face, p.templ.c_str(), p.controls[i].x,
                         p.controls[i].y, p.controls[i].z, p.tangents[i].x, p.tangents[i].y, p.tangents[i].z, p.depth[i]);
                }
            }
        }
    }
    writeTriObj(dir + "/yarn_tubes.obj", yarnTubes(e, 0.5 * st.diameter, 8));
    {
        std::ofstream f(dir + "/yarn_faces.csv", std::ios::binary);
        f << "face,kind,strip,slot,template,crossings_ok,expected,found,unexpected,min_depth_gap_d,min_clearance_mm,contact_samples,min_bend_mm,own_plane_ok,found_own_plane,"
             "own_plane_tilt_deg,bottom_mm,top_mm,wale_min_mm,wale_max_mm,shear_deg,non_planar_mm,head_outside_mm,row_outside_mm\n";
        for (const YarnFaceCheck& fc : c.faces)
            putf(f, "%d,%s,%d,%d,%s,%d,%s,%s,%d,%.4f,%.4f,%d,%.4f,%d,%s,%.2f,%.3f,%.3f,%.3f,%.3f,%.1f,%.3f,%.3f,%.3f\n", fc.face, toString(s.kind[fc.face]),
                 s.strip.empty() ? -1 : s.strip[fc.face], s.slot.empty() ? -1 : s.slot[fc.face], fc.templ.c_str(), fc.crossingsOk ? 1 : 0, fc.expected.c_str(), fc.found.c_str(),
                 fc.unexpectedCrossings, fc.minDepthGap < 1e299 ? fc.minDepthGap : -1.0, fc.minClearance < 1e299 ? fc.minClearance : -1.0, fc.contactSamples,
                 fc.minRadius < 1e299 ? fc.minRadius : -1.0, fc.ownPlaneOk ? 1 : 0, fc.foundOwnPlane.c_str(), fc.ownPlaneTilt, fc.bottomLength, fc.topLength, fc.waleMin, fc.waleMax,
                 fc.shear, fc.nonPlanar, fc.headOutside, fc.rowOutside);
    }
    {   // the diagnosis of a crossing order: every head-row crossing of every face, in the face's reference frame
        std::ofstream f(dir + "/yarn_crossings.csv", std::ios::binary);
        f << "face,kind,strip,slot,template,crossings_ok,head,row_piece,along_row,along_head,letter,depth_gap_d,u_mm,v_mm\n";
        for (const YarnFaceCheck& fc : c.faces)
            for (const YarnCrossing& x : fc.crossings)
                putf(f, "%d,%s,%d,%d,%s,%d,%d,%d,%.4f,%.4f,%c,%.4f,%.3f,%.3f\n", fc.face, toString(s.kind[fc.face]), s.strip.empty() ? -1 : s.strip[fc.face],
                     s.slot.empty() ? -1 : s.slot[fc.face], fc.templ.c_str(), fc.crossingsOk ? 1 : 0, x.head, x.rowPiece, x.alongRow, x.alongHead, x.letter, x.depthGap, x.u, x.v);
    }
    {
        std::ofstream f(dir + "/yarn_contacts.csv", std::ios::binary);
        f << "piece_a,kind_a,face_a,index_a,piece_b,kind_b,face_b,index_b,min_distance_mm,samples\n";
        for (const YarnCheck::Contact& ct : c.contacts) {
            const YarnPiece& pa = e.pieces[ct.pieceA];
            const YarnPiece& pb = e.pieces[ct.pieceB];
            putf(f, "%d,%s,%d,%d,%d,%s,%d,%d,%.4f,%d\n", ct.pieceA, toString(pa.kind), pa.face, pa.index, ct.pieceB, toString(pb.kind), pb.face, pb.index, ct.distance, ct.samples);
        }
    }
    {
        std::ofstream f(dir + "/yarn_check.txt", std::ios::binary);
        f << title << "\nsettings: diameter " << st.diameter << " mm, course ports 0.5 -/+ " << st.courseOffset << ", wale ports " << st.walePort << " (pole " << st.poleWalePort
          << "), cap height " << st.capHeight << " h, templates " << st.templateSet << ", placement " << toString(st.placement) << ", pole join " << (st.joinPoleRound ? "on" : "off") << ", binding thread " << (st.bindingThread ? "on" : "off") << "\n";
        for (const std::string& l : e.log) f << l << "\n";
        f << c.text();
        f << "defects (" << c.defects.size() << "):\n";
        for (const std::string& d : c.defects) f << "  " << d << "\n";
    }
    {
        std::ofstream f(dir + "/yarn_summary.csv", std::ios::binary);
        f << "key,value\n";
        putf(f, "status,%s\npieces,%d\nstrands,%d\nopen_strands,%d\nclosed_strands,%d\nyarn_ends,%d\nports_shared,%d\nports_unmatched,%d\nports_overused,%d\n", toString(c.status),
             c.pieces, c.strands, c.openStrands, c.closedStrands, c.yarnEnds, c.portsSharedOk, c.portsUnmatched, c.portsOverused);
        putf(f, "max_port_gap_mm,%.3g\nmax_tangent_jump_deg,%.3g\nmax_port_off_surface_mm,%.3g\nfaces_checked,%d\nfaces_crossings_ok,%d\nfaces_without_template,%d\n", c.maxPortGap,
             c.maxTangentJump, c.maxPortOffSurface, c.facesChecked, c.facesCrossingsOk, c.facesWithoutTemplate);
        putf(f, "interlocks_expected,%d\ninterlocks_verified,%d\nunexpected_crossings,%d\nmin_depth_gap_d,%.4f\nlength_mm,%.1f\nsamples,%d\n", c.interlocksExpected,
             c.interlocksVerified, c.unexpectedCrossings, c.minDepthGap < 1e299 ? c.minDepthGap : -1.0, c.totalLength, c.samples);
        putf(f, "min_clearance_mm,%.4f\ncontact_samples,%d\noverlap_samples,%d\nsevere_samples,%d\nmin_bend_mm,%.4f\nkink_samples,%d\ncaps,%d\npole_joined,%d\ndefects,%zu\n",
             c.minClearance < 1e299 ? c.minClearance : -1.0, c.contactSamples, c.overlapSamples, c.severeSamples, c.minRadius < 1e299 ? c.minRadius : -1.0, c.kinkSamples, e.caps,
             e.poleJoined ? 1 : 0, c.defects.size());
        for (int k = 0; k < 5; ++k)
            putf(f, "bend_tight_%s,%d\nsamples_%s,%d\n", toString(static_cast<YarnPieceKind>(k)), c.kinkByKind[k], toString(static_cast<YarnPieceKind>(k)), c.samplesByKind[k]);
    }
}

// ----------------------------------------------------------------------------------------------------------- fixtures
// Builds, relaxes (CP3 settings) and returns the fixture's typed mesh, virtual faces and reference.
bool buildFixture(const YarnFixture& fx, double w, double h, TypedStitchMesh& m, std::vector<std::vector<int>>& faces, TriMesh& ref, RelaxProduct& relax, std::string* err) {
    TriMesh unusedRef;
    std::vector<DeclaredBoundary> decl;
    std::vector<double> rho, z;
    if (fx.pole) {
        // the arm's pole-start fixture: rings k h around the pole in a plane; the reference a disc ending at the last ring
        if (!poleDiscFixture(fx.counts, fx.shaping, h, m, unusedRef, decl, err)) return false;
        const int K = static_cast<int>(fx.counts.size());
        for (int k = 0; k <= K; ++k) {
            rho.push_back(k == 0 ? 1e-3 : k * h);
            z.push_back(0);
        }
        ref = revolutionThroughRings(rho, z, 192, 4);
        // the plane z = 0 seen from +z: the fixture's faces are counter-clockwise from +z, as the reference's must be
        bool flip = false;
        for (const Vec3& n : ref.vnormal)
            if (n.z < 0) flip = true;
        if (flip) {
            for (auto& t : ref.tri) std::swap(t[1], t[2]);
            ref.computeNormals();
        }
    } else {
        if (!ringTubeFixture(fx.counts, fx.shaping, 20, h, m, unusedRef, decl, err, fx.boundaryShaping)) return false;
        for (size_t k = 0; k < fx.counts.size(); ++k) {
            rho.push_back(fx.counts[k] * w / (2 * kPi));
            z.push_back(k == 0 ? 0.0 : z.back() + std::sqrt(std::max(0.0, h * h - (rho[k] - rho[k - 1]) * (rho[k] - rho[k - 1]))));
        }
        ref = revolutionThroughRings(rho, z, 192, 4);
        const MeshQuery q(ref);
        for (size_t k = 0; k < m.rings.size(); ++k)
            for (int v : m.rings[k]) {
                const double phi = std::atan2(m.mesh.pos[v].y, m.mesh.pos[v].x);
                m.mesh.pos[v] = q.closest(Vec3{rho[k] * std::cos(phi), rho[k] * std::sin(phi), z[k]}).p;
            }
    }
    VirtualShift vs;
    if (!buildVirtualShift(m, ShiftReading::PathSideDown, vs, err)) return false;
    faces = vs.corners;
    std::vector<RelaxFace> rf;
    for (int f = 0; f < m.mesh.faces(); ++f) rf.push_back(RelaxFace{vs.corners[f], m.role[f]});
    std::vector<char> border(m.mesh.pos.size(), 0);
    for (const std::vector<int>& l : polygonBorderLoops(m.mesh))
        for (int v : l) border[v] = 1;
    RelaxGuides guides;
    const RelaxGuides* g = nullptr;
    if (fx.pole) {
        guides.lineOf.assign(m.mesh.pos.size(), -1);
        guides.pinned.assign(m.mesh.pos.size(), 0);
        guides.pinned[m.rings[0][0]] = 1;
        g = &guides;
    }
    MeshRelaxRequest r;
    r.w = w;
    r.h = h;
    if (!relaxStitchMesh(ref, triangleBoundaryLoops(ref), rf, m.mesh.pos, border, r, relax, err, g)) return false;
    m.mesh.pos = relax.pos;
    return true;
}

std::vector<YarnFixture> fixtureCases() {
    return {
        {"K_patch", "regular stockinette: a tube of 12 x 3 rows, cast-on below, bind-off above", false, {12, 12, 12, 12}, {{}, {}, {}}},
        {"K_interior", "regular stockinette for the interior comparison (template calibration): a tube of 24 x 6 rows; strips 2-3 are two rows from either boundary",
         false, {24, 24, 24, 24, 24, 24, 24}, {{}, {}, {}, {}, {}, {}}},
        {"INC_isolated", "one kp in an interior strip (12 -> 13)", false, {12, 12, 13, 13}, {{}, {6}, {}}},
        {"DEC_isolated", "one d12k in an interior strip (13 -> 12)", false, {13, 13, 12, 12}, {{}, {6}, {}}},
        {"kp_cast_on", "boundary shaping: kp + cast-on (12 -> 13 in the first strip)", false, {12, 13, 13, 13}, {{6}, {}, {}}, true},
        {"d12k_bind_off", "boundary shaping: d12k + bind-off (13 -> 12 in the last strip)", false, {13, 13, 13, 12}, {{}, {}, {6}}, true},
        {"pole_start", "gathered pole start: pole, ring of 5, 5 -> 7 (kp at slots 1, 3), 7 -> 7, bind-off", true, {5, 7, 7}, {{1, 3}, {}}},
        {"pole_start_adjacent", "gathered pole start + adjacent shaping: 5 -> 8 with kp at slots 1, 2, 3", true, {5, 8, 8}, {{1, 2, 3}, {}}},
        {"INC_every_other", "dense: five kp on every other stitch (12 -> 17)", false, {12, 12, 17, 17}, {{}, {1, 3, 5, 7, 9}, {}}},
        {"INC_adjacent_pairs", "dense: two pairs of neighbouring kp (16 -> 20)", false, {16, 16, 20, 20}, {{}, {2, 3, 9, 10}, {}}},
        {"DEC_every_other", "dense: five d12k on every other stitch (17 -> 12)", false, {17, 17, 12, 12}, {{}, {1, 3, 5, 7, 9}, {}}},
        {"neg_flipped_head", "negative control: K_patch with one head's depth flipped (strip 1, slot 4): the crossing check must fail there", false, {12, 12, 12, 12},
         {{}, {}, {}}, false, true, 4, -1},
        {"neg_missing_face", "negative control: K_patch with one face left without yarn (strip 1, slot 4): unmatched ports, extra strands", false, {12, 12, 12, 12},
         {{}, {}, {}}, false, true, -1, 4},
    };
}

int commandYarnFixtures(const Opts& a) {
    const std::string out = a.get("--out");
    fs::create_directories(out);
    const double w = 5.5, h = 3.4;
    const YarnSettings base = settingsFrom(a);
    const std::vector<YarnFixture> cases = lab::yarnFixtureCases();
    std::ofstream csv(out + "/yarn_fixtures.csv", std::ios::binary);
    csv << "case,negative,relax_status,relax_steps,status,strands,open,closed,yarn_ends,ports_unmatched,ports_overused,faces_ok,faces,interlocks_verified,interlocks_expected,"
           "unexpected,min_depth_gap_d,max_port_gap_mm,max_tangent_jump_deg,min_clearance_mm,contact_pct,overlap_samples,severe_samples,min_bend_mm,kinks,as_expected\n";
    std::ofstream rep(out + "/yarn_fixtures.txt", std::ios::binary);
    putf(rep, "gauge %.1f x %.1f, yarn diameter %.2f mm; fixtures relaxed with the CP3 model (pole start: pole pinned) before the yarn is embedded\n", w, h, base.diameter);
    for (const YarnFixture& fx : cases) {
        TypedStitchMesh m;
        std::vector<std::vector<int>> faces;
        TriMesh ref;
        RelaxProduct relax;
        std::string err;
        if (!buildFixture(fx, w, h, m, faces, ref, relax, &err)) {
            putf(rep, "%s: not built: %s\n", fx.name.c_str(), err.c_str());
            continue;
        }
        YarnSettings st = base;
        for (int f = 0; f < m.mesh.faces(); ++f) {
            if (!m.strip.empty() && m.strip[f] == 1 && m.slot[f] == fx.flipHeadsOfSlot) st.flipHeadsOfFace = f;
            if (!m.strip.empty() && m.strip[f] == 1 && m.slot[f] == fx.skipSlot) st.skipFace = f;
        }
        YarnEmbedding e;
        if (!embedYarn(m, faces, ref, w, h, st, e, &err)) {
            putf(rep, "%s: no yarn: %s\n", fx.name.c_str(), err.c_str());
            continue;
        }
        const YarnCheck c = checkYarn(m, faces, ref, e, st);
        const std::string dir = out + "/" + fx.name;
        writeYarn(dir, m, e, c, st, fx.name + ": " + fx.purpose);
        lab::writePolygonObj(dir + "/stitch_mesh.obj", m);
        writeTriObj(dir + "/reference.obj", ref);
        const bool asExpected = fx.negative ? c.status == Status::Fail : c.status == Status::Pass;
        putf(csv, "%s,%d,%s,%d,%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.3f,%.3g,%.3g,%.4f,%.2f,%d,%d,%.4f,%d,%d\n", fx.name.c_str(), fx.negative ? 1 : 0, toString(relax.status), relax.steps,
             toString(c.status), c.strands, c.openStrands, c.closedStrands, c.yarnEnds, c.portsUnmatched, c.portsOverused, c.facesCrossingsOk, c.facesChecked,
             c.interlocksVerified, c.interlocksExpected, c.unexpectedCrossings, c.minDepthGap < 1e299 ? c.minDepthGap : -1.0, c.maxPortGap, c.maxTangentJump,
             c.minClearance < 1e299 ? c.minClearance : -1.0, c.samples ? 100.0 * c.contactSamples / c.samples : 0.0, c.overlapSamples, c.severeSamples,
             c.minRadius < 1e299 ? c.minRadius : -1.0, c.kinkSamples, asExpected ? 1 : 0);
        putf(rep, "\n== %s (%s)%s: relaxation %s (%d steps); %s", fx.name.c_str(), fx.purpose.c_str(), asExpected ? "" : "  ** NOT AS EXPECTED **", toString(relax.status),
             relax.steps, c.text().c_str());
        for (size_t i = 0; i < c.defects.size() && i < 8; ++i) putf(rep, "   defect: %s\n", c.defects[i].c_str());
        if (c.defects.size() > 8) putf(rep, "   ... %zu defects in all (yarn_check.txt)\n", c.defects.size());
    }
    std::printf("yarn-fixtures: written to %s\n", out.c_str());
    return 0;
}

// --------------------------------------------------------------------------------------------------- one frozen mesh
int commandYarnEmbed(const Opts& a) {
    const std::string caseDir = a.get("--case"), which = a.get("--which", "after"), out = a.get("--out");
    const double w = a.number("--w", 5.5), h = a.number("--h", 3.4);
    TypedStitchMesh s;
    std::string saved;
    if (!lab::readTypedMeshFiles(caseDir + "/" + which, s, saved)) {
        std::fprintf(stderr, "yarn-embed: cannot read the typed mesh in %s/%s\n", caseDir.c_str(), which.c_str());
        return 1;
    }
    TriMesh ref;
    lab::readReferenceObj(caseDir + "/reference.obj", ref);
    if (ref.tri.empty()) {
        std::fprintf(stderr, "yarn-embed: no reference.obj in %s\n", caseDir.c_str());
        return 1;
    }
    char digest[32];
    std::snprintf(digest, sizeof digest, "%016llx", static_cast<unsigned long long>(typedDigest(s)));
    VirtualShift vs;
    std::string err;
    if (!buildVirtualShift(s, ShiftReading::PathSideDown, vs, &err)) {
        std::fprintf(stderr, "yarn-embed: %s\n", err.c_str());
        return 1;
    }
    const YarnSettings st = settingsFrom(a);
    YarnEmbedding e;
    if (!embedYarn(s, vs.corners, ref, w, h, st, e, &err)) {
        std::fprintf(stderr, "yarn-embed: %s\n", err.c_str());
        return 1;
    }
    const YarnCheck c = checkYarn(s, vs.corners, ref, e, st);
    fs::create_directories(out);
    writeYarn(out, s, e, c, st, "frozen mesh " + caseDir + "/" + which + " (typed digest " + digest + ", saved " + saved + ")");
    lab::writePolygonObj(out + "/stitch_mesh.obj", s);
    std::ofstream f(out + "/yarn_check.txt", std::ios::app | std::ios::binary);
    if (saved.empty())
        putf(f, "frozen mesh: typed digest as read %s (the run saved none; the driver records the input files' sha256)\n", digest);
    else
        putf(f, "frozen mesh: typed digest as read %s, as saved %s: %s\n", digest, saved.c_str(), saved == digest ? "unchanged" : "DIFFERENT");
    std::printf("yarn-embed: %s", c.text().c_str());
    std::printf("frozen mesh digest %s (saved %s)\n", digest, saved.c_str());
    return 0;
}

} // namespace

int yarnCommand(const std::string& cmd, int argc, char** argv) {
    const Opts a = parse(argc, argv);
    if (cmd == "yarn-fixtures") return commandYarnFixtures(a);
    if (cmd == "yarn-embed") return commandYarnEmbed(a);
    return 2;
}

namespace knit::lab {

std::vector<YarnFixture> yarnFixtureCases() { return fixtureCases(); }
bool buildYarnFixture(const YarnFixture& fx, double w, double h, TypedStitchMesh& m, std::vector<std::vector<int>>& faces, TriMesh& ref, RelaxProduct& relax, std::string* err) {
    return buildFixture(fx, w, h, m, faces, ref, relax, err);
}
YarnSettings yarnSettingsFrom(const std::map<std::string, std::string>& kv) {
    Opts o;
    o.kv = kv;
    return settingsFrom(o);
}
void writeYarnFiles(const std::string& dir, const TypedStitchMesh& s, const YarnEmbedding& e, const YarnCheck& c, const YarnSettings& st, const std::string& title) {
    writeYarn(dir, s, e, c, st, title);
}
bool readReferenceObj(const std::string& path, TriMesh& ref) { return readObjTriangles(path, ref); }
void writeTriangleObj(const std::string& path, const TriMesh& m) { writeTriObj(path, m); }

} // namespace knit::lab
