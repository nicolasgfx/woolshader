#pragma once
// The first yarn-generation experiment: yarn centrelines embedded
// in a typed stitch mesh, face by face, and checked as yarn. No thickness, contact or relaxation yet: this is the
// construction the next checkpoint's yarn solver will start from.
//
// Face semantics (Yuksel et al. 2012 section 4, Fig. 3): a face holds the top of one loop and the bottom of the next. Per
// face, from canonical templates (yarn_embed.cpp, proposed):
//   head      per bottom edge (not cast-on): the older loop, entering and leaving through that edge's two ports
//   row       the face's own row yarn: from the W- port up through each top edge's two ports and out of the W+ port
//             (K: two pieces; kp: three, the second the window between the knit and the purl slot)
//   cap       per bind-off top edge: the last row's loop closing above the boundary (proposed)
// plus, once per component: the binding thread (proposed: the yarn continues from the helix's end through every cap once)
// and the pole-round join (proposed: the gathered pole round opens at the path vertex and continues into strip 1).
//
// Ports. A course edge has two, at 0.5 -/+ courseOffset along the edge from its W- side vertex; a wale edge one, at walePort
// from its bottom vertex. Both faces of an edge place a port from the same vertex identities, so they compute the same
// point (projected onto the reference) and the same tangent (perpendicular to the edge in the reference's tangent plane):
// position and direction continuity at ports hold by construction. Interior controls are canonical (u, v, z) points:
// mean value coordinates on the canonical polygon (Floater 2003) give the foot on the face's corners, projected onto the
// reference, offset z D along its normal (Y12 section 4: "mean value coordinates" plus "an offset value in the surface
// normal direction"). Each piece is a cubic Hermite curve through its controls: the port tangents shared, the interior
// ones Catmull-Rom.
//
// The faces are the virtual (helical) polygons of the shift (virtual_shift.h): the strand follows the helix. The start
// step's collapsed W- corner is where the yarn starts (or where the pole-round join arrives).
#include <cstdint>
#include <string>
#include <vector>

#include "model/mesh.h"
#include "model/status.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"

namespace knit {

// Where a face's yarn is placed between its corners.
enum class KnotPlacement {
    // The first-yarn rule (kept for comparison): the mean value foot on the face's corners is moved to the closest point of the reference and offset
    // along the reference normal there. On a face that folds across a crease of the reference (the row-start face of a
    // strip whose stitch count changes spans two rows) the closest points bunch up and the head leaves its face.
    ReferenceProjection,
    // Yuksel et al. 2012 section 4 as written: the mean value foot on the face's corners, offset along the corner normals
    // interpolated with the same weights. Ports the same way on their edge. The yarn follows the face, so the map from
    // the template to the face stays one to one however the reference bends under it.
    FaceInterpolation,
};
const char* toString(KnotPlacement p);

struct YarnSettings {
    double diameter = 1.16;       // D [mm]: z offsets are in units of D; the tubes' radius is D / 2. Production knits this gauge at 1.16 (thickness_presets.cpp: the contact diameter)
    double courseOffset = 0.2;    // course ports at 0.5 -/+ this along the edge from its W- side vertex
    double walePort = 0.35;       // wale ports at this from the bottom vertex
    double poleWalePort = 0.5;    // ... on a gathered pole start's wale edges (from the pole)
    double poleTailPort = 0.25;   // the pole round's yarn end on its join edge (from the pole)
    double capHeight = 0.45;      // bind-off caps rise this x h above their edge
    int samplesPerSpan = 8;       // samples per Hermite span (checks, tubes)
    // Interior tangents: false (frozen): Catmull-Rom, 0.5 (P[i+1] - P[i-1]) on both spans; true (a variant tested after the run):
    // Catmull-Rom's direction scaled on each span by that span's chord (fewer hooks where controls are unevenly spaced)
    bool spanTangents = false;
    bool joinPoleRound = true;    // the gathered pole round continues into strip 1 (else it stays a closed strand)
    bool bindingThread = true;    // the yarn continues from the helix's end through every cap (else the caps stay live)
    // The template set: "baseline" (the first-yarn templates, frozen) or "calibrated"
    // (the candidate). Calibrated changes the K template's knots (depth +-0.65 D, production's clasp proportions, longer depth
    // flips), kp's and d12k's depth (x 1.3; d12k's two parent heads separated) and the shared ports of every edge the rule in
    // embedYarn admits: wale ports calibratedWaleDepth x D along the normal, course ports calibratedCourseDepth x D, course
    // tangents leaning calibratedLean degrees towards the W- side. Both faces of an edge derive the same port. Pole faces, kp's
    // purl slot and W+, and the B edges of kp and d12k keep the frozen ports.
    // "mixed": calibrated K and kp, baseline d12k with frozen ports on every edge of a d12k face.
    std::string templateSet = "baseline";
    KnotPlacement placement = KnotPlacement::FaceInterpolation;
    double calibratedWaleDepth = -0.5, calibratedCourseDepth = 0.5, calibratedLean = 35.0;
    // negative controls only (the checks must fail with the defect back): flip the depth of this face's heads; leave this
    // face without yarn
    int flipHeadsOfFace = -1;
    int skipFace = -1;
};

enum class YarnPieceKind : char { Head, Row, Cap, Binding, PoleJoin };
const char* toString(YarnPieceKind k);

// A port: where pieces meet. Course ports (kind 0) by (W- side vertex, other vertex, index 0 / 1); wale ports (kind 1) by
// (bottom, top); a collapsed corner (kind 2) by its vertex; a yarn end that meets nothing (kind 3) by a unique number.
struct YarnPort {
    int kind = 3, a = -1, b = -1, index = 0;
    bool operator<(const YarnPort& o) const;
    bool operator==(const YarnPort& o) const { return kind == o.kind && a == o.a && b == o.b && index == o.index; }
};

struct YarnPiece {
    int face = -1;                // the face it lies in (-1: the binding thread)
    YarnPieceKind kind = YarnPieceKind::Row;
    int index = 0;                // head: which bottom edge (0, 1); row: which segment; cap: which top edge
    std::string templ;            // the template it came from
    YarnPort start, end;          // in yarn order
    std::vector<Vec3> controls;   // Hermite knots [mm]
    std::vector<Vec3> tangents;   // at the knots [mm]
    std::vector<double> depth;    // z of each knot [D]
    std::vector<Vec3> samples;    // the curve, samplesPerSpan per span
};

struct YarnStrand {
    std::vector<int> pieces;      // in yarn order
    bool closed = false;
    std::vector<Vec3> samples;    // concatenated (the shared port sample once)
    std::vector<int> sampleFace;  // per sample: its piece's face
    std::vector<int> samplePiece;
    std::vector<double> arc;      // arc length per sample [mm]
};

struct YarnEmbedding {
    std::vector<YarnPiece> pieces;
    std::vector<YarnStrand> strands;
    std::vector<std::string> faceTemplate;   // per face ("" when the face has none)
    std::vector<int> facesWithoutTemplate;
    bool poleJoined = false;
    int caps = 0;
    std::vector<std::string> log;
};

// faces: the virtual corners (VirtualShift::corners) with the stored kinds, roles and boundary roles. reference: the surface
// the mesh lies on (mm, vertex normals). False only for an unusable input (no faces, sizes disagree).
bool embedYarn(const TypedStitchMesh& s, const std::vector<std::vector<int>>& faces, const TriMesh& reference, double w, double h, const YarnSettings& settings,
               YarnEmbedding& out, std::string* err);

// One projected crossing of a face's row yarn with one of its heads, in the face's reference frame: what the letters of
// YarnFaceCheck::found are made of. Diagnosis of a wrong crossing order; no verdict reads it.
struct YarnCrossing {
    int head = 0;           // which head (the Head piece's YarnPiece::index)
    int rowPiece = 0;       // which row segment (the Row piece's YarnPiece::index)
    double alongRow = 0;    // position along that row segment: 0 at its start port, 1 at its end port
    double alongHead = 0;   // position along the head: 0 to 1
    char letter = 'B';      // the row in front of the head (F) or behind it (B), along the frame's normal
    double depthGap = 0;    // |height difference| [D]
    double u = 0, v = 0;    // where, in the face frame [mm] (origin: the face's centre on the reference, u along its first edge)
};

// Per face: the crossings its yarn makes in projection on the face's own plane (the local frame at its centroid, the
// reference normal), against what its construction requires.
struct YarnFaceCheck {
    int face = -1;
    std::string templ;
    std::string expected, found;  // e.g. "head0:BFFB head1:BFFB h1/h0:FF" against the computed sequences
    bool crossingsOk = true;
    // The same rule in the frame of the yarn as it is (the yarn-settle continuation): the plane fitting the face's head and row
    // samples best, instead of the reference's tangent plane. After the yarn moved against the reference, the reference frame
    // can misrepresent the 3D state; equal verdicts in both frames are consistent evidence, a difference is shown for
    // inspection. Neither frame proves that no strand passed through another (the motion certificate is that evidence).
    bool ownPlaneOk = false;
    std::string foundOwnPlane;
    double ownPlaneTilt = 0;      // angle between the two frames' normals [deg]
    int unexpectedCrossings = 0;  // crossings between pieces of the face that its construction does not have
    double minDepthGap = 1e300;   // smallest |z difference| at its intended crossings [D]
    double minClearance = 1e300;  // nearest other yarn from its samples [mm]
    double minRadius = 1e300;     // tightest bend of its samples [mm]
    int contactSamples = 0;       // samples whose nearest other yarn is closer than D
    // Diagnosis only (no verdict reads these): the head-row crossings one by one, and the shape of the virtual face the
    // template was mapped into, both in the reference frame.
    std::vector<YarnCrossing> crossings;
    double bottomLength = 0, topLength = 0;   // the bottom edges and the top edges, each summed [mm]
    double waleMin = 0, waleMax = 0;          // the shorter and the longer wale edge [mm]
    double shear = 0;             // angle between the bottom-to-top direction (edge midpoints) and the bottom edges' normal [deg]
    double nonPlanar = 0;         // largest distance of a corner from the frame's plane [mm]
    double headOutside = 0;       // how far a head sample lies outside the face's polygon in the frame [mm] (0: all inside)
    double rowOutside = 0;        // the same for the row samples
};

struct YarnCheck {
    Status status = Status::Unresolved;
    // the strands
    int pieces = 0, strands = 0, openStrands = 0, closedStrands = 0, yarnEnds = 0;
    int portsSharedOk = 0, portsUnmatched = 0, portsOverused = 0;   // ports with one piece in and one out; with one side only (not a yarn end); with more
    double maxPortGap = 0;        // between consecutive pieces of a strand [mm]
    double maxTangentJump = 0;    // between consecutive pieces of a strand [deg]
    double maxPortOffSurface = 0; // ports to the reference [mm]
    // crossings and interlocks
    int facesChecked = 0, facesCrossingsOk = 0, facesWithoutTemplate = 0;
    int facesOwnPlaneOk = 0;      // faces whose crossings are as required in the yarn's own plane (YarnFaceCheck::ownPlaneOk)
    int interlocksOwnPlane = 0;   // interlocks with the required order in that plane
    int interlocksExpected = 0, interlocksVerified = 0;   // heads (older loops through a bottom edge) with the required order
    int unexpectedCrossings = 0;
    double minDepthGap = 1e300;
    // geometry of the centrelines
    double minClearance = 1e300;  // between yarn more than 3 D apart along the yarn [mm]
    int samples = 0, contactSamples = 0, overlapSamples = 0, severeSamples = 0;   // nearest other yarn < D, < 0.8 D, < 0.5 D
    double minRadius = 1e300;     // tightest bend [mm]
    int kinkSamples = 0;          // bends tighter than D / 2
    int kinkByKind[5] = {0, 0, 0, 0, 0};    // per YarnPieceKind
    int samplesByKind[5] = {0, 0, 0, 0, 0};
    double totalLength = 0;       // [mm]
    struct Contact { int pieceA = -1, pieceB = -1; double distance = 0; int samples = 0; };
    std::vector<Contact> contacts;      // piece pairs closer than 0.8 D, the closest first (pieceA < pieceB)
    std::vector<YarnFaceCheck> faces;
    std::vector<std::string> defects;   // one line each, the worst first
    std::string text() const;
};

YarnCheck checkYarn(const TypedStitchMesh& s, const std::vector<std::vector<int>>& faces, const TriMesh& reference, const YarnEmbedding& e, const YarnSettings& settings);

// Tubes of radius r around every strand (rotation-minimising frames, `sides` vertices per ring).
TriMesh yarnTubes(const YarnEmbedding& e, double r, int sides);

// The helper the yarn-settle adapter (yarn_curves.h) shares with the embedding: a piece's Hermite curve (per samples per
// span, the last knot once; spanScaled as YarnSettings::spanTangents).
std::vector<Vec3> sampleYarnHermite(const std::vector<Vec3>& knots, const std::vector<Vec3>& tangents, int per, bool spanScaled);

} // namespace knit
