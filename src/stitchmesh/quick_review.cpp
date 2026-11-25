// The quick review while drawing (quick_review.h).
#include "stitchmesh/quick_review.h"

#include <chrono>
#include <cmath>
#include <map>

namespace knit {

namespace {

// Where two segments of one trace cross inside the same triangle (in its plane): the first such place, if any. Segment i
// runs from point i to point i + 1 inside triangle points[i].tri; neighbouring segments share a point and are skipped.
bool selfCrossing(const TriMesh& m, const MarkingTrace& tr, bool closed, int& first, int& second, Vec3& at) {
    std::map<int, std::vector<int>> segmentsOf;
    const int n = static_cast<int>(tr.points.size());
    for (int i = 0; i + 1 < n; ++i)
        if (tr.points[i].tri >= 0) segmentsOf[tr.points[i].tri].push_back(i);
    for (const std::pair<const int, std::vector<int>>& entry : segmentsOf) {
        const std::vector<int>& seg = entry.second;
        if (seg.size() < 2) continue;
        const std::array<int, 3>& t = m.tri[entry.first];
        const Vec3 o = m.pos[t[0]];
        const Vec3 e1 = normalized(m.pos[t[1]] - o);
        const Vec3 e2 = normalized(cross(cross(e1, m.pos[t[2]] - o), e1));
        auto flat = [&](const Vec3& x) { return Vec3{dot(x - o, e1), dot(x - o, e2), 0}; };
        for (size_t a = 0; a < seg.size(); ++a)
            for (size_t b = a + 1; b < seg.size(); ++b) {
                const int i = seg[a], j = seg[b];
                if (std::abs(i - j) <= 1 || (closed && i == 0 && j == n - 2)) continue;   // neighbours share a point
                const Vec3 p = flat(tr.points[i].pos), r = flat(tr.points[i + 1].pos) - p;
                const Vec3 q = flat(tr.points[j].pos), s = flat(tr.points[j + 1].pos) - q;
                const double den = r.x * s.y - r.y * s.x;
                if (std::fabs(den) < 1e-18) continue;
                const Vec3 d = q - p;
                const double u = (d.x * s.y - d.y * s.x) / den, v = (d.x * r.y - d.y * r.x) / den;
                if (u <= 0 || u >= 1 || v <= 0 || v >= 1) continue;
                first = i;
                second = j;
                at = tr.points[i].pos + u * (tr.points[i + 1].pos - tr.points[i].pos);
                return true;
            }
    }
    return false;
}

// A problem or an unfinished state among the first `count` findings (those of the resolution: the marking cannot be
// evaluated as drawn) on this marking.
bool problemOn(const MarkingReview& r, size_t count, int markingId) {
    for (size_t i = 0; i < count && i < r.findings.size(); ++i) {
        const MarkingFinding& f = r.findings[i];
        if (f.markingId == markingId && (f.level == FindingLevel::Problem || f.level == FindingLevel::Incomplete)) return true;
    }
    return false;
}

} // namespace

MarkingReview quickReview(const QuickRequest& q, std::vector<int>* markingIds, std::vector<MarkingTrace>* traces) {
    const std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
    MarkingReview review;
    review.depth = ReviewDepth::Quick;
    review.measures.secondsKnown = Known::Current;
    if (!q.reference || !q.markings || !q.pieces || q.piece < 0 || q.piece >= static_cast<int>(q.pieces->pieces.size())) return review;
    // 1. resolution, per-marking checks and the choice: the build's own step 1
    MarkedBuildRequest rq;
    rq.pieces = q.pieces;
    rq.piece = q.piece;
    rq.knitComponent = q.knitComponent;
    rq.markings = q.markings;
    rq.mmPerUnit = q.mmPerUnit;
    rq.w = q.w;
    rq.h = q.h;
    rq.previewOffered = q.previewOffered;
    ResolvedMarkings rm;
    resolvePieceMarkings(rq, q.reference->refined, q.reference->nb, review, rm);
    const size_t resolutionFindings = review.findings.size();   // a self-crossing (below) still lets the rows be counted
    // 2. a path or curve that crosses itself (an opening's loop never does)
    const double toFile = 1.0 / q.mmPerUnit;
    for (const ResolvedMarking& u : rm.used) {
        if (!u.resolved || !u.borderLoop.empty() || u.marking->shape == MarkingShape::Point || !markingComplete(*u.marking)) continue;
        int first = -1, second = -1;
        Vec3 at;
        if (!selfCrossing(q.reference->refined.mesh, u.trace, u.marking->shape == MarkingShape::ClosedCurve, first, second, at)) continue;
        MarkingFinding f;
        f.markingId = u.marking->id;
        f.level = FindingLevel::Problem;
        f.cause = FindingCause::Marking;
        f.message = "This marking crosses itself here; move a point so that it does not loop back over itself.";
        f.firstTracePoint = first;
        f.lastTracePoint = second + 1;
        f.at.push_back(at * toFile);
        review.findings.push_back(f);
    }
    reviewMeaningChoice(rm, review);
    // 3. the rows, only from a field of exactly this start, finish, gauge and reference
    const MeaningChoice& c = rm.choice;
    const bool complete = c.start >= 0 && c.finish >= 0 && c.rowStart >= 0;
    if (complete) {
        const int startId = rm.used[c.start].marking->id, finishId = rm.used[c.finish].marking->id, pathId = rm.used[c.rowStart].marking->id;
        const uint64_t key = markedFieldKey(rm, q.markings->digest, toString(q.preparation), q.piece, q.mmPerUnit, q.h);
        const bool current = q.field && q.field->key == key && q.field->reference && q.field->reference->key == q.reference->key;
        const bool usable = !problemOn(review, resolutionFindings, startId) && !problemOn(review, resolutionFindings, finishId) &&
                            !problemOn(review, resolutionFindings, pathId);
        if (current && usable) {
            const std::vector<RowCrossing> rows = rowCrossings(*q.field, rm.used[c.rowStart].trace);
            crossingFindings(rows, q.field->route, startId, pathId, toFile, review);
            int events = 0;
            for (const RowCrossing& rc : rows) events += rc.split || rc.count != 1 ? 1 : 0;
            review.measures.rows = q.field->rows;
            review.measures.contourEvents = events;
            review.measures.rowsKnown = review.measures.contourEventsKnown = Known::Current;
        } else {
            review.measures.rowsKnown = review.measures.contourEventsKnown = current || !q.field ? Known::No : Known::Outdated;
            MarkingFinding f;
            f.level = FindingLevel::Ok;
            f.cause = FindingCause::Input;
            f.message = !usable  ? "Row checks wait until the start, finish and row-start path can be traced (see the problems above)."
                        : q.field ? "Row checks are pending: the start, finish, gauge or reference changed since the last layout. Build the layout (or let it update) to check the rows."
                                  : "Row checks need a layout of this start and finish: build the layout (or let it update) to check the rows.";
            review.findings.push_back(f);
        }
    }
    if (markingIds) *markingIds = rm.markingIds;
    if (traces) *traces = rm.traces;
    review.measures.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    return review;
}

} // namespace knit
