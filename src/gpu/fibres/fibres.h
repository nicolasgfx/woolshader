#pragma once
// GPU fibres: the flat parameter block and the per-candidate records of the generator, shared by
// fibres.cu and fibre_candidate.cuh.
#include <cstdint>

#include "fibres/fibre.h"
#include "fibres/fibre_record.h"
#include "gpu/geometry/device_geometry.h"

namespace knit::gpu {

// Maximal vertex count of a polyline fibre: root + anchor + nv (nv <= maxSegments). The host
// refuses settings with 2 + max(maxSegments, max(2, minSegments)) > kMaxFibreVerts.
constexpr int kMaxFibreVerts = 50;

// What became of a candidate: knit::FibreOutcome (fibres/fibre.h), the CPU twin's list, stored as a byte per
// candidate. NotProcessed is never left in a finished run.
__host__ __device__ constexpr uint8_t outcomeCode(FibreOutcome o) { return static_cast<uint8_t>(o); }
// The key is outside [keyLo, keyHi): not a candidate of this run, not counted. Outside FibreOutcome's range.
constexpr uint8_t kOutOfWindowCode = outcomeCode(FibreOutcome::Count);

struct FibreParamsFlat {
    // FibreParams
    double density;
    int budget;              // FibreFrame::candidates given (s = (i + jitter) L / n)
    float frontCos, keepProbability, curvatureWeight;
    float length, lengthMean, lengthMin, lengthMax, lengthSpread, thickness, thicknessSigma, taper, rootDepth;
    float helixDeg, helixSigmaDeg, liftMinDeg, liftMaxDeg, bowTheta, bowRadial;
    float hugHeight, plyPitch;  // FibreParams::hugHeight, plyPitch
    float segLenTarget, maxCurveDev, minTurnDeg, maxTurnDeg, attachMax, sampleStepFactor, serializedMargin,
        exportChordFactor, rootTurnDeg, minEmergeDeg, plyOuterCos;
    float randomFraction, randomElevMinDeg, randomElevMaxDeg, randomBend;
    float randomSpreadDeg;  // FibreParams::randomSpreadDeg (180: any direction)
    float rLengthMean, rLengthMin, rLengthMax, rLengthSpread, rThickness, rThicknessSigma;  // populationParams of Random
    float rCurlRadius, rCurlLength;  // FibreParams::randomCurlRadius, randomCurlLength
    float crimpAmplitude, crimpWavelength;  // FibreParams::crimpAmplitude, crimpWavelength
    int handedness, minSegments, maxSegments;
    int hasRandomOverrides, lipschitz;
    int light;  // FibreParams::lightChecks
    int drape;  // FibreParams::drape (light checks only)
    double randomFractionD;  // exact population draw (u17 < fraction)
    uint64_t seed;           // fibreSeed (piece-mixed)
    float keyLo, keyHi;      // budget order: only candidates with float key in [keyLo, keyHi) (0, 2: all)
    // pattern / frame
    float r;                 // pp.radius() of the supporting strands
    float ownWindow;         // 3 D
    float kappaRef;          // 1 / (4 r)
    int surf;                // frame->outward with exposureLength > 0
    float exposureLength;
    float backDepth;         // curved: FibreFrame::backDepth
    float zTop, zBottom;     // flat: SurfaceIndex::zTop, zBottom
    int hidden;              // frame->hidden present
    int nHidden;
    float hidScale, hidOffset[3];  // common frame of the hidden meshes: q = p * scale + offset
    int quantize;
    float qInvSigma, qSigma, qCentre[3];  // the export's rounding: library -> scene -> float -> library
    int parents;             // FibreFrame::plies (plies)
};

// One supporting strand's candidate range: candidates [first, first + n) of the run.
struct FibreStrand {
    long long first;
    int n;
    int yarn;               // index into the supporting curve set
};

// Kernel A survivors / kernel B input: the candidate's global index and anchor data.
struct FibreSeed {
    long long cand;         // global candidate index
    int strand;             // into FibreStrand[]
    int i;                  // candidate index on the strand
    double s;               // anchor arc length
    float theta;
    int parentYarn;
};

struct FibreViews {
    CurveSetView cs;            // supporting strands (yarns or plies)
    SurfaceView sm;
    SurfaceIndexView index;
    CurveSetView parentCs;      // plies: the parent yarns
    SurfaceView parentSm;
    const int* plyYarnOf;       // per ply: parent yarn
    const double* plyOffset;    // per ply: parent arc length at u = 0
    const double* plyScale;     // per ply: parent arc length per unit u
    MeshView outward;           // library mm (interpolated vertex normals)
    const MeshView* hiddenMeshes;  // device array of nHidden meshes (common frame)
    const IslandWeight* islands;   // per control, at CurveRecord::ctrlOff of cs; null: none
};

} // namespace knit::gpu
