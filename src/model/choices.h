#pragma once
// Choices a settings document holds (settings/scene_settings.h) and a component below settings carries out:
//   HandKnit         knitting, stitch/hand_tension: how unevenly a hand knits
//   FibreIslands     fibres, fibre_islands: patches of different fuzz
//   PlyProfileModel  knitting, tube/ply_profile: the cross-section of a ply
// They live here and not with that component, because the viewer reads settings' headers and may see neither
// knitting nor fibres (src/ARCHITECTURE.md).

namespace knit {

// How unevenly a hand knits (knit.hand.*). All zero: the template, exactly.
struct HandKnit {
    int seed = 1;
    double unevenness = 0;  // [% of the stitch] how far a stitch's tension strays
    double columnShift = 0; // [mm] how far a column of stitches slides against its neighbours
    double columnRows = 6;  // [rows] how many rows a column keeps about the same slide
    double loopRate = 0;    // raised loops per cm2
    double loopHeight = 0;  // [mm] the tallest raised loop
    double waveHeight = 0;  // [mm] slow rise and fall of the fabric
    double spotRate = 0;    // thick or thin spots per cm2
    double spotChange = 0;  // [% of D] the largest change of diameter
    double spotLength = 12; // [mm] along the yarn
    bool inUse() const { return unevenness > 0 || columnShift > 0 || loopRate > 0 || waveHeight > 0 || spotRate > 0; }
    bool operator==(const HandKnit&) const = default;
};

// Patches of different fibres on the even fuzz (fibres.islands.*). rate 0: none.
struct FibreIslands {
    double rate = 0;      // islands per cm2
    double size = 12;     // [mm] across
    double strength = 0;  // [%] how far an island's centre departs from the fuzz around it
    double fuzzy = 1, matted = 1, worn = 1;  // relative shares of the three kinds
    bool inUse() const { return rate > 0 && strength > 0; }
    bool operator==(const FibreIslands&) const = default;
};

// The cross-section of one ply in the shading profile (knitting/tube/ply_profile.h).
enum class PlyProfileModel { Circles, Swept };

} // namespace knit
