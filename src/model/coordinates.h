#pragma once
// Where a point is. The library works in Z-up millimetres; a scene has its own units and centre
// (p_lib = sigma * (p_scene - centre)); the exported scene is Y-up; a source file keeps its own units and up axis.
// Import, export, the settings' scene mapping and the strand and piece writers use these; the solver does not.
#include "model/vec.h"

namespace knit {

// Library millimetres to scene units, and on to the exported scene.
struct SceneMapping {
    double sigma = 1;   // [mm per scene unit]
    Vec3 centre;        // [scene units]
    Vec3 toScene(const Vec3& q) const { return q / sigma + centre; }
    // A library point in the exported scene (Y-up), and the same rounded through float32 as the files hold it.
    Vec3 toExport(const Vec3& p) const {
        const Vec3 q = toScene(p);
        return {q.x, q.z, -q.y};
    }
    Vec3 toExportF(const Vec3& p) const {
        const Vec3 q = toExport(p);
        return {double(float(q.x)), double(float(q.y)), double(float(q.z))};
    }
};

// A source file's up axis (OBJ convention Y, or Z), and a library point (Z-up mm) back in the file's own
// units and axes: p_file = (p_lib - offset) / scale, then Z-up to the file's up axis.
enum class UpAxis { Y, Z };
inline Vec3 toFileCoordinates(double scale, const Vec3& offset, UpAxis up, const Vec3& lib) {
    const Vec3 z = (lib - offset) * (1.0 / scale);  // Z-up file units
    return up == UpAxis::Y ? Vec3{z.x, z.z, -z.y} : z;
}

} // namespace knit
