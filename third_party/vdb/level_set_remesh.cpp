// The level set remesh through OpenVDB. The calls follow OpenVDB's documentation (the OpenVDB Cookbook; tools/MeshToVolume.h,
// tools/VolumeToMesh.h): meshToLevelSet with a linear transform of one voxel size and the default narrow band, then
// volumeToMesh at the zero level without adaptivity.
#include "level_set_remesh.h"

#include <exception>

#include <openvdb/openvdb.h>
#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>

namespace knit {

bool remeshThroughLevelSet(const LevelSetInput& in, LevelSetOutput& out, std::string* err) {
    out = LevelSetOutput{};
    try {
        openvdb::initialize();
        std::vector<openvdb::Vec3s> points(in.positions.size());
        for (size_t i = 0; i < in.positions.size(); ++i) points[i] = {in.positions[i][0], in.positions[i][1], in.positions[i][2]};
        std::vector<openvdb::Vec3I> triangles(in.triangles.size());
        for (size_t i = 0; i < in.triangles.size(); ++i) triangles[i] = {in.triangles[i][0], in.triangles[i][1], in.triangles[i][2]};

        const openvdb::math::Transform::Ptr transform = openvdb::math::Transform::createLinearTransform(in.voxelSize);
        const openvdb::FloatGrid::Ptr levelSet = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*transform, points, triangles);

        std::vector<openvdb::Vec3s> meshPoints;
        std::vector<openvdb::Vec3I> meshTriangles;
        std::vector<openvdb::Vec4I> meshQuads;
        const double isovalue = 0.0, adaptivity = 0.0;
        openvdb::tools::volumeToMesh(*levelSet, meshPoints, meshTriangles, meshQuads, isovalue, adaptivity);

        out.positions.resize(meshPoints.size());
        for (size_t i = 0; i < meshPoints.size(); ++i) out.positions[i] = {meshPoints[i].x(), meshPoints[i].y(), meshPoints[i].z()};
        out.quads.resize(meshQuads.size());
        for (size_t i = 0; i < meshQuads.size(); ++i) out.quads[i] = {meshQuads[i][0], meshQuads[i][1], meshQuads[i][2], meshQuads[i][3]};
        out.triangles.resize(meshTriangles.size());
        for (size_t i = 0; i < meshTriangles.size(); ++i) out.triangles[i] = {meshTriangles[i][0], meshTriangles[i][1], meshTriangles[i][2]};
        return true;
    } catch (const std::exception& e) {
        if (err) *err = std::string("the level set remesh failed: ") + e.what();
    } catch (...) {
        if (err) *err = "the level set remesh failed with an unknown error";
    }
    out = LevelSetOutput{};
    return false;
}

} // namespace knit
