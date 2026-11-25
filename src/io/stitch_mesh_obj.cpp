// A built stitch mesh as an OBJ file, and read back (stitch_mesh_obj.h).
#include "io/stitch_mesh_obj.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace knit {

namespace {

// stitchmesh_lab's names for the kinds (kindGroup in cli/stitchmesh_lab.cpp), in StitchKind's order.
const char* const kKindName[6] = {"K", "INC_kp", "DEC_d12k", "special_quad", "unresolved", "pole_start"};

// The row and kind a group name states; false for a group that is not a row group.
bool parseGroupName(const std::string& name, int& row, StitchKind& kind) {
    if (name.rfind("row_", 0) != 0) return false;
    const size_t underscore = name.find('_', 4);
    if (underscore == std::string::npos) return false;
    const std::string number = name.substr(4, underscore - 4), kindName = name.substr(underscore + 1);
    row = number == "none" ? -1 : std::atoi(number.c_str());
    for (int k = 0; k < 6; ++k)
        if (kindName == kKindName[k]) {
            kind = static_cast<StitchKind>(k);
            return true;
        }
    return false;
}

} // namespace

std::string stitchGroupName(int row, StitchKind kind) {
    char number[16];
    if (row < 0) std::snprintf(number, sizeof number, "none");
    else std::snprintf(number, sizeof number, "%03d", row);
    return std::string("row_") + number + "_" + kKindName[static_cast<int>(kind)];
}

bool writeStitchMeshObj(const std::string& path, const StitchMeshData& mesh, const std::vector<std::string>& comments, const std::string& objectName, std::string* err) {
    const StitchMeshUse use = stitchMeshUse(mesh);
    if (!use.drawable) {
        if (err) *err = use.reason;
        return false;
    }
    std::ofstream f(path, std::ios::binary);
    if (!f) {
        if (err) *err = "cannot write " + path;
        return false;
    }
    for (const std::string& c : comments) f << "# " << c << "\n";
    f << "o " << (objectName.empty() ? std::string("stitch_mesh") : objectName) << "\n";
    char line[128];
    for (const Vec3& p : mesh.pos) {
        std::snprintf(line, sizeof line, "v %.17g %.17g %.17g\n", p.x, p.y, p.z);
        f << line;
    }
    std::string group;
    for (int face = 0; face < mesh.faces(); ++face) {
        const int row = face < static_cast<int>(mesh.row.size()) ? mesh.row[face] : -1;
        const StitchKind kind = face < static_cast<int>(mesh.kind.size()) ? mesh.kind[face] : StitchKind::Regular;
        const std::string name = stitchGroupName(row, kind);
        if (name != group) {
            group = name;
            f << "g " << group << "\n";
        }
        f << "f";
        for (int k = mesh.faceStart[face]; k < mesh.faceStart[face + 1]; ++k) f << " " << mesh.corner[k] + 1;
        f << "\n";
    }
    f.close();
    if (!f) {
        if (err) *err = "writing " + path + " failed";
        return false;
    }
    return true;
}

bool readPolygonObj(const std::string& path, StitchMeshData& out, std::string* err) {
    out = StitchMeshData{};
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        if (err) *err = "cannot read " + path;
        return false;
    }
    int row = -1;
    StitchKind kind = StitchKind::Regular;
    std::string line;
    int lineNumber = 0;
    out.faceStart.push_back(0);
    while (std::getline(f, line)) {
        ++lineNumber;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream in(line);
        std::string tag;
        in >> tag;
        if (tag == "v") {
            Vec3 p;
            if (!(in >> p.x >> p.y >> p.z)) {
                if (err) *err = path + ":" + std::to_string(lineNumber) + ": a vertex needs three numbers";
                return false;
            }
            out.pos.push_back(p);
        } else if (tag == "g") {
            std::string name;
            in >> name;
            row = -1;
            kind = StitchKind::Regular;
            parseGroupName(name, row, kind);
        } else if (tag == "f") {
            std::string corner;
            while (in >> corner) {
                const int index = std::atoi(corner.c_str());   // the vertex index before any '/'
                out.corner.push_back(index < 0 ? static_cast<int>(out.pos.size()) + index : index - 1);
            }
            out.faceStart.push_back(static_cast<int>(out.corner.size()));
            out.kind.push_back(kind);
            out.row.push_back(row);
            out.finding.push_back(0);
            out.deviationMm.push_back(0);
        }
    }
    if (out.faceStart.size() == 1) out.faceStart.clear();
    return true;
}

} // namespace knit
