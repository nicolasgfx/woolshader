#pragma once
// What the two source files of stitchmesh_lab share: the arm commands (stitchmesh_arm.cpp) and the typed-mesh writers
// and readers of the earlier checkpoints (stitchmesh_lab.cpp), so both write the same files.
#include <map>
#include <string>
#include <vector>

#include "model/mesh.h"

#include "stitchmesh/relax_mesh.h"
#include "stitchmesh/typed_mesh.h"
#include "stitchmesh/virtual_shift.h"
#include "stitchmesh/yarn_embed.h"

// The imported-arm commands: import, propose, markings, arm. Returns the exit code, 2 for an unknown command.
int armCommand(const std::string& cmd, int argc, char** argv);
// The smooth-reference commands (stitchmesh_smooth.cpp): prep-compare, smooth-fixture, shaping-fixtures. 2 for an unknown command.
int smoothCommand(const std::string& cmd, int argc, char** argv);
// The yarn-generation commands (stitchmesh_yarn.cpp): yarn-fixtures, yarn-embed. 2 for an unknown command.
int yarnCommand(const std::string& cmd, int argc, char** argv);
// The yarn-settle commands (stitchmesh_settle.cpp): yarn-settle, yarn-settle-fixtures. 2 for an unknown command.
int settleCommand(const std::string& cmd, int argc, char** argv);

namespace knit::lab {

// The typed mesh's files as the build command writes them (candidate.obj, face_types.csv, rows.csv, ...).
void writeTypedMeshFiles(const std::string& dir, const TypedStitchMesh& s, const TypedMeshCheck& c, const VirtualShift* vs, const VirtualShiftCheck* vc);
// Reads them back (the reload): false when a file is missing or malformed.
bool readTypedMeshFiles(const std::string& dir, TypedStitchMesh& s, std::string& savedDigest);
void writeCellsFile(const std::string& path, const std::vector<CellMeasure>& cells);
void writePolygonObj(const std::string& path, const TypedStitchMesh& s);
void writeRelaxHistory(const std::string& path, const RelaxProduct& p);
std::string relaxSummary(const RelaxProduct& p);
std::string relaxSettings(const MeshRelaxRequest& r);


// The yarn commands' pieces (stitchmesh_yarn.cpp) the yarn-settle commands reuse, so both build and write the same cases.
struct YarnFixture {
    std::string name, purpose;
    bool pole = false;
    std::vector<int> counts;
    std::vector<std::vector<int>> shaping;
    bool boundaryShaping = false;
    bool negative = false;
    int flipHeadsOfSlot = -1, skipSlot = -1;   // negative controls: by (strip 1, slot)
};
std::vector<YarnFixture> yarnFixtureCases();
// Builds, relaxes (CP3 settings) and returns the fixture's typed mesh, virtual faces and reference.
bool buildYarnFixture(const YarnFixture& fx, double w, double h, TypedStitchMesh& m, std::vector<std::vector<int>>& faces, TriMesh& ref, RelaxProduct& relax, std::string* err);
YarnSettings yarnSettingsFrom(const std::map<std::string, std::string>& kv);   // --diameter, --join-pole, --binding, --tangents
// Centrelines, pieces, controls, tubes (radius D / 2, 8 sides), the check and its per-face, contact and summary tables.
void writeYarnFiles(const std::string& dir, const TypedStitchMesh& s, const YarnEmbedding& e, const YarnCheck& c, const YarnSettings& st, const std::string& title);
bool readReferenceObj(const std::string& path, TriMesh& ref);   // v / f lines; vertex normals computed
void writeTriangleObj(const std::string& path, const TriMesh& m);

} // namespace knit::lab
