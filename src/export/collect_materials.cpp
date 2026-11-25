// PartMaterials: every material of a scene, and which one each part takes; see collect_scene.h.
#include "export/collect_scene.h"
#include "model/log.h"

#include <array>
#include <cctype>
#include <cstdio>

#include "scene/run_log.h"

namespace knit {

std::string materialIdOf(std::string id) {
    for (char& c : id)
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') c = '_';
    return id;
}

namespace {

// A material id with a colour in it: "<prefix><rrggbb>".
std::string colourId(const char* prefix, uint32_t key) {
    char b[24];
    std::snprintf(b, sizeof b, "%s%06x", prefix, key);
    return std::string(b);
}

std::string partYarnMaterialId(const ScenePart& r) { return materialIdOf("yarn_" + r.name); }

} // namespace

std::string stitchColourId(uint32_t key) { return colourId("yarn_c", key); }

PartMaterials::PartMaterials(const Scene& t, const SceneSettings& s, ExportScene& scene, std::string& log)
    : t_(t), s_(s), scene_(scene), log_(log) {}

// A flat colour with no sheen: the table entry of an unknitted surface.
std::string PartMaterials::addPlainColour(const std::string& id, const Vec3& colour) {
    MaterialSpec m;
    m.color = colour;
    m.sheen = 0;
    scene_.baseMaterials.push_back({id, m});
    return id;
}

std::string PartMaterials::yarnOf(int part) const {
    const auto it = yarn_.find(part);
    return it != yarn_.end() ? it->second : std::string("wool_white");
}

// One colour per knitted part (material.yarn_colour_source = object): yarn entries after the two
// legacy ones, so they get the textured and fibre variants too. Parts of exactly the same colour
// share one entry (lossless): named after their source material when its colour is uniform (the
// teddy's nine pieces of a bear: one fabric), else after the first such part.
void PartMaterials::addObjectYarnMaterials() {
    std::map<std::array<double, 3>, std::string> idOfColour;
    std::set<std::string> ids;
    for (size_t i = 0; i < t_.parts.size(); ++i) {
        const ScenePart& r = t_.parts[i];
        if (!r.hasStitches()) continue;
        const Vec3 c = partYarnColour(t_, r, s_);
        const std::array<double, 3> key{c.x, c.y, c.z};
        auto it = idOfColour.find(key);
        if (it == idOfColour.end()) {
            const SourcePrimitive& g = t_.source.prims[r.prim];
            std::string id = g.triColor.empty() ? materialIdOf("yarn_" + g.material) : partYarnMaterialId(r);
            if (ids.count(id)) id = partYarnMaterialId(r);
            ids.insert(id);
            MaterialSpec m = s_.yarnWhite;
            m.color = c;
            scene_.baseMaterials.push_back({id, m});
            it = idOfColour.emplace(key, id).first;
        }
        yarn_[static_cast<int>(i)] = it->second;
    }
}

// Texture colourwork (material.yarn_colour_source = texture): one yarn entry per palette colour of
// the scene (stitch_colours.h); a part's dominant colour stands in for it where one colour is needed
// (the material id of its stray fibres). Without a texture on any knitted part: one colour per part.
void PartMaterials::addStitchColourMaterials() {
    const StageClock clock;
    SceneStitchColours colours = sceneStitchColours(t_);
    if (colours.palette.empty()) {
        appendFormat(log_, "colourwork: no knitted part has a base colour texture: one colour per part\n");
        addObjectYarnMaterials();
        return;
    }
    for (uint32_t k : colours.palette) {
        colourIndex[k] = static_cast<int>(colourKeys.size());
        colourKeys.push_back(k);
        MaterialSpec m = s_.yarnWhite;
        m.color = stitchColourLinear(k);
        scene_.baseMaterials.push_back({stitchColourId(k), m});
    }
    long stitches = 0;
    for (size_t i = 0; i < t_.parts.size(); ++i) {
        StitchColours& c = colours.parts[i];
        if (c.key.empty()) continue;
        std::map<uint32_t, int> n;
        for (uint32_t k : c.key) ++n[k];
        uint32_t dom = c.key.front();
        for (const std::pair<const uint32_t, int>& colour : n)
            if (colour.second > n[dom]) dom = colour.first;
        yarn_[static_cast<int>(i)] = stitchColourId(dom);
        dominant[static_cast<int>(i)] = dom;
        stitches += static_cast<long>(c.key.size());
        stitchColours[static_cast<int>(i)] = std::move(c);
    }
    appendFormat(log_, "colourwork: %zu yarn colours over %ld stitches (%d off the surface: their part's mean colour); %s\n", colourKeys.size(),
            stitches, colours.stitchesOffSurface, clock.line("colours").c_str());
}

void PartMaterials::addPlainMaterials() {
    for (const SourcePrimitive& g : t_.source.prims) {
        if (!plain.count(g.material)) plain[g.material] = addPlainColour(materialIdOf("plain_" + g.material), g.baseColor);
    }
}

void PartMaterials::addFurMaterials() {
    if (s_.pileColourSource == "single") {
        scene_.baseMaterials.push_back({"fur", s_.pileMaterial});
        return;
    }
    // Per part: the source "object", and with "texture" the parts without a texture.
    SceneSettings obj = s_;
    obj.yarnColourSource = "object";
    for (size_t i = 0; i < t_.parts.size(); ++i) {
        MaterialSpec m = s_.pileMaterial;
        m.color = partYarnColour(t_, t_.parts[i], obj);
        fur[static_cast<int>(i)] = "fur_" + partYarnMaterialId(t_.parts[i]).substr(5);  // "yarn_<part>" -> "fur_<part>"
        scene_.baseMaterials.push_back({fur[static_cast<int>(i)], m});
    }
}

void PartMaterials::addSceneMaterials() {
    // Knit skipped (knit.enabled = false): no yarn or wool fibre materials (nothing uses them).
    // ExportScene starts with the legacy table: replaced, or emptied.
    if (knitWorkflow(s_)) scene_.baseMaterials = {{"wool_white", s_.yarnWhite}, {"wool_brown", s_.yarnBrown}};
    else scene_.baseMaterials.clear();
    if (knitWorkflow(s_) && s_.yarnColourSource == "object") addObjectYarnMaterials();
    if (knitWorkflow(s_) && s_.yarnColourSource == "texture") addStitchColourMaterials();
    scene_.yarnEntries = static_cast<int>(scene_.baseMaterials.size());
    if (knitWorkflow(s_) && s_.filling) scene_.baseMaterials.push_back({"stuffing", s_.stuffing});
    addPlainMaterials();
    if (furWorkflow(s_)) addFurMaterials();
    scene_.fibreLook.hairMaterials = s_.fibreDirected.model == MaterialSpec::Model::Fibre;
    scene_.fibreLook.hair = s_.fibreDirected;
    scene_.fibreLook.random = s_.fibreRandom;
    if (knitWorkflow(s_) && !s_.fibreInheritsYarn) scene_.fibreLook.orientedOverride = s_.fibreDirected;
    scene_.fibreLook.randomInheritsYarn = knitWorkflow(s_) && s_.fibreRandomInheritsYarn;
    scene_.fibreLook.hairKeepsYarnColour = s_.fibreHairColour == "yarn";
}

std::string PartMaterials::texturedPlainMaterial(const SourcePrimitive& g, const std::string& texturePath) {
    const std::string key = plain[g.material] + "|" + texturePath;
    const auto known = texMaterial_.find(key);
    if (known != texMaterial_.end()) return known->second;
    // The texture's multiplier as triColor applies it: glTF's base colour factor; an
    // import's texture alone (its baseColor is the mean colour, texture included).
    const std::string id = "tex" + std::to_string(texMaterial_.size()) + "_" + plain[g.material].substr(6);  // "plain_<material>"
    addPlainColour(id, t_.source.imported ? Vec3{1, 1, 1} : g.baseColor);
    scene_.albedoTextures[id] = texturePath;
    texMaterial_.emplace(key, id);
    return id;
}

std::string PartMaterials::meanColourMaterial(const std::string& name, const Vec3& colour) {
    return addPlainColour(materialIdOf("mean_" + name), colour);
}

std::string PartMaterials::furColourMaterial(uint32_t colour) {
    const std::string id = colourId("fur_c", colour);
    if (furColours_.insert(colour).second) {
        MaterialSpec m = s_.pileMaterial;
        m.color = stitchColourLinear(colour);
        scene_.baseMaterials.push_back({id, m});
    }
    return id;
}

} // namespace knit
