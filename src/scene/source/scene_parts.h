#pragma once
// loadPartsToKnit and loadBaseParts: from a SceneEntry to its parts, with nothing knitted.
//
// Reads the GLB, the stored import or the sewn pieces, and makes one ScenePart per component or
// chart. Called by the knit stage, the saved-knit restore, the fibre job and knit_check. An empty
// result or false comes with the reason when the scene or its settings cannot be used.
#include <optional>
#include <string>

#include <memory>

#include "scene/scene.h"
#include "scene/source/guided_source.h"
#include "scene/source/sewn_pieces.h"

namespace knit {

// What knitting needs beside the parts.
struct KnitSource {
    StitchPattern pattern = StitchPattern::Stockinette;
    SpiralSource pieces;  // the spiral layout's pieces as the knitting places them; empty on UV charts
    std::shared_ptr<const GuidedInputs> guided;  // the guided route's file, markings and piece (knit.route = guided); null otherwise
};
// The scene's parts for knitting, by its layout (a knitted scene of sewn pieces as spirals,
// sewn_pieces.h; an import on the guided route as its one marked piece, guided_source.h; every other
// on its mesh), with the knitting settings checked.
std::optional<KnitSource> loadPartsToKnit(const SceneEntry& entry, const SceneSettings& settings, const SourceOptions& options, Scene& out,
                                          std::string* err);
// The meshes as parts, with neither knitting charts nor a valid knitting setup required (a mesh without
// uv loads): what fur grows on (fur_stage.h), and a scene whose Knit step is skipped.
bool loadBaseParts(const SceneEntry& entry, const SceneSettings& settings, const SourceOptions& options, Scene& out, std::string* err);

// Content hash of a scene's inputs: the GLB file, or the pieces manifest and every piece file.
uint64_t sceneAssetHash(const SceneEntry& entry, std::string* err = nullptr);

} // namespace knit
