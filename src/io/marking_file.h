#pragma once
// The surface markings file, "knit-markings v5" (surface_marking.h; v1 to v4 files are read and rewritten as v5): plain text, one statement per line, so a set
// can be read, diffed and written by hand or by another program (a proposal written by the stitch-mesh lab). Two
// users: the viewer's Mark step and stitchmesh_lab. Unknown keys refuse the file; barycentric weights round-trip bit
// for bit (%.17g).
//
//   knit-markings v5                   (v4: without knit_component; v3: also without setup and basis; v2: also without the preparation, meaning facets;
//                                      v1: also without settings)
//   source = scenes/teddy-bears/obj/bears.obj
//   digest = 0123456789abcdef          meshDigest of the file's mesh as read
//   next_id = 4
//   preparation = facets               v3: facets | loop_limit, the reference the mesher prepares (prepared_reference.h)
//   marking 3                          a block per marking, keys in any order, defaults when absent
//   meaning = row_start                start | finish | row_start | boundary | cut | attachment
//   shape = path                       point | path | closed
//   origin = user                      user | proposal
//   proposal = none                    none | offered | accepted | rejected
//   proposal_id = 0
//   group = 0
//   setup = 2                          v4, written when > 0: proposals offered together as one compatible setup
//   basis = 0123456789abcdef           v4, written when non-zero: what a proposal was computed for (proposalBasis)
//   knit_component = 1                 v5, written when > 0: the knitting component of its piece (0: the first)
//   note = inner side                  free text to the end of the line (note, reason, unsuitable)
//   setting ring1 = auto               v2: a construction setting of this marking (key without spaces)
//   point 61234 0.25 0.25 0.5          file triangle and barycentric weights, in order
//   end
#include <string>

#include "model/surface_marking.h"

namespace knit {

std::string writeMarkings(const MarkingSet& s);
[[nodiscard]] bool readMarkings(const std::string& text, MarkingSet& out, std::string* err);

// The same through a file (written beside the target and renamed over it).
[[nodiscard]] bool saveMarkingFile(const std::string& path, const MarkingSet& s, std::string* err);
[[nodiscard]] bool loadMarkingFile(const std::string& path, MarkingSet& out, std::string* err);

} // namespace knit
