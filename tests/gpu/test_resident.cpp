// Resident meshes leave the device when the job queue moves to another scene (gpu::releaseResident), and stay
// while the jobs stay on one scene. No device needed: a plain Resident stands for a mesh.
#include <chrono>
#include <memory>
#include <string>
#include <thread>

#include "gpu_common.h"
#include "gpu/device/runtime.h"
#include "jobs/scene_builder.h"

using namespace knit;

namespace {

void runKnitJob(SceneBuilder& b, const std::string& scene, const std::string& importSource) {
    auto s = std::make_shared<SceneSettings>();
    s->scene = scene;
    s->importSource = importSource;
    const int published = b.snapshot().published;
    b.startKnit(KnitRequest{s});
    const auto t0 = std::chrono::steady_clock::now();
    while ((b.snapshot().busy || b.snapshot().published == published) &&
           std::chrono::steady_clock::now() - t0 < std::chrono::seconds(20))
        std::this_thread::sleep_for(std::chrono::milliseconds(2));
    CHECK(b.snapshot().published == published + 1);
}

void putResident(uint64_t key) {
    auto r = std::make_shared<gpu::Resident>();
    r->bytes = 64;
    gpu::residentPut(key, r);
}

} // namespace

TEST(jobs_release_resident_meshes_on_another_scene) {
    SceneBuilder b([](const JobRequest&, const JobControl&, std::string*) { return std::make_shared<BuiltScene>(); });
    runKnitJob(b, "teddy", "");
    putResident(1);
    CHECK(gpu::memStats().resident == 64);
    runKnitJob(b, "teddy", "");  // the same scene: its meshes stay
    CHECK(gpu::memStats().resident == 64);
    CHECK(gpu::residentGet(1) != nullptr);
    runKnitJob(b, "wool-swatch", "");  // another scene: they go
    CHECK(gpu::memStats().resident == 0);
    CHECK(gpu::residentGet(1) == nullptr);
    putResident(2);
    runKnitJob(b, "wool-swatch", "");
    CHECK(gpu::memStats().resident == 64);
    runKnitJob(b, "import", "a.obj");  // two imported files are two scenes
    putResident(3);
    runKnitJob(b, "import", "b.obj");
    CHECK(gpu::memStats().resident == 0);
}
