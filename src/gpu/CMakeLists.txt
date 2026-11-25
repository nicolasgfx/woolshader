# The CUDA backend (ARCHITECTURE.md): the kernels (single precision) and the host
# adapters. The scenes compute on it; CUDA is required. Public: the four headers the program calls
# (device, knitting port, surfaces, fibre generator); everything else, CUDA included, is private.
set(KNIT_CUDA_ARCH 89 CACHE STRING "CUDA architecture of the GPU backend (89: RTX 40xx)")
knit_component(gpu
    SOURCES
        device/primitives.cu device/residency.cpp device/runtime.cu
        fibres/fibres.cu
        geometry/arc_length.cu geometry/device_geometry.cu
        port/knit_numerics.cpp port/ncp_call_file.cpp
        solver/ncp_control.cu solver/ncp_evaluate.cu solver/ncp_live.cu solver/ncp_pcg.cu solver/ncp_refresh.cu
        solver/ncp_setup.cu solver/ncp_step.cu
        surface/lattice.cu surface/mesh_surface.cu surface/piece_surface.cu surface/plies.cu
        surface/scene_queries.cu surface/surface_build.cu
        validation/clearance.cu validation/curvature.cu validation/signature.cu
    PUBLIC_HEADERS
        device.h fibre_generator.h knitting.h surfaces.h
    PRIVATE_HEADERS
        device/kernel.cuh device/primitives.h device/runtime.h
        fibres/fibre_candidate.cuh fibres/fibres.h
        geometry/device_geometry.h geometry/math.cuh geometry/surface.cuh
        port/ncp_call_file.h
        solver/ncp.h solver/ncp_device.cuh
        surface/placement.h
        validation/checks.h
    PUBLIC_DEPS model knitting fibres
    PRIVATE_LIBS CUDA::cudart_static
)
set_target_properties(knit_gpu PROPERTIES CUDA_STANDARD 20 CUDA_STANDARD_REQUIRED ON CUDA_ARCHITECTURES "${KNIT_CUDA_ARCH}")
# No fast math (the first comparison is made without it); IEEE division / square root.
target_compile_options(knit_gpu PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo -Xcompiler=/utf-8 -Xcompiler=/fp:precise -Xcompiler=/Zc:preprocessor
                       --diag-suppress=177>)
