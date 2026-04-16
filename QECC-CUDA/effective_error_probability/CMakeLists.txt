cmake_minimum_required(VERSION 3.21)
project(Eff_Err_Prob LANGUAGES CXX CUDA)

# # Force NVCC to accept newer MSVC (trial only, may be unstable)
# set(CMAKE_CUDA_HOST_COMPILER "cl.exe")
# set(CMAKE_CUDA_FLAGS_INIT "--allow-unsupported-compiler")


# ------------------------
# Language & Standards
# ------------------------
# Force C++20 for all targets
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CUDA standard
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# ------------------------
# CUDA Architectures
# RTX 4060m -> Ada (compute capability 8.9)
# Optionally add "native" for forward-compatibility
# ------------------------
set(CMAKE_CUDA_ARCHITECTURES 89)

# ------------------------
# MSVC runtime (dynamic /MD)
# ------------------------
if (MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")
endif()

# ------------------------
# CUDA Toolkit (cuBLAS, cuRAND, etc.)
# ------------------------
find_package(CUDAToolkit REQUIRED)

# ------------------------
# Sources
# ------------------------
file(GLOB_RECURSE SOURCES src/*.cpp src/*.cu)

add_executable(main ${SOURCES})
target_compile_options(main PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--allow-unsupported-compiler>
)
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/include)

# ------------------------
# Linking CUDA libraries (extend if needed)
# ------------------------
target_link_libraries(main PRIVATE CUDA::cublas)

# ------------------------
# Properties for CUDA compilation
# ------------------------
set_target_properties(main PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin
)

# ------------------------
# Debug configuration for kernels
# (-G = generate debug info)
# ------------------------
target_compile_options(main PRIVATE
    $<$<CONFIG:Debug>:-G>
)

# ------------------------
# Additional nice-to-haves
# ------------------------

# # Let NVCC be more permissive with MSVC versions (avoids "unsupported compiler" warnings)
# set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --allow-unsupported-compiler")

# Treat .cuh as header files for IDEs
set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CUDA)