cmake_minimum_required(VERSION 3.18)
project(MyCudaOnlyProject LANGUAGES CUDA)

# Set CUDA standard and runtime preferences
set(CMAKE_CUDA_STANDARD 20)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_ARCHITECTURES 89)  # For RTX 4060m (Ada, SM_89)

# Optional: runtime library for MSVC to avoid /MT vs /MD mismatch
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreadedDLL")  # Use /MD

# Optional: If you need cuBLAS, cuRAND, etc.
find_package(CUDAToolkit REQUIRED)

# Find all CUDA source files
file(GLOB_RECURSE SOURCES src/*.cu)

# Create executable
add_executable(main ${SOURCES})

target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/include)


target_link_libraries(main PRIVATE CUDA::cublas)

# Optional: If you use separate compilation of device code
set_target_properties(main PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)
