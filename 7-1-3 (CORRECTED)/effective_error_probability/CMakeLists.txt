cmake_minimum_required(VERSION 3.24)

# Set build type to Release by default if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

project(generate_reencoding_errors LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optimization: Enable Interprocedural Optimization (LTO/IPO)
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)
if(ipo_supported)
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
    message(STATUS "LTO/IPO is enabled.")
endif()

# Find OpenMP
find_package(OpenMP REQUIRED)

# Collect sources (Explicitly listing is better, but keeping GLOB for your workflow)
file(GLOB_RECURSE SOURCES "src/*.cpp")
file(GLOB_RECURSE HEADERS "include/*.hpp")

add_executable(main ${SOURCES} ${HEADERS})

# Link OpenMP and set include directories
target_link_libraries(main PRIVATE OpenMP::OpenMP_CXX)
target_include_directories(main PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Optimization: Add architecture-specific flags for GCC
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    # -O3 is handled by CMAKE_BUILD_TYPE=Release automatically
    # -march=native tells GCC to optimize for YOUR specific CPU
    target_compile_options(main PRIVATE -march=native)
endif()