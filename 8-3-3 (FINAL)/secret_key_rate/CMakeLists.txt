cmake_minimum_required(VERSION 3.24)

# Silence the FindBoost deprecation warning
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)
endif()

project(secret_key_rate LANGUAGES CXX)

# Standard C++20 for MSVC
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- 1. THE FIX: Modern Boost Discovery ---
# We look for the specific component config file provided by vcpkg.
find_package(boost_program_options CONFIG REQUIRED)

# --- 2. THE OPTIMIZATION: IPO / Whole Program Optimization ---
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT error)

if(ipo_supported)
    message(STATUS "IPO/LTO (Whole Program Optimization) enabled.")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
else()
    message(STATUS "IPO not supported: ${error}")
endif()

# --- 3. Source Collection ---
file(GLOB SOURCES "src/*.cpp")
set(EXTERNAL_SOURCES "../effective_error_probability/src/AuxMath.cpp")

add_executable(main ${SOURCES} ${EXTERNAL_SOURCES})

target_include_directories(main PRIVATE 
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "../effective_error_probability/include"
)

# --- 4. Linking ---
# vcpkg provides the target as Boost::program_options or boost_program_options::boost_program_options
target_link_libraries(main PRIVATE Boost::program_options)

# --- 5. Compiler Specific Performance Tweaks ---
if(MSVC)
    # /O2: Speed, /Oi: Intrinsic functions, /Ot: Favor fast code, /GL: Whole program optimization
    target_compile_options(main PRIVATE /O2 /Oi /Ot /W4)
    # Enable AVX2 if your processor supports it (standard for 2026 hardware)
    target_compile_options(main PRIVATE /arch:AVX2)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(main PRIVATE -O3 -march=native -flto)
    target_link_options(main PRIVATE -O3 -flto)
endif()