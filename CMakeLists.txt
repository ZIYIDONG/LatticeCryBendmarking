cmake_minimum_required(VERSION 3.22)
project(LatticeCryBenchmarking)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2")

set(INC_DIR ${CMAKE_SOURCE_DIR}/include_plain-LWE)

# ── 参数集选择 ──
#    cmake .. -DUSE_128BIT_PARAMS=ON   → 128-bit 安全参数 (n=512, q=134219777)
#    cmake ..                           → 演示参数 (n=8, q=257)
option(USE_128BIT_PARAMS "Use 128-bit security parameters" OFF)
if(USE_128BIT_PARAMS)
    add_compile_definitions(LATTICE_128BIT)
    message(STATUS ">>> Building with 128-bit security parameters (n=512, q=134219777)")
endif()

# ─── 辅助宏 ───
function(add_demo name)
    add_executable(${name} ${ARGN})
    target_include_directories(${name} PRIVATE ${INC_DIR})
endfunction()

# Build all demo sources as a static library (exclude main and runner)
file(GLOB DEMO_SRCS
    ${CMAKE_SOURCE_DIR}/src_plain-LWE/*_plain-LWE.cpp
)
list(FILTER DEMO_SRCS EXCLUDE REGEX "main_plain-LWE\\.cpp$")
list(FILTER DEMO_SRCS EXCLUDE REGEX "runner_plain-LWE\\.cpp$")
add_library(demos STATIC ${DEMO_SRCS})
target_include_directories(demos PRIVATE ${INC_DIR})

# Single executable
add_executable(LatticeCryBenchmarking src_plain-LWE/main_plain-LWE.cpp)
target_include_directories(LatticeCryBenchmarking PRIVATE ${INC_DIR})
target_link_libraries(LatticeCryBenchmarking PRIVATE demos)
