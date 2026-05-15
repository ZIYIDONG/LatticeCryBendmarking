cmake_minimum_required(VERSION 3.22)
project(LatticeCryBenchmarking)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2")

set(INC_DIR ${CMAKE_SOURCE_DIR}/include_plain-LWE)

# ── 安全级别选择 ──
#    cmake ..                            → Demo (default, n=64, ~40-bit)
#    cmake .. -DSECURITY_LEVEL=1         → NIST L1 (n=512, ~128-bit)
#    cmake .. -DSECURITY_LEVEL=3         → NIST L3 (n=768, ~192-bit)
#    cmake .. -DSECURITY_LEVEL=5         → NIST L5 (n=1024,~256-bit)
set(SECURITY_LEVEL "" CACHE STRING "Security level: 1, 3, or 5 (empty=demo, n=64)")
if(SECURITY_LEVEL STREQUAL "5")
    add_compile_definitions(LATTICE_LEVEL_5)
    message(STATUS ">>> Building NIST Level 5 (n=1024, q~2^32)")
elseif(SECURITY_LEVEL STREQUAL "3")
    add_compile_definitions(LATTICE_LEVEL_3)
    message(STATUS ">>> Building NIST Level 3 (n=768, q~2^32)")
elseif(SECURITY_LEVEL STREQUAL "1")
    add_compile_definitions(LATTICE_LEVEL_1)
    message(STATUS ">>> Building NIST Level 1 (n=512, q~2^27)")
else()
    message(STATUS ">>> Building Demo (n=64, q=8191, ~40-bit)")
endif()

# ─── 辅助宏 ───
function(add_demo name)
    add_executable(${name} ${ARGN})
    target_include_directories(${name} PRIVATE ${INC_DIR})
endfunction()

# Build all demo sources as a static library
file(GLOB DEMO_SRCS
    ${CMAKE_SOURCE_DIR}/src_plain-LWE/*_plain-LWE.cpp
)
list(FILTER DEMO_SRCS EXCLUDE REGEX "main_plain-LWE\\.cpp$")
list(FILTER DEMO_SRCS EXCLUDE REGEX "runner_plain-LWE\\.cpp$")
add_library(demos STATIC ${DEMO_SRCS})
target_include_directories(demos PRIVATE ${INC_DIR})

add_executable(LatticeCryBenchmarking src_plain-LWE/main_plain-LWE.cpp)
target_include_directories(LatticeCryBenchmarking PRIVATE ${INC_DIR})
target_link_libraries(LatticeCryBenchmarking PRIVATE demos)
