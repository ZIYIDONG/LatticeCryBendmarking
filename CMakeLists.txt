cmake_minimum_required(VERSION 3.22)
project(LatticeCryBenchmarking)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2")

set(INC_DIR ${CMAKE_SOURCE_DIR}/include_plain-LWE)

# ── 安全级别选择 ──
#    cmake .. -DSECURITY_LEVEL=1  → NIST Level 1 (AES-128 等效, n=512)
#    cmake .. -DSECURITY_LEVEL=3  → NIST Level 3 (AES-192 等效, n=768)
#    cmake .. -DSECURITY_LEVEL=5  → NIST Level 5 (AES-256 等效, n=1024)
#    cmake ..                     → Demo          (n=8,  仅算法验证)
set(SECURITY_LEVEL "" CACHE STRING "NIST security level: 1, 3, or 5 (empty=demo)")
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
    add_compile_definitions(LATTICE_DEMO)
    message(STATUS ">>> Building Demo (n=8, q=257)")
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
