cmake_minimum_required(VERSION 3.22)
project(LatticeCryBenchmarking)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2")

set(INC_DIR ${CMAKE_SOURCE_DIR}/include)

# Option: unified parameters
option(UNIFIED_PARAMS "Use unified parameters (generate params_cfg.h from CMake defaults)" ON)
# Default unified parameter values (can be overridden via -D on cmake command line)
# Standard 128-bit parameters (commented out)
# set(MP12_N   512 CACHE STRING "MP12 n")
# set(MP12_Q   134219777 CACHE STRING "MP12 q")
# set(MP12_B   27 CACHE STRING "MP12 b")
# set(MID_LAMBDA 128 CACHE STRING "MID lambda")
# set(MID_N_ID   3 CACHE STRING "MID N_id")
# set(MID_D      1 CACHE STRING "MID d")
# set(MP12_SIGMA 1.0 CACHE STRING "MP12 sigma")
# set(UNIENC_SIGMA 3.2 CACHE STRING "UniEnc LWE sigma")

# Small/demo parameters (default)
set(MP12_N   8 CACHE STRING "MP12 n")
set(MP12_Q   257 CACHE STRING "MP12 q")
set(MP12_B   2 CACHE STRING "MP12 b")
set(MP12_SIGMA 1.0 CACHE STRING "MP12 sigma")
set(UNIENC_SIGMA 3.2 CACHE STRING "UniEnc LWE sigma")
set(MID_LAMBDA 8 CACHE STRING "MID lambda")
set(MID_N_ID   3 CACHE STRING "MID N_id")
set(MID_D      1 CACHE STRING "MID d")

# Configure params header into build directory when enabled
if(UNIFIED_PARAMS)
    set(UNIFIED_PARAMS_ENABLED 1)
    configure_file(${CMAKE_SOURCE_DIR}/include/params_cfg.h.in ${CMAKE_BINARY_DIR}/generated/params_cfg.h @ONLY)
    include_directories(${CMAKE_BINARY_DIR}/generated)
else()
    set(UNIFIED_PARAMS_ENABLED 0)
endif()

# ─── 辅助宏: 少写重复代码 ───
function(add_demo name)
    add_executable(${name} ${ARGN})
    target_include_directories(${name} PRIVATE ${INC_DIR})
endfunction()

# Build all demo sources as a static library (exclude main.cpp and runner.cpp)
file(GLOB DEMO_SRCS
    ${CMAKE_SOURCE_DIR}/src/*.cpp
)
list(REMOVE_ITEM DEMO_SRCS ${CMAKE_SOURCE_DIR}/src/main.cpp ${CMAKE_SOURCE_DIR}/src/runner.cpp)
add_library(demos STATIC ${DEMO_SRCS})
target_include_directories(demos PRIVATE ${INC_DIR} ${CMAKE_BINARY_DIR}/generated)

# Single executable: main.cpp is the sole entry point and links the demos library
add_executable(LatticeCryBenchmarking src/main.cpp)
target_include_directories(LatticeCryBenchmarking PRIVATE ${INC_DIR} ${CMAKE_BINARY_DIR}/generated)
target_link_libraries(LatticeCryBenchmarking PRIVATE demos)

# Note: runner.cpp is kept in the tree but not built. If you prefer the runner as the entry,
# change the executable source above to src/runner.cpp and adjust accordingly.

