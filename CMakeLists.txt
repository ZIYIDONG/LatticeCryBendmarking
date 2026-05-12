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

# ─── 主 MP12 程序 (包含 Test 1–11) ───
add_demo(LatticeCryBenchmarking
    src/main.cpp
    src/mp12deltrapgen.cpp
)

# ─── Powersof 系列 ───
add_demo(test_powersof           src/test_powersof.cpp)
add_demo(test_powersof_modswitch src/test_powersof_modswitch.cpp)

# ─── FRD 系列 ───
add_demo(test_frd  src/test_frd.cpp)
add_demo(debug_frd src/debug_frd.cpp)

# ─── matops 基准测试 ───
add_demo(bench_matops src/bench_matops.cpp)

# ─── expand/extend 基准测试 ───
add_demo(test_expand src/test_expand.cpp)

# ─── eval 测试 ───
add_demo(test_eval src/test_eval.cpp)

# ─── 多身份 FHE 阈值解密 ───
add_demo(test_decrypt  src/test_decrypt.cpp)
add_demo(bench_decrypt src/bench_decrypt.cpp)
