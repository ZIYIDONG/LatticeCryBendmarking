cmake_minimum_required(VERSION 3.22)
project(LatticeCryBenchmarking)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2")

set(INC_DIR ${CMAKE_SOURCE_DIR}/include)

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
