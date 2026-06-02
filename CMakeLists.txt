cmake_minimum_required(VERSION 3.22)
project(LatticeCryBenchmarking)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2")

set(INC_DIR ${CMAKE_SOURCE_DIR}/include_plain-LWE)
set(INC_RLWEOR_NTRU_DIR ${CMAKE_SOURCE_DIR}/include_RLWEorNTRU)

# ── 安全级别选择 (NIST IR 8413) ──
#
#   统一 CMake 标志 -DSECURITY_LEVEL=N，两模块各自映射到 NIST 级别:
#
#   标志           plain-LWE (KEM/加密)       IBAGS (签名)
#   ────────────  ──────────────────────    ──────────────────
#   Demo (默认)    n=64   Demo              n=64   Demo
#   SEC_LEVEL=1    n=512  NIST Level 1      n=512  NIST Level 2  (签名无 Level 1)
#   SEC_LEVEL=2    (同 Level 1)             n=512  NIST Level 2
#   SEC_LEVEL=3    n=768  NIST Level 3      n=1024 NIST Level 3  (Dilithium3)
#   SEC_LEVEL=4    (同 Level 3)             n=1024 NIST Level 4* 
#   SEC_LEVEL=5    n=1024 NIST Level 5      n=1024 NIST Level 5  (Dilithium5)
#
#   参考: Kyber (L1/L3/L5), Dilithium (L2/L3/L5)
set(SECURITY_LEVEL "" CACHE STRING "NIST security level: 1-5 (empty=demo)")
if(SECURITY_LEVEL STREQUAL "5")
    add_compile_definitions(LATTICE_LEVEL_5)
    message(STATUS ">>> NIST Level 5 — plain-LWE: n=1024 | IBAGS: n=1024 (Dilithium5)")
elseif(SECURITY_LEVEL STREQUAL "4")
    add_compile_definitions(LATTICE_LEVEL_5)
    message(STATUS ">>> NIST Level 4 — plain-LWE: n=1024 | IBAGS: n=1024")
elseif(SECURITY_LEVEL STREQUAL "3")
    add_compile_definitions(LATTICE_LEVEL_3)
    message(STATUS ">>> NIST Level 3 — plain-LWE: n=768 | IBAGS: n=1024 (Dilithium3)")
elseif(SECURITY_LEVEL STREQUAL "2")
    add_compile_definitions(LATTICE_LEVEL_1)
    message(STATUS ">>> NIST Level 2 — plain-LWE: n=512 | IBAGS: n=512 (Dilithium2)")
elseif(SECURITY_LEVEL STREQUAL "1")
    add_compile_definitions(LATTICE_LEVEL_1)
    message(STATUS ">>> NIST Level 1 — plain-LWE: n=512 | IBAGS: n=512 (Level 2)")
else()
    message(STATUS ">>> Demo (n=64) — plain-LWE q=8191, IBAGS q=7681")
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

# ============================================================================
# IBAGS (RLWE/NTRU) — 核心库
# ============================================================================
set(IBAGS_CORE_SRCS
    src_RLWEorNTRU/params.cpp
    src_RLWEorNTRU/secure_memory.cpp
    src_RLWEorNTRU/mod_reduce.cpp
    src_RLWEorNTRU/poly.cpp
    src_RLWEorNTRU/poly_ntt.cpp
    src_RLWEorNTRU/xof.cpp
    src_RLWEorNTRU/csprng.cpp
    src_RLWEorNTRU/encode.cpp
    src_RLWEorNTRU/transcript.cpp
    src_RLWEorNTRU/hash.cpp
    src_RLWEorNTRU/gauss.cpp
    src_RLWEorNTRU/reject.cpp
    src_RLWEorNTRU/ntru_secret.cpp
    src_RLWEorNTRU/ntru_trapgen.cpp
)
add_library(ibags_core STATIC ${IBAGS_CORE_SRCS})
target_include_directories(ibags_core PUBLIC ${INC_RLWEOR_NTRU_DIR})
target_compile_options(ibags_core PRIVATE -Wall -Wextra)

# OpenSSL for SHAKE256 XOF
find_package(OpenSSL REQUIRED)
target_link_libraries(ibags_core PUBLIC OpenSSL::Crypto)

# ============================================================================
# IBAGS 统一测试入口 (串联 core / poly / ntru_trapgen)
# ============================================================================

# 统一测试: 三个套件合并为单一可执行文件
add_executable(test_ibags_all
    src_RLWEorNTRU/test_ibags_all.cpp
    src_RLWEorNTRU/test_ibags_core.cpp
    src_RLWEorNTRU/test_ibags_poly.cpp
    src_RLWEorNTRU/test_ntru_trapgen.cpp
)
target_include_directories(test_ibags_all PRIVATE ${INC_RLWEOR_NTRU_DIR})
target_compile_options(test_ibags_all PRIVATE -Wall -Wextra)
target_link_libraries(test_ibags_all PRIVATE ibags_core)

# ============================================================================
# IBAGS 基准测试
# ============================================================================
add_executable(bench_ibags
    src_RLWEorNTRU/bench_ibags.cpp
    src_RLWEorNTRU/bench_utils.cpp
)
target_include_directories(bench_ibags PRIVATE ${INC_RLWEOR_NTRU_DIR})
target_compile_options(bench_ibags PRIVATE -Wall -Wextra)
target_link_libraries(bench_ibags PRIVATE ibags_core)


