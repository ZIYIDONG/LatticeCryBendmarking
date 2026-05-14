#pragma once
#include "mp12_plain-LWE.h"
#include "decrypt_plain-LWE.h"

/**
 * UnifiedParams — 零中间文件参数提供者
 *
 * 切换方式 (CMakeLists.txt):
 *   cmake ..                                  → 演示参数 (n=8,  q=257,  b=2)
 *   cmake .. -DUSE_128BIT_PARAMS=ON           → 128-bit (n=512, q=134219777, b=27)
 *
 * 原理: CMake 通过 add_compile_definitions(LATTICE_128BIT)
 *       将宏注入所有编译单元的预处理器，无需 configure_file
 */
namespace unified {

/* ════════════════════════════════════════════════════
   MP12 核心参数
   ════════════════════════════════════════════════════ */
inline mp12::Params default_mp12_params() {
#ifdef LATTICE_128BIT
    constexpr int  N = 512;
    constexpr long Q = 134219777;
    constexpr int  B = 27;
#else
    constexpr int  N = 8;
    constexpr long Q = 257;
    constexpr int  B = 2;
#endif
    return mp12::Params::make(N, Q, B);
}

/* ════════════════════════════════════════════════════
   多身份阈值参数
   ════════════════════════════════════════════════════ */
inline cryptolib::MIDParams default_midparams_128(int d = 1, int N_id = 3) {
    auto p = default_mp12_params();
#ifdef LATTICE_128BIT
    constexpr int  LAMBDA = 128;
    constexpr int  B_CHI  = 1;
#else
    constexpr int  LAMBDA = 8;
    constexpr int  B_CHI  = 1;
#endif
    return cryptolib::MIDParams::make(p.n, d, p.q, N_id, p.b, /*lambda=*/LAMBDA, /*B_chi=*/B_CHI);
}

} // namespace unified
