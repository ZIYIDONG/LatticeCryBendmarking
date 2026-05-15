#pragma once
#include "mp12_plain-LWE.h"
#include "decrypt_plain-LWE.h"

/**
 * UnifiedParams — compile-time parameter selector
 *
 * Security levels (CMakeLists.txt):
 *   cmake ..                                     → Demo       (n=8,   q=257,       b=2)
 *   cmake .. -DSECURITY_LEVEL=1                  → NIST L1    (n=512, q=134219777, b=27)
 *   cmake .. -DSECURITY_LEVEL=3                  → NIST L3    (n=768, q=4294967291, b=32)
 *   cmake .. -DSECURITY_LEVEL=5                  → NIST L5    (n=1024,q=4294967311, b=32)
 *
 * CMake injects: LATTICE_DEMO / LATTICE_LEVEL_1 / LATTICE_LEVEL_3 / LATTICE_LEVEL_5
 */
namespace unified {

/* ════════════════════════════════════════════════════
   MP12 核心参数
   ════════════════════════════════════════════════════ */
inline mp12::Params default_mp12_params() {
#if defined(LATTICE_LEVEL_5)
    constexpr int  N = 1024;
    constexpr long Q = 4294967311L;   // ~2^32, prime
    constexpr int  B = 32;
#elif defined(LATTICE_LEVEL_3)
    constexpr int  N = 768;
    constexpr long Q = 4294967291L;   // ~2^32, prime
    constexpr int  B = 32;
#elif defined(LATTICE_LEVEL_1)
    constexpr int  N = 512;
    constexpr long Q = 134219777;     // ~2^27, prime
    constexpr int  B = 27;
#else  // LATTICE_DEMO (default)
    constexpr int  N = 8;
    constexpr long Q = 257;           // prime
    constexpr int  B = 2;
#endif
    return mp12::Params::make(N, Q, B);
}

/* ════════════════════════════════════════════════════
   多身份阈值参数
   ════════════════════════════════════════════════════ */
inline cryptolib::MIDParams default_midparams_128(int d = 1, int N_id = 3) {
    auto p = default_mp12_params();
#if defined(LATTICE_LEVEL_5) || defined(LATTICE_LEVEL_3) || defined(LATTICE_LEVEL_1)
    constexpr int LAMBDA = 128;
    constexpr int B_CHI  = 1;
#else
    constexpr int LAMBDA = 8;
    constexpr int B_CHI  = 1;
#endif
    return cryptolib::MIDParams::make(p.n, d, p.q, N_id, p.b, /*lambda=*/LAMBDA, /*B_chi=*/B_CHI);
}

} // namespace unified
