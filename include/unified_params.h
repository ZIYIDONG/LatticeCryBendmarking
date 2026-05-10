#pragma once
#include "mp12.h"
#include "decrypt.h"

#if __has_include("params_cfg.h")
#  include "params_cfg.h"
#  define UNIFIED_CFG_AVAILABLE 1
#else
#  define UNIFIED_CFG_AVAILABLE 0
#endif

namespace unified {

inline mp12::Params default_mp12_params_128() {
#if UNIFIED_CFG_AVAILABLE
    return mp12::Params::make(MP12_N, MP12_Q, MP12_B);
#else
    // Conservative 128-bit demo parameters (fallback)
    return mp12::Params::make(512, 134219777, 27);
#endif
}

inline cryptolib::MIDParams default_midparams_128(int d = 1, int N_id = 3) {
#if UNIFIED_CFG_AVAILABLE
    mp12::Params p = default_mp12_params_128();
    return cryptolib::MIDParams::make(p.n, d, p.q, N_id, p.b, /*lambda=*/MID_LAMBDA, /*B_chi=*/1);
#else
    auto p = default_mp12_params_128();
    // lambda set to 128 for 128-bit security; B_chi left small for demo
    return cryptolib::MIDParams::make(p.n, d, p.q, N_id, p.b, /*lambda=*/128, /*B_chi=*/1);
#endif
}

} // namespace unified
