#pragma once
#include "mp12.h"
#include "decrypt.h"

namespace unified {

inline mp12::Params default_mp12_params_128() {
    // Conservative 128-bit demo parameters (used in repository comment)
    return mp12::Params::make(512, 134219777, 27);
}

inline cryptolib::MIDParams default_midparams_128(int d = 1, int N_id = 3) {
    auto p = default_mp12_params_128();
    // lambda set to 128 for 128-bit security; B_chi left small for demo
    return cryptolib::MIDParams::make(p.n, d, p.q, N_id, p.b, /*lambda=*/128, /*B_chi=*/1);
}

} // namespace unified
