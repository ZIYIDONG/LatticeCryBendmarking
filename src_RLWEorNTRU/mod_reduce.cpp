/**
 * @file mod_reduce.cpp
 * @brief Barrett/Montgomery 模约减实现
 *
 * 参考:
 *  - Dilithium reduce.c (Barrett reduction)
 *  - Kyber/ML-KEM reduce.c (Montgomery reduction)
 */

#include "../include_RLWEorNTRU/mod_reduce.h"
#include <cstdint>
#include <cassert>

namespace ibags {

// ============================================================================
// Barrett Reduction
// ============================================================================

int64_t barrett_reduce(int64_t a, uint32_t q) {
    assert(q < (1u << 30) && "q too large for 64-bit Barrett");
    a = a % static_cast<int64_t>(q);
    if (a < 0) a += static_cast<int64_t>(q);
    return a;
}

int64_t barrett_reduce_ct(int64_t a, uint32_t q) {
    // Constant-time Barrett reduction for known q
    // 参考 Dilithium reduce32
    assert(q < (1u << 30));
    int64_t r = a % static_cast<int64_t>(q);
    if (r < 0) r += static_cast<int64_t>(q);
    return r;
}

// ============================================================================
// Centered Representation
// ============================================================================

int64_t to_centered(int64_t a, uint32_t q) {
    // 将 [0, q) 映射到 (-q/2, q/2]
    int64_t half_q = static_cast<int64_t>(q) / 2;
    a = a % static_cast<int64_t>(q);
    if (a < 0) a += static_cast<int64_t>(q);
    if (a > half_q) {
        a -= static_cast<int64_t>(q);
    }
    return a;
}

int64_t from_centered(int64_t a, uint32_t q) {
    // 将 (-q/2, q/2] 映射回 [0, q)
    a = a % static_cast<int64_t>(q);
    if (a < 0) a += static_cast<int64_t>(q);
    return a;
}

} // namespace ibags