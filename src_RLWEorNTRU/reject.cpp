/**
 * @file reject.cpp
 * @brief 拒绝采样实现 — 从 XOF 产生均匀分布元素
 *
 * 两个原语:
 *   1. hash_to_uniform_coeff_reject: Z_q 拒绝采样
 *   2. coeff_reject_small:          小系数拒绝采样 [-B, B]
 */

#include "../include_RLWEorNTRU/reject.h"
#include "../include_RLWEorNTRU/xof.h"
#include "../include_RLWEorNTRU/params.h"

#include <cstdint>
#include <cstring>
#include <stdexcept>

namespace ibags {

// ============================================================================
// hash_to_uniform_coeff_reject — Z_q 拒绝采样
// ============================================================================

int64_t hash_to_uniform_coeff_reject(Xof& xof, int64_t q, int q_shift) {
    // 掩码: 覆盖 log2(q) 向上取整的位宽
    // q_shift = ceil(log2(q)), 最大 32 bits (兼容 q~2^32)
    if (q_shift <= 0 || q_shift > 32) {
        throw std::invalid_argument(
            "hash_to_uniform_coeff_reject: q_shift must be in [1,32]");
    }
    if (q <= 0) {
        throw std::invalid_argument(
            "hash_to_uniform_coeff_reject: q must be > 0");
    }

    // 构造掩码: mask = 2^{q_shift} - 1
    uint64_t mask = (q_shift >= 64) ? UINT64_MAX : ((1ULL << q_shift) - 1ULL);

    for (;;) {
        uint64_t r = xof.squeeze_u64();
        uint64_t t = r & mask;
        if (t < static_cast<uint64_t>(q)) {
            return static_cast<int64_t>(t);
        }
        // else: reject and resample
    }
}

// ============================================================================
// hash_to_uniform — 完整多项式
// ============================================================================

Poly hash_to_uniform(Xof& xof, const Params& pp) {
    if (pp.n <= 0) return Poly(0, {});

    // ceil(log2(q)) — 保守选择
    int q_shift = 0;
    int64_t tmp = pp.q;
    while (tmp > 0) {
        q_shift++;
        tmp >>= 1;
    }
    // For exact power-of-two q, subtract 1 to get correct bitwidth
    // e.g., q=8191 (2^13 - 1) → q_shift=13 correctly
    //      q=65536 (2^16) → q_shift=17, but we keep 16 since mask=2^{16}-1
    // Let's verify correctness: if q==0, already caught above
    if ((q_shift > 0) && (static_cast<int64_t>(1ULL << (q_shift - 1)) >= pp.q)) {
        // Actually, the loop above gives ceil(log2(q+1))?
        // Recalculate cleanly
        q_shift = 0;
        tmp = pp.q - 1; // largest value
        while (tmp > 0) {
            q_shift++;
            tmp >>= 1;
        }
    }

    // Recalculate q_shift properly:
    // We need smallest k such that 2^k > q (strict >)
    // Then mask = 2^k - 1. But if q is power-of-two, k = log2(q) + 1
    // Simpler: k = 0; while ((1ULL << k) <= q) k++;
    q_shift = 0;
    while ((1ULL << q_shift) <= static_cast<uint64_t>(pp.q)) {
        q_shift++;
    }
    // Now mask = 2^{q_shift} - 1 > q-1, which ensures uniform rejection with < q checks

    std::vector<int64_t> coeffs(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        coeffs[i] = hash_to_uniform_coeff_reject(xof, pp.q, q_shift);
    }
    return Poly(pp.n, std::move(coeffs));
}

// ============================================================================
// coeff_reject_small — 挑战采样 [-B, B]
// ============================================================================

int64_t coeff_reject_small(Xof& xof, int B, int mask_bits) {
    // B = 挑战界 (θ)
    // mask_bits = ceil(log2(2*B + 1))
    if (B <= 0) return 0;
    if (mask_bits <= 0 || mask_bits > 16) {
        throw std::invalid_argument(
            "coeff_reject_small: mask_bits must be in [1,16]");
    }

    uint64_t mask = (1ULL << mask_bits) - 1ULL;
    int max_val = 2 * B;

    for (;;) {
        uint64_t r = xof.squeeze_u64();
        uint64_t t = r & mask;
        if (t <= static_cast<uint64_t>(max_val)) {
            return static_cast<int64_t>(t) - B;
        }
        // else: reject and resample
    }
}

} // namespace ibags