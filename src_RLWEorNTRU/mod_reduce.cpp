/**
 * @file mod_reduce.cpp
 * @brief Barrett (乘-移位) / Montgomery (REDC) 模约减实现
 *
 * Barrett 算法 (Handbook of Applied Cryptography §14.3.3):
 *   对模数 q，预计算 mu = floor(2^k / q)，取 k = 64。
 *   给定 a < q * 2^32，计算:
 *     q_est = floor(a * mu / 2^64)
 *     r = a - q_est * q
 *     if (r >= q) r -= q
 *   此时 r = a mod q 且 r ∈ [0, q)。
 *   q_est 与真实商 floor(a/q) 相差至多 1，保证一次修正足矣。
 *
 * Montgomery REDC (Handbook of Applied Cryptography §14.3.2):
 *   选 R = 2^32，预计算 q_inv = -q^{-1} mod R。
 *   给定 T < q * R，计算:
 *     m = (T mod R) * q_inv mod R
 *     t = (T + m * q) / R
 *     if (t >= q) t -= q
 *   此时 t = T * R^{-1} mod q。
 *
 * 正确性要求: 所有 q < 2^30。
 *   Barrett: |a| < q * 2^32 ⇒ 中间乘法用 unsigned __int128 安全。
 *   Montgomery: T < q * R = q * 2^32 ⇒ m*q < 2^62, T + m*q < 2^63, 均 < 2^64。
 *   乘法 a*b (用于 montgomery_mul): a,b < q < 2^30 ⇒ a*b < 2^60，用 __int128 安全。
 *
 * 参考:
 *  - Dilithium reduce.c / reduce.h
 *  - Kyber/ML-KEM reduce.c
 *  - PQClean 各方案的 reduce 实现
 */

#include "../include_RLWEorNTRU/mod_reduce.h"
#include <cstdint>
#include <cassert>
#include <cstring>

// __int128 支持检查 (GCC/Clang, 非 MSVC)
#if !defined(__SIZEOF_INT128__)
#  error "mod_reduce requires __int128 support (GCC/Clang on 64-bit)"
#endif

namespace ibags {

// ============================================================================
// Barrett 约减 — 实现
// ============================================================================

BarrettConst make_barrett(uint32_t q) {
    assert(q > 0 && q < (1u << 30) && "q must be in (0, 2^30)");
    uint64_t mu = (static_cast<unsigned __int128>(1) << 64) / q;
    return {q, mu};
}

int64_t barrett_reduce(int64_t a, const BarrettConst& bc) {
    // 对 |a| 做 Barrett 约减，再根据符号调整。
    // Barrett 正确性要求 |a| < q * 2^32。若 |a| 超界则回退到 % q。
    constexpr uint64_t BOUND = (uint64_t)1 << 32;
    bool neg = (a < 0);
    uint64_t abs_a;
    if (neg) {
        abs_a = static_cast<uint64_t>(-(a + 1)) + 1;
    } else {
        abs_a = static_cast<uint64_t>(a);
    }

    // 若值过大 (如环约减的 conv 累加值), 回退到硬件除法
    if (abs_a >= BOUND * bc.q) {
        int64_t r = static_cast<int64_t>(abs_a % bc.q);
        if (neg && r != 0) r = static_cast<int64_t>(bc.q) - r;
        return r;
    }

    // Barrett 商估计: q_est = floor(abs_a * mu / 2^64)
    uint64_t q_est = static_cast<uint64_t>(
        (static_cast<unsigned __int128>(abs_a) * bc.mu) >> 64);

    uint64_t r = abs_a - q_est * bc.q;

    if (r >= bc.q) {
        r -= bc.q;
    }

    if (neg && r != 0) {
        r = bc.q - r;
    }
    return static_cast<int64_t>(r);
}

int64_t barrett_reduce(int64_t a, uint32_t q) {
    // 便捷版：每次构造 Barrett 常数，仅用于非热点路径
    return barrett_reduce(a, make_barrett(q));
}

// ============================================================================
// Constant-time Barrett 约减
// ============================================================================

int64_t barrett_reduce_ct(int64_t a, const BarrettConst& bc) {
    // 常数时间实现：所有分支用算术掩码消除。
    // 使用 |a| + 符号恢复方案，正确处理 INT64_MIN。
    const uint64_t sign_mask = static_cast<uint64_t>(a >> 63);
    // 常数时间绝对值: abs(x) = (x ^ mask) - mask, 对 INT64_MIN 也安全
    const uint64_t abs_a = (static_cast<uint64_t>(a) ^ sign_mask) - sign_mask;

    // Barrett on abs_a
    const uint64_t q_est = static_cast<uint64_t>(
        (static_cast<unsigned __int128>(abs_a) * bc.mu) >> 64);
    uint64_t r = abs_a - q_est * bc.q;

    // 条件减法: r >= bc.q ? r -= bc.q
    // mask_ge = 0  if r < q,  ~0 if r >= q
    const uint64_t mask_ge = ~((r - bc.q) >> 63);
    r -= mask_ge & bc.q;

    // 符号恢复: 若原值为负且 r != 0, 则 r = q - r
    const uint64_t r_nonzero = 1 ^ ((r | -r) >> 63);  // r != 0 → 1
    const uint64_t adjust = sign_mask & r_nonzero;
    // r = adjust ? (q - r) : r
    r = (bc.q - r) * adjust + r * (1 - adjust);

    return static_cast<int64_t>(r);
}

int64_t barrett_reduce_ct(int64_t a, uint32_t q) {
    return barrett_reduce_ct(a, make_barrett(q));
}

// ============================================================================
// Montgomery 约减 — 实现
// ============================================================================

MontgomeryConst make_montgomery(uint32_t q) {
    assert(q > 0 && q < (1u << 30) && "q must be in (0, 2^30)");
    assert((q & 1) && "q must be odd for Montgomery reduction");

    // q_inv = -q^{-1} mod 2^32
    // 使用 Newton 迭代: q_inv = q_inv * (2 - q * q_inv) mod 2^k
    // 初始: q_inv = q (mod 2^3 成立因为 q*q ≡ 1 mod 8 对奇数 q 总成立)
    uint32_t q_inv = q;
    for (int i = 0; i < 4; ++i) {
        q_inv = q_inv * (2u - q * q_inv);
    }
    // 取负: q_inv = -q^{-1} mod 2^32
    q_inv = -q_inv;

    // R^2 mod q = (2^32)^2 mod q = 2^64 mod q
    uint64_t r2 = (static_cast<unsigned __int128>(1) << 64) % q;

    return {q, r2, q_inv};
}

int64_t to_montgomery(int64_t a, const MontgomeryConst& mc) {
    // a_mont = a * R mod q
    // 计算 a * R^2，再用 REDC 除掉一个 R 因子:
    //   REDC(a * R^2) = a * R^2 * R^{-1} = a * R (mod q)
    unsigned __int128 T = static_cast<unsigned __int128>(a) * mc.r_sq_mod_q;
    return montgomery_reduce(static_cast<int64_t>(static_cast<uint64_t>(T)), mc);
}

int64_t from_montgomery(int64_t a_mont, const MontgomeryConst& mc) {
    // a = a_mont * R^{-1} mod q = REDC(a_mont)
    return montgomery_reduce(a_mont, mc);
}

int64_t montgomery_reduce(int64_t T, const MontgomeryConst& mc) {
    // 标准 REDC: 输入 T < q * R，输出 T * R^{-1} mod q
    // 要求 T 在 [0, 2^63) 范围内 (int64 非负一半)
    // NTT 场景中 T ∈ [0, q*R) ⊂ [0, 2^62)，未触及 int64 上限。
    uint64_t uT = static_cast<uint64_t>(T);

    // m = (T mod R) * q_inv mod R
    uint32_t m = static_cast<uint32_t>(uT) * mc.q_inv;

    // t = (T + m * q) / R
    // m < 2^32, q < 2^30 ⇒ m * q < 2^62, T + m*q < 2^63 < 2^64 ✓
    uint64_t t = (uT + static_cast<uint64_t>(m) * mc.q) >> 32;

    if (t >= mc.q) {
        t -= mc.q;
    }
    return static_cast<int64_t>(t);
}

int64_t montgomery_mul(int64_t a_mont, int64_t b_mont, const MontgomeryConst& mc) {
    // Montgomery 乘法: a_mont * b_mont * R^{-1} mod q
    // a_mont = a_true * R mod q, b_mont = b_true * R mod q
    // REDC(a_mont * b_mont) = a_true * b_true * R^2 * R^{-1} = a_true * b_true * R mod q ✓

    uint64_t ua = static_cast<uint64_t>(a_mont);
    uint64_t ub = static_cast<uint64_t>(b_mont);
    unsigned __int128 prod = static_cast<unsigned __int128>(ua) * ub;

    // REDC(prod)
    uint64_t prod_lo = static_cast<uint64_t>(prod);
    uint32_t m = static_cast<uint32_t>(prod_lo) * mc.q_inv;

    unsigned __int128 t128 = prod + static_cast<unsigned __int128>(m) * mc.q;
    uint64_t t = static_cast<uint64_t>(t128 >> 32);

    if (t >= mc.q) {
        t -= mc.q;
    }
    return static_cast<int64_t>(t);
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
