#pragma once
/**
 * @file mod_reduce.h
 * @brief Barrett (乘-移位) / Montgomery (REDC) 模约减与 centered representation
 *
 * Barrett 算法:
 *   预计算 mu = floor(2^64 / q)，约减时 t = floor(a * mu / 2^64)，
 *   余数 r = a - t*q，最多一次条件减法。正确性要求 |a| < q * 2^32。
 *
 * Montgomery 算法:
 *   选择 R = 2^32。预计算 q_inv = -q^{-1} mod R，将整数映射到 Montgomery 域
 *   a_mont = a * R mod q。Montgomery 乘法计算 a_mont * b_mont * R^{-1} mod q。
 *
 * Dilithium 风格常数时间:
 *   使用符号掩码消除分支，对 secret-dependent 数据路径保证无 timing leak。
 *
 * 参考:
 *  - Dilithium reduce.c / reduce.h (signed Barrett & Montgomery)
 *  - Kyber/ML-KEM reduce.c (Montgomery reduction)
 *  - Handbook of Applied Cryptography §14.3.3 (Barrett), §14.3.2 (Montgomery)
 */

#include <cstdint>

namespace ibags {

// ============================================================================
// Barrett 约减 — 预计算常数 + 约减函数
// ============================================================================

/// Barrett 预计算常数
struct BarrettConst {
    uint32_t q;    ///< 模数 (q < 2^30)
    uint64_t mu;   ///< floor(2^64 / q)
};

/// 从模数 q 构造 Barrett 常数
/// 前置条件: 0 < q < 2^30
BarrettConst make_barrett(uint32_t q);

/// Barrett 约减: a mod q → [0, q)
/// 正确性要求: |a| < q * 2^32（通常远小于此，NTT 中间值必定满足）
/// 参数 a 可为负值
int64_t barrett_reduce(int64_t a, const BarrettConst& bc);

/// Barrett 约减 (无预计算便捷版): a mod q → [0, q)
/// 每次调用内部构造 BarrettConst，仅用于非热点路径或测试
/// 热点路径请使用 barrett_reduce(a, bc) 并复用 BarrettConst
int64_t barrett_reduce(int64_t a, uint32_t q);

/// 常数时间 Barrett 约减: a mod q → [0, q)
/// 不使用数据依赖分支，适用于 secret-dependent 输入
int64_t barrett_reduce_ct(int64_t a, const BarrettConst& bc);

/// 常数时间 Barrett 约减 (无预计算便捷版)
int64_t barrett_reduce_ct(int64_t a, uint32_t q);

// ============================================================================
// Montgomery 约减 — 预计算常数 + 约减函数
// ============================================================================

/// Montgomery 预计算常数 (R = 2^32)
struct MontgomeryConst {
    uint32_t q;           ///< 模数 (q < 2^30)
    uint64_t r_sq_mod_q;  ///< R^2 mod q, 用于进入 Montgomery 域
    uint32_t q_inv;       ///< -q^{-1} mod 2^32 (即 q * q_inv ≡ -1 mod 2^32)
};

/// 从模数 q 构造 Montgomery 常数
/// 前置条件: 0 < q < 2^30 且 q 为奇数
MontgomeryConst make_montgomery(uint32_t q);

/// 将普通整数 a [0, q) 映射到 Montgomery 域: a * R mod q
int64_t to_montgomery(int64_t a, const MontgomeryConst& mc);

/// 将 Montgomery 域整数 a_mont 映射回普通域: a_mont * R^{-1} mod q
int64_t from_montgomery(int64_t a_mont, const MontgomeryConst& mc);

/// Montgomery 约减 (REDC): 给定 T < q * R，计算 T * R^{-1} mod q
/// 这是 Montgomery 乘法的核心原语
int64_t montgomery_reduce(int64_t T, const MontgomeryConst& mc);

/// Montgomery 乘法: 给定两个 Montgomery 域整数 a_mont, b_mont，
/// 返回 (a_true * b_true) * R mod q
/// 等效于 REDC(a_mont * b_mont)
int64_t montgomery_mul(int64_t a_mont, int64_t b_mont, const MontgomeryConst& mc);

// ============================================================================
// Centered Representation
// ============================================================================

/// 将 [0, q) 映射到 (−q/2, q/2]
int64_t to_centered(int64_t a, uint32_t q);

/// 将 (−q/2, q/2] 映射回 [0, q)
int64_t from_centered(int64_t a, uint32_t q);

} // namespace ibags
