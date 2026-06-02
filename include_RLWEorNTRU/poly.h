#pragma once
/**
 * @file poly.h
 * @brief R_p = Z_p[Y]/(Y^n + 1) 多项式类型
 *
 * 设计要求:
 *  - Poly 持有系数数组和维度 n，来源于 Params
 *  - 所有算术操作显式传入 Params
 *  - 支持 canonical ([0, q)) 和 centered ((−q/2, q/2]) 两种表示
 *  - 操作完成后输出 canonical
 *
 * 参考:
 *  - NFLlib 的 polynomial abstraction
 *  - PQClean 中 poly 结构体设计
 *  - Dilithium/Falcon reference 中 poly 类型划分
 */

#include "params.h"
#include "mod_reduce.h"

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>

namespace ibags {

// ============================================================================
// Poly — R_p = Z_p[Y]/(Y^n + 1) 多项式
// ============================================================================

class Poly {
public:
    // ── 构造 ──

    /// 零多项式，维度来自 Params.n
    explicit Poly(int n);

    /// 从系数列表构造，自动检查 size == n（Debug 断言）
    Poly(int n, std::vector<int64_t> coeffs);

    /// 零多项式（工厂函数）
    static Poly zero(const Params& pp);

    /// 从 int64 系数构造（不做规约，调用者负责 canonical）
    static Poly from_coeffs(const Params& pp, std::vector<int64_t> coeffs);

    /// 从 canonical [0, q) 系数构造
    static Poly from_canonical(const Params& pp, std::vector<int64_t> coeffs);

    // ── 访问 ──

    int n() const noexcept { return n_; }
    const std::vector<int64_t>& coeffs() const noexcept { return coeffs_; }
    int64_t operator[](int i) const;
    int64_t& operator[](int i);

    /// 原始系数指针（用于 NTT 等）
    const int64_t* data() const noexcept { return coeffs_.data(); }
    int64_t* data() noexcept { return coeffs_.data(); }

    // ── 状态 ──

    /// 判断所有系数是否在 [0, q) 内
    bool is_canonical(const Params& pp) const;

    /// 判断所有系数是否在 (−q/2, q/2] 内
    bool is_centered(const Params& pp) const;

    /// 将所有系数规约到 [0, q)
    void reduce_mod_q(const Params& pp);

    /// 将所有系数规约到 (−q/2, q/2]，返回新 Poly
    Poly normalize_centered(const Params& pp) const;

    // ── 移动 / 拷贝 ──
    Poly(const Poly&) = default;
    Poly& operator=(const Poly&) = default;
    Poly(Poly&&) noexcept = default;
    Poly& operator=(Poly&&) noexcept = default;

private:
    int n_;
    std::vector<int64_t> coeffs_;
};

// ============================================================================
// 自由函数 — 模约减 / 规范化
// ============================================================================

/// 将 poly 所有系数规约到 [0, q)（原地）
void poly_reduce(Poly& poly, const Params& pp);

/// 环约减: 将长度 ≤ 2n-1 的 conv 结果折叠到 R_p
/// 输入长度为 2n-1（来自卷积），输出 canonical Poly
Poly poly_ring_reduce_raw(const Params& pp, const std::vector<int64_t>& conv);

// ============================================================================
// 自由函数 — 算术运算
// ============================================================================

/// polynomial addition: out = a + b  (mod q, Y^n+1)
Poly poly_add(const Poly& a, const Poly& b, const Params& pp);

/// polynomial subtraction: out = a − b  (mod q, Y^n+1)
Poly poly_sub(const Poly& a, const Poly& b, const Params& pp);

/// additive inverse: out = −a  (mod q, Y^n+1)
Poly poly_neg(const Poly& a, const Params& pp);

/// scalar multiplication: out = s * a  (mod q, Y^n+1)
/// 逐系数乘以标量 s 后 Barrett 约减，O(n)
Poly poly_mul_scalar(const Poly& a, int64_t scalar, const Params& pp);

/// Naive O(n^2) negacyclic convolution (仅 PoC/test)
/// 使用 Y^n = -1 规则
Poly poly_mul_naive(const Poly& a, const Poly& b, const Params& pp);

// ============================================================================
// 自由函数 — 范数 / 比较
// ============================================================================

/// 计算 centered representation 下的最大系数绝对值
/// 防止 abs(INT64_MIN) 溢出
uint64_t poly_norm_inf(const Poly& a, const Params& pp);

/// 判断 ||a||∞ ≤ bound
/// 用 constant-time 风格实现，避免早停泄露
bool poly_norm_bound_check(const Poly& a, uint64_t bound, const Params& pp);

/// constant-time 多项式相等比较
/// 不使用 early return，比较所有系数
bool poly_equal_ct(const Poly& a, const Poly& b);

// ============================================================================
// 自由函数 — 求逆
// ============================================================================

/// 在 R_p 中计算可逆多项式的逆
/// 满足 a * inv ≡ 1 (mod q, Y^n+1)
/// 使用扩展欧几里得算法（NTRU 风格）
/// 检测不可逆并返回错误
Status poly_inv(const Poly& a, const Params& pp, Poly* inv);

/// 检测多项式在 R_p = Z_q[x]/(x^n+1) 中是否可逆
/// 不计算逆，仅返回 true/false
/// 可用于 TrapGen 快速可逆性检查
[[nodiscard]] bool poly_is_invertible_mod_q(const Poly& a, const Params& pp);

} // namespace ibags
