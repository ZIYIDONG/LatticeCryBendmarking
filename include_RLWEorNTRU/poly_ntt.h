#pragma once
/**
 * @file poly_ntt.h
 * @brief NTT-domain 多项式表示及 NTT 变换
 *
 * 设计原则:
 *  - NttPoly 是 NTT 域的独立类型，不与普通 Poly 混用
 *  - NttTable 预计算所有 NTT 所需数据 (twiddle factors, bitrev, 约减常数)
 *    调用方创建一次、多次复用，避免每调用重新计算
 *  - 仅对 NTT-friendly 的 q 有效 (q ≡ 1 mod 2n)
 *  - 支持正向 NTT、逆向 INTT、逐点乘法
 *
 * 参考:
 *  - Dilithium ntt.c / invntt.c
 *  - PQClean ML-DSA NTT 实现
 *  - Kyber ntt.c
 */

#include "params.h"
#include "errors.h"
#include "mod_reduce.h"

#include <cstdint>
#include <vector>

namespace ibags {

// ============================================================================
// NttPoly — NTT 域多项式
// ============================================================================

class NttPoly {
public:
    explicit NttPoly(int n);
    explicit NttPoly(const Params& pp);

    int n() const noexcept { return n_; }
    const int64_t* data() const noexcept { return coeffs_.data(); }
    int64_t* data() noexcept { return coeffs_.data(); }
    int64_t operator[](int i) const;
    int64_t& operator[](int i);

    NttPoly(const NttPoly&) = default;
    NttPoly& operator=(const NttPoly&) = default;
    NttPoly(NttPoly&&) noexcept = default;
    NttPoly& operator=(NttPoly&&) noexcept = default;

private:
    int n_;
    std::vector<int64_t> coeffs_;
};

// ============================================================================
// NttTable — NTT 预计算表 (一次计算，多次复用)
// ============================================================================

/**
 * @struct NttTable
 * @brief 持有特定参数集下 NTT/INTT 所有预计算数据。
 *
 * 字段:
 *  - zetas[i]  = ψ^i          (i = 0..2n-1), ψ 为 primitive 2n-th root of unity
 *  - zetas_inv[i] = ψ^{-i}    (i = 0..2n-1)
 *  - bitrev[i]   = bit-reversal permutation of i  (i = 0..n-1)
 *  - n_inv       = n^{-1} mod q  (用于 INTT 缩放)
 *  - barrett     = Barrett 约减预计算常数
 *  - montgomery  = Montgomery 约减预计算常数 (可选, 用于优化的逐点乘法)
 */
struct NttTable {
    int n;
    uint32_t q;
    int log_n;

    std::vector<int64_t> zetas;       ///< ψ^i, i = 0..2n-1
    std::vector<int64_t> zetas_inv;   ///< ψ^{-i}, i = 0..2n-1
    std::vector<int>     bitrev;      ///< bit-reversal permutation, 0..n-1
    int64_t              n_inv;       ///< n^{-1} mod q

    BarrettConst         barrett;     ///< Barrett 约减常数
    MontgomeryConst      montgomery;  ///< Montgomery 约减常数

    /**
     * @brief 工厂: 从参数集构造完整的 NTT 预计算表。
     *
     * 前置条件: is_ntt_friendly(pp) == true (q ≡ 1 mod 2n)
     * 若参数不满足 NTT 友好条件则抛出 std::invalid_argument。
     *
     * 算法:
     *  1. 寻找 primitive 2n-th root of unity ψ
     *  2. 计算 zetas = ψ^i, zetas_inv = ψ^{-i}
     *  3. 预计算 bit-reversal 索引表
     *  4. 计算 n^{-1} mod q
     *  5. 预计算 Barrett/Montgomery 约减常数
     */
    static NttTable create(const Params& pp);
};

// ============================================================================
// NTT / INTT 接口
// ============================================================================

/// 正向 NTT: normal domain → NTT domain
/// 使用预计算的 NttTable
Status poly_ntt(const class Poly& a, const NttTable& tbl, NttPoly* out);

/// 逆 NTT: NTT domain → normal domain (canonical)
Status poly_invntt(const NttPoly& a, const NttTable& tbl, class Poly* out);

/// NTT 域逐点乘法
/// 输入必须是 NTT-domain type，输出仍为 NTT-domain
Status poly_pointwise_mul_ntt(const NttPoly& a, const NttPoly& b,
                              const NttTable& tbl, NttPoly* out);

/// 判断 q 是否满足 NTT-friendly 条件: q ≡ 1 (mod 2n)
bool is_ntt_friendly(const Params& pp);

/// NTT-based 多项式乘法: out = a * b in Z_q[X]/(X^n+1)
/// 完整管线: NTT(a) + NTT(b) + pointwise_mul + INTT
Status poly_mul_ntt(const class Poly& a, const class Poly& b,
                    const NttTable& tbl, class Poly* out);

/// 检查 NTT 实现的 roundtrip 一致性
/// INTT(NTT(a)) == a (带 tolerance)
Status ntt_roundtrip_test(const class Poly& a, const NttTable& tbl);

} // namespace ibags
