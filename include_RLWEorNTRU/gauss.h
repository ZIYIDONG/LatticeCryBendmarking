#pragma once
/**
 * @file gauss.h
 * @brief 离散高斯采样 — 基于 SHAKE256 的整数高斯采样
 *
 * 算法:
 *   - GaussSamp(R, M): 确定性伪随机离散高斯采样器
 *   - 使用 XOF (SHAKE256) 拒绝采样生成 E_bits × M 个样本
 *   - 输出 SmallPoly（系数域在 [-B, B]）
 *
 * 核心原语:
 *   1. gauss_sample_coeff(xof, sigma): 从 XOF 采样单个离散高斯整数
 *   2. gauss_sample_poly(xof, pp): 采样整体系数多项式
 *   3. gauss_sample_seeds(pp, R, M): 为 M 个签名者生成确定性种子流
 *
 * 设计原则:
 *   - 确定性: same R, same pp → same output (用于签名者一致性)
 *   - 无外部随机依赖: 所有随机性来自 XOF
 *   - 紧凑输出: SmallPoly 使用 int16_t 系数存储
 *
 * 参考:
 *   - Falcon 的 fpr_sampler (Gaussian_sampler)
 *   - Dilithium 的 poly_uniform_gamma1 (拒绝采样)
 *   - GPV08 离散高斯采样器 (Klein 算法)
 */

#include "xof.h"
#include "csprng.h"
#include "domain.h"
#include "params.h"
#include "poly.h"
#include "errors.h"

#include <cstdint>
#include <cstddef>
#include <vector>

namespace ibags {

// ============================================================================
// GaussFunc — 小系数多项式 (紧致表示)
// ============================================================================

/**
 * @brief 小系数多项式 — 使用 int16_t 存储
 *
 * 系数范围受限于 [−B, B]，通常 B ≤ 2^10
 * 用于 NTRU 密钥/签名，需先提升到 Poly (int64_t) 才能做环运算
 */
class SmallPoly {
public:
    SmallPoly() : n_(0) {}

    explicit SmallPoly(int n) : n_(n), coeffs_(n, 0) {}

    SmallPoly(int n, std::vector<int16_t> coeffs)
        : n_(n), coeffs_(std::move(coeffs)) {}

    static SmallPoly zero(int n) {
        return SmallPoly(n);
    }

    int n() const noexcept { return n_; }
    const std::vector<int16_t>& coeffs() const noexcept { return coeffs_; }
    int16_t operator[](int i) const;
    int16_t& operator[](int i);

    /// 提升到 Poly (int64_t)
    Poly to_poly() const;

    SmallPoly(const SmallPoly&) = default;
    SmallPoly& operator=(const SmallPoly&) = default;
    SmallPoly(SmallPoly&&) noexcept = default;
    SmallPoly& operator=(SmallPoly&&) noexcept = default;

private:
    int n_;
    std::vector<int16_t> coeffs_;
};

// ============================================================================
// 离散高斯采样原语
// ============================================================================

/**
 * @brief 从 XOF 采样单个离散高斯整数 (sigma 给定)
 *
 * 实现: 基于 sigma 的表驱动拒绝采样
 * 每个样本消耗若干个 squeeze_u64 调用 (熵效率高)
 *
 * @param xof     已初始化的 XOF (可处于 absorb 或 squeeze 模式)
 * @param sigma   高斯宽度
 * @param q       模数 (用于模约减)
 * @return        整数样本 ∈ [−cutoff, cutoff]
 */
int64_t gauss_sample_coeff(Xof& xof, double sigma, int q);

/**
 * @brief 从 XOF 采样整体系数多项式 (coefficient-wise)
 *
 * 为 pp.n 个系数分别调用 gauss_sample_coeff
 *
 * @param xof  已初始化的 XOF
 * @param pp   参数集
 * @return     小系数多项式
 */
SmallPoly gauss_sample_poly(Xof& xof, const Params& pp);

// ============================================================================
// GaussFunc: 确定性种子流生成
// ============================================================================

/**
 * @brief GaussFunc 函数: 从主种子生成 M 个签名者的种子列表
 *
 * 协议用法:
 *   每个签名者 j 获得 seed[j] = GaussFunc(R, M)[j]
 *   签名者使用 seed[j] 派生自己的 alpha_j 承诺
 *
 * 实现:
 *   XOF(domain=GAUSS_FUNCTION, seed=R) → squeeze M × seed_len bytes
 *
 * 确定性保证:
 *   相同 (R, M, pp) → 相同输出（所有诚实签名者一致）
 *
 * @param pp           参数集
 * @param R            共享主种子 (32 字节)
 * @param M            签名者总数
 * @param out_seeds    输出种子列表 (每个 32 字节)
 */
void gauss_gen_seeds(const Params& pp,
                     ByteSpan R,
                     int M,
                     std::vector<std::vector<uint8_t>>& out_seeds);

/**
 * @brief 从单个种子膨胀为签名者 j 的 alpha_j 承诺多项式
 *
 * 算法:
 *   1. XOF(domain=GAUSS_FUNCTION, seed=seed_j)
 *   2. 采样小系数多项式 alpha_j = gauss_sample_poly(xof, pp)
 *
 * @param pp      参数集
 * @param seed_j  签名者 j 的种子 (32 字节)
 * @return        alpha_j 小系数多项式
 */
SmallPoly gauss_expand_alpha(const Params& pp, ByteSpan seed_j);

// ============================================================================
// 辅助函数: build_gauss_xof
// ============================================================================

/**
 * @brief 构造用于高斯采样的 XOF，吸收 domain + seed
 *
 * @param domain  domain label
 * @param seed    种子字节
 * @return        已吸收 domain + seed 的 XOF (未 finalized)
 */
Xof build_gauss_xof(std::string_view domain, ByteSpan seed);

} // namespace ibags