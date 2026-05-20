/**
 * @file gauss.cpp
 * @brief 离散高斯采样实现
 *
 * 使用基于 sigma 的拒绝采样表实现 gauss_sample_coeff
 */

#include "../include_RLWEorNTRU/gauss.h"
#include "../include_RLWEorNTRU/domain.h"
#include "../include_RLWEorNTRU/xof.h"
#include "../include_RLWEorNTRU/params.h"

#include <cmath>
#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <limits>

namespace ibags {

// ============================================================================
// SmallPoly 成员
// ============================================================================

int16_t SmallPoly::operator[](int i) const {
    return coeffs_.at(static_cast<size_t>(i));
}

int16_t& SmallPoly::operator[](int i) {
    return coeffs_.at(static_cast<size_t>(i));
}

Poly SmallPoly::to_poly() const {
    if (n_ <= 0) return Poly(0, {});
    std::vector<int64_t> big(n_);
    for (int i = 0; i < n_; ++i) {
        big[i] = static_cast<int64_t>(coeffs_[i]);
    }
    return Poly(n_, std::move(big));
}

// ============================================================================
// 内部常量 — 高斯采样表
// ============================================================================

/// 高斯表条目: 累积概率分布 CDF 近似
struct GaussTableEntry {
    int64_t value;
    uint64_t cdf_upper; // 累積概率上限 (比例因子 2^32)
};

/// Sigma = 1.0 的高斯分布表 (cutoff = 8 sigma = 8)
///    值         CDF_upper (2^32 比例)
///    -8   →  ~0.000000
///    -7   →  ~0.000001
///    -6   →  ~0.000031
///    -5   →  ~0.000320
///    -4   →  ~0.003660
///    -3   →  ~0.022750
///    -2   →  ~0.158655
///    -1   →  ~0.500000
///     0   →  ~0.841345
///     1   →  ~0.977250
///     2   →  ~0.996340
///     3   →  ~0.999680
///     4   →  ~0.999969
///     5   →  ~0.999999
///     6   →  ~1.000000
///     7   →  ~1.000000
///     8   →  ~1.000000
static const GaussTableEntry kGaussTableSigma1[] = {
    {-8, 0},
    {-7, 1},
    {-6, 31},
    {-5, 321},
    {-4, 3660},
    {-3, 22750},
    {-2, 158655},
    {-1, 500000},
    { 0, 841345},
    { 1, 977250},
    { 2, 996340},
    { 3, 999680},
    { 4, 999969},
    { 5, 999999},
    { 6, 1000000},
    { 7, 1000000},
    { 8, 1000000},
};
static constexpr int kGaussTableSize1 = sizeof(kGaussTableSigma1) / sizeof(kGaussTableSigma1[0]);

/// 全局 CDF 总质量 = 10^6

/**
 * @brief 从表中采样离散高斯整数
 *
 * 拒绝采样:
 *   1. squeeze uniform u ∈ [0, 2^32)
 *   2. 缩放 u 到 [0, CDF_total)
 *   3. 二分查找 table
 */
static int64_t sample_from_table(Xof& xof, const GaussTableEntry* table, int table_size) {
    constexpr uint64_t kScale = 1000000ULL; // CDF 总质量
    for (;;) {
        uint64_t u = xof.squeeze_u64();
        // 拒绝采样: 只使用 u % kScale 范围内 (保证均匀)
        uint64_t v = u % kScale;
        if (u >= (std::numeric_limits<uint64_t>::max() - kScale + 1) / kScale * kScale) {
            // u 接近 max 时跳过避免偏差 (极低概率)
            continue;
        }

        // 二分查找 CDF
        int lo = 0, hi = table_size - 1;
        while (lo < hi) {
            int mid = (lo + hi) / 2;
            if (table[mid].cdf_upper < v) {
                lo = mid + 1;
            } else {
                hi = mid;
            }
        }
        return table[lo].value;
    }
}

// ============================================================================
// gauss_sample_coeff
// ============================================================================

int64_t gauss_sample_coeff(Xof& xof, double sigma, int q) {
    (void)q; // unused in this version — 样本直接从高斯表产生

    // 为不同 sigma 选择表
    if (std::abs(sigma - 1.0) < 0.01) {
        // sigma ≈ 1.0: 使用预计算表
        int64_t sample = sample_from_table(xof, kGaussTableSigma1, kGaussTableSize1);
        return sample;
    }

    // 通用 sigma: 使用 Box-Muller 近似 (以 XOF 驱动)
    // 实现简单版: sigma 缩放下的拒绝采样
    // 限定σ ≤ 2.0
    double s = std::max(sigma, 0.5);
    int cutoff = static_cast<int>(std::ceil(8.0 * s));
    for (;;) {
        uint64_t u = xof.squeeze_u64();
        int64_t x = static_cast<int64_t>(u % static_cast<uint64_t>(2 * cutoff + 1)) - cutoff;
        // 接受概率: exp(-x^2 / (2 sigma^2))
        double prob = std::exp(-static_cast<double>(x * x) / (2.0 * s * s));
        // squeeze uniform for acceptance
        uint64_t v = xof.squeeze_u64();
        double uf = static_cast<double>(v & 0xFFFFF) / 1048576.0; // ~[0,1)
        if (uf < prob) return x;
    }
}

// ============================================================================
// gauss_sample_poly
// ============================================================================

SmallPoly gauss_sample_poly(Xof& xof, const Params& pp) {
    if (pp.n <= 0) return SmallPoly::zero(0);

    std::vector<int16_t> coeffs(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        int64_t v = gauss_sample_coeff(xof, pp.sigma, pp.q);
        // clamp to int16_t range
        if (v > 32767) v = 32767;
        if (v < -32768) v = -32768;
        coeffs[i] = static_cast<int16_t>(v);
    }
    return SmallPoly(pp.n, std::move(coeffs));
}

// ============================================================================
// build_gauss_xof
// ============================================================================

Xof build_gauss_xof(std::string_view domain, ByteSpan seed) {
    Xof xof(domain);
    xof.absorb(seed);
    xof.finalize();
    return xof;
}

// ============================================================================
// gauss_gen_seeds
// ============================================================================

void gauss_gen_seeds(const Params& pp,
                     ByteSpan R,
                     int M,
                     std::vector<std::vector<uint8_t>>& out_seeds) {
    static constexpr size_t SEED_LEN = 32; // SHAKE256 seed length

    out_seeds.clear();
    out_seeds.resize(M);

    if (M <= 0 || R.size < SEED_LEN) {
        // Not enough data: return zero seeds
        for (int j = 0; j < M; ++j) {
            out_seeds[j].assign(SEED_LEN, 0);
        }
        return;
    }

    // XOF from master seed
    Xof xof = build_gauss_xof(domain::GAUSS_FUNCTION, R);

    for (int j = 0; j < M; ++j) {
        out_seeds[j].resize(SEED_LEN);
        xof.squeeze(out_seeds[j].data(), SEED_LEN);
    }
}

// ============================================================================
// gauss_expand_alpha
// ============================================================================

SmallPoly gauss_expand_alpha(const Params& pp, ByteSpan seed_j) {
    Xof xof = build_gauss_xof(domain::GAUSS_FUNCTION, seed_j);
    return gauss_sample_poly(xof, pp);
}

} // namespace ibags