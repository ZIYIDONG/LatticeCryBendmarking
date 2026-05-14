#pragma once
/**
 * extend.h — Extend(pp, U, id_i, id_j) → X_j
 *
 *   X_j = GSW.LComb(U, b_j - b_i)
 *
 *   其中:
 *     U[r][s] = V^{(r+1,s+1)} ∈ Z_q^N    (m×m 个 LWE 密文,N=(d+1)n+1)
 *     b_k 是 A_{id_k} = [B_k; b_k] 的下半部分,长度 m
 *     输出 X_j ∈ Z_q^{N×m}
 *
 *   语义: X_j 是  (b_j - b_i)^T · R  的同态密文,其中 R 是 UniEnc 用的掩盖矩阵
 */

#include "matops_plain-LWE.h"
#include <vector>
#include <stdexcept>

namespace cryptolib {

using matops::Mat;
using matops::Vec;

// UniEnc 输出的 m×m 个 LWE 密文,U[r][s] 是一个 N 维列向量
using UniEncU = std::vector<std::vector<Vec>>;

/* ═══════════════════════════════════════════════
   GSW.LComb — 同态线性组合
   X[:, s] = Σ_r  coeffs[r] · U[r][s]    (mod q)
   输出 X ∈ Z_q^{N × m}
   ═══════════════════════════════════════════════ */
inline Mat gsw_lcomb(const UniEncU& U, const Vec& coeffs, long q)
{
    const size_t m = coeffs.size();
    if (U.size() != m)
        throw std::invalid_argument("gsw_lcomb: |U| rows != m");
    for (const auto& row : U)
        if (row.size() != m)
            throw std::invalid_argument("gsw_lcomb: U not square m×m");
    if (U[0][0].empty())
        throw std::invalid_argument("gsw_lcomb: V is empty vector");

    const size_t N = U[0][0].size();
    Mat X(N, Vec(m, 0));

    for (size_t s = 0; s < m; ++s) {
        for (size_t r = 0; r < m; ++r) {
            long c = coeffs[r];
            if (c == 0) continue;                     // 跳过零, 显著加速
            const Vec& V = U[r][s];
            if (V.size() != N)
                throw std::runtime_error("gsw_lcomb: V[r][s] dim mismatch");
            for (size_t t = 0; t < N; ++t)
                X[t][s] = (X[t][s] + c * V[t]) % q;
        }
    }

    // 规范化到 [0, q)
    for (auto& row : X)
        for (auto& x : row) if (x < 0) x += q;

    return X;
}

/* ═══════════════════════════════════════════════
   Extend(U, b_i, b_j) → X_j
   ═══════════════════════════════════════════════ */
inline Mat extend(const UniEncU& U,
                  const Vec& b_i,
                  const Vec& b_j,
                  long q)
{
    if (b_i.size() != b_j.size())
        throw std::invalid_argument("extend: |b_i| != |b_j|");

    // diff = b_j - b_i  (mod q)
    const size_t m = b_i.size();
    Vec diff(m);
    for (size_t k = 0; k < m; ++k) {
        long d = (b_j[k] - b_i[k]) % q;
        if (d < 0) d += q;
        diff[k] = d;
    }

    return gsw_lcomb(U, diff, q);
}

} // namespace cryptolib
