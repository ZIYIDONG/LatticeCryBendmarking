#pragma once
/**
 * eval.h — Eval(pp, C, (Ĉ₁, …, Ĉ_ℓ)) 同态加法 / 同态乘法
 *
 *   AddEval:   Ĉ⁺ = Ĉ₁ + Ĉ₂                       (mod q)
 *   MultEval:  Ĉˣ = Ĉ₁ · G⁻¹(Ĉ₂)                  (mod q)
 *
 *   G⁻¹ 把每个元素按 base-b 展开成 k = ⌈log_b q⌉ 个"位":
 *     G⁻¹: Z_q^{r × c}  →  {0,1,…,b-1}^{(r·k) × c}
 *     满足 G · G⁻¹(X) ≡ X   (mod q)
 *
 *   其中 G = I_r ⊗ (1, b, b², …, b^{k-1})  ∈ Z_q^{r × (r·k)}
 *
 *   尺寸兼容:
 *     mult_eval 要求 cols(C₁) = rows(C₂) · k,
 *     对满足 m = R·k 的 ABE/HIBE 参数自动成立.
 */

#include "matops.h"
#include <vector>
#include <stdexcept>

namespace cryptolib {

using matops::Mat;
using matops::Vec;

/* ─────────────────────────────────────────────
   §1 辅助: 计算 k = ⌈log_b q⌉
   ───────────────────────────────────────────── */
inline int eval_compute_k(long q, int b) {
    if (b < 2) throw std::invalid_argument("eval_compute_k: b must be >= 2");
    int k = 0;
    long p = 1;
    while (p < q) { p *= b; ++k; }
    return k;
}

/* ─────────────────────────────────────────────
   §2 gadget_inverse — G⁻¹ 运算
   输入  X ∈ Z_q^{r × c}
   输出  Y ∈ Z^{(r·k) × c}, 每个元素 ∈ {0,…,b-1}
   ───────────────────────────────────────────── */
inline Mat gadget_inverse(const Mat& X, long q, int b = 2) {
    if (X.empty() || X[0].empty())
        throw std::invalid_argument("gadget_inverse: empty matrix");
    const int k = eval_compute_k(q, b);
    const size_t r = X.size();
    const size_t c = X[0].size();
    Mat Y(r * k, Vec(c, 0));

    for (size_t i = 0; i < r; ++i) {
        for (size_t s = 0; s < c; ++s) {
            long v = X[i][s] % q;
            if (v < 0) v += q;
            for (int j = 0; j < k; ++j) {
                Y[i * k + j][s] = v % b;
                v /= b;
            }
        }
    }
    return Y;
}

/* ─────────────────────────────────────────────
   §3 build_gadget — 构造 G (用于测试和验证)
   G = I_r ⊗ (1, b, b², …, b^{k-1})
   形状: r × (r·k)
   ───────────────────────────────────────────── */
inline Mat build_gadget(size_t r, long q, int b = 2) {
    const int k = eval_compute_k(q, b);
    Mat G(r, Vec(r * k, 0));
    for (size_t i = 0; i < r; ++i) {
        long pw = 1;
        for (int j = 0; j < k; ++j) {
            G[i][i * k + j] = pw % q;
            pw *= b;
        }
    }
    return G;
}

/* ─────────────────────────────────────────────
   §4 AddEval — 同态加法
   Ĉ⁺ = Ĉ₁ + Ĉ₂  (mod q)
   ───────────────────────────────────────────── */
inline Mat add_eval(const Mat& C1, const Mat& C2, long q) {
    if (C1.size() != C2.size() || C1[0].size() != C2[0].size())
        throw std::invalid_argument("add_eval: shape mismatch");
    const size_t r = C1.size();
    const size_t c = C1[0].size();

    Mat out(r, Vec(c));
    for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < c; ++j) {
            long v = (C1[i][j] + C2[i][j]) % q;
            if (v < 0) v += q;
            out[i][j] = v;
        }
    return out;
}

/* ─────────────────────────────────────────────
   §5 MultEval — 同态乘法
   Ĉˣ = Ĉ₁ · G⁻¹(Ĉ₂)  (mod q)
   ───────────────────────────────────────────── */
inline Mat mult_eval(const Mat& C1, const Mat& C2, long q, int b = 2) {
    Mat Ginv = gadget_inverse(C2, q, b);     // (rows(C2)·k) × cols(C2)

    if (C1.empty() || C1[0].empty())
        throw std::invalid_argument("mult_eval: C1 empty");
    if (C1[0].size() != Ginv.size())
        throw std::invalid_argument(
            "mult_eval: cols(C1) != rows(G⁻¹(C2)); "
            "检查参数是否满足 m = R·k");

    const size_t r   = C1.size();
    const size_t mid = Ginv.size();
    const size_t c   = Ginv[0].size();

    Mat out(r, Vec(c, 0));
    // i-k-j 顺序, 缓存友好 (和 matops::mat_mul 保持风格一致)
    for (size_t i = 0; i < r; ++i) {
        for (size_t l = 0; l < mid; ++l) {
            long a = C1[i][l];
            if (a == 0) continue;
            const Vec& grow = Ginv[l];
            for (size_t j = 0; j < c; ++j)
                out[i][j] = (out[i][j] + a * grow[j]) % q;
        }
    }

    for (auto& row : out)
        for (auto& x : row) if (x < 0) x += q;

    return out;
}

} // namespace cryptolib
