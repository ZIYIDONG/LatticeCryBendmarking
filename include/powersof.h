#pragma once
/**
 * Powers-of-b  与  BitDecomp_b  实现
 *
 * 数学背景:
 *   设 q 为模数, b 为基(通常 2), k = ceil(log_b q)
 *   gadget 向量 g = (1, b, b^2, ..., b^{k-1})
 *
 *   Powersof_b(a)     = a · g  =  (a, ba, b²a, ..., b^{k-1}a)  mod q
 *   BitDecomp_b(a)    = (a_0, a_1, ..., a_{k-1})  s.t.  Σ a_i · b^i = a (mod q)
 *
 *   核心对偶恒等式 (GSW / 密钥切换):
 *     <BitDecomp_b(a), Powersof_b(b)> = a · b  (mod q)
 *
 * 用途:
 *   - GSW 同态加密的密钥切换
 *   - BGV/BFV 的 relinearization key
 *   - MP12 gadget trapdoor 的对偶视角
 *   - GGH15 graded encoding
 */

#include <vector>
#include <cstdint>
#include <cmath>
#include <stdexcept>

namespace cryptolib {

using Vec = std::vector<long>;
using Mat = std::vector<Vec>;

inline long mod_q(long x, long q) {
    return ((x % q) + q) % q;
}

/* ══════════════════════════════════════════════════
   §1  辅助: 计算 k = ceil(log_b q)
   ══════════════════════════════════════════════════ */
inline int compute_k(long q, int b) {
    if (q <= 1 || b < 2)
        throw std::invalid_argument("Need q > 1 and b >= 2");
    int k = 0;
    long pw = 1;
    while (pw < q) { pw *= b; k++; }
    return k;
}

/* ══════════════════════════════════════════════════
   §2  Powersof_b  ——  标量版本
   ══════════════════════════════════════════════════ */
/**
 * Powersof_b(a) = (a, ba, b²a, ..., b^{k-1}a)  mod q
 *
 * 输入: a ∈ Z_q, 基 b, 模数 q
 * 输出: 长度 k = ceil(log_b q) 的向量
 *
 * 复杂度: O(k)
 */
inline Vec powers_of_b_scalar(long a, int b, long q) {
    int k = compute_k(q, b);
    Vec result(k);
    long val = mod_q(a, q);
    for (int j = 0; j < k; j++) {
        result[j] = val;
        val = mod_q(val * b, q);   // 下一项 = 当前项 · b
    }
    return result;
}

/* ══════════════════════════════════════════════════
   §3  Powersof_b  ——  向量版本
   ══════════════════════════════════════════════════ */
/**
 * Powersof_b(v) ∈ Z_q^{nk}
 * 把每个 v_i 单独展开为 (v_i, b·v_i, ..., b^{k-1}·v_i)
 * 然后顺序拼接
 *
 * 这正是 MP12 中向量 v 与 gadget 矩阵 G^T 相乘的结果:
 *   Powersof_b(v) = G^T · v   (其中 G = I_n ⊗ g^T)
 */
inline Vec powers_of_b_vec(const Vec& v, int b, long q) {
    int n = (int)v.size();
    int k = compute_k(q, b);
    Vec result(n * k);
    for (int i = 0; i < n; i++) {
        long val = mod_q(v[i], q);
        for (int j = 0; j < k; j++) {
            result[i * k + j] = val;
            val = mod_q(val * b, q);
        }
    }
    return result;
}

/* ══════════════════════════════════════════════════
   §4  Powersof_b  ——  矩阵版本(逐列展开)
   ══════════════════════════════════════════════════ */
/**
 * 对矩阵 M ∈ Z_q^{n×m} 的每一列应用 Powersof_b
 * 输出维度: nk × m
 *
 * 用于密钥切换:把密钥矩阵 sk 展开成长形式,
 * 与 BitDecomp 后的密文相乘可恢复原内积
 */
inline Mat powers_of_b_mat(const Mat& M, int b, long q) {
    if (M.empty()) return {};
    int n = (int)M.size();
    int m = (int)M[0].size();
    int k = compute_k(q, b);
    Mat result(n * k, Vec(m, 0));
    for (int j = 0; j < m; j++) {
        for (int i = 0; i < n; i++) {
            long val = mod_q(M[i][j], q);
            for (int e = 0; e < k; e++) {
                result[i * k + e][j] = val;
                val = mod_q(val * b, q);
            }
        }
    }
    return result;
}

/* ══════════════════════════════════════════════════
   §5  BitDecomp_b  ——  对偶分解
   ══════════════════════════════════════════════════ */
/**
 * BitDecomp_b(a) = (a_0, a_1, ..., a_{k-1})
 *   满足  Σ_{j=0}^{k-1} a_j · b^j  ≡  a  (mod q)
 *   且每个 a_j ∈ [0, b)  (非平衡版本)
 *
 * 这是 Powersof_b 的对偶,合在一起满足:
 *   <BitDecomp_b(x), Powersof_b(y)> = x·y  (mod q)
 */
inline Vec bit_decomp_scalar(long a, int b, long q) {
    int k = compute_k(q, b);
    Vec digits(k);
    long val = mod_q(a, q);
    for (int j = 0; j < k; j++) {
        digits[j] = val % b;
        val /= b;
    }
    return digits;
}

/**
 * 向量版 BitDecomp_b: 对 v ∈ Z_q^n 的每个分量分解,顺序拼接
 * 输出维度: nk
 */
inline Vec bit_decomp_vec(const Vec& v, int b, long q) {
    int n = (int)v.size();
    int k = compute_k(q, b);
    Vec result(n * k);
    for (int i = 0; i < n; i++) {
        long val = mod_q(v[i], q);
        for (int j = 0; j < k; j++) {
            result[i * k + j] = val % b;
            val /= b;
        }
    }
    return result;
}

/* ══════════════════════════════════════════════════
   §6  平衡 BitDecomp  ——  数字落在 [-b/2, b/2)
   ══════════════════════════════════════════════════ */
/**
 * 平衡分解: 每个数字 ∈ [-b/2, b/2),范数更小,常用于陷门采样
 * 这与 mp12.h 中的 sample_g 是同一个分解
 */
inline Vec bit_decomp_balanced(long a, int b, long q) {
    int k = compute_k(q, b);
    Vec digits(k);
    long val = mod_q(a, q);
    for (int j = 0; j < k; j++) {
        long d = val % b;
        if (d > b / 2) d -= b;          // 中心化到 [-b/2, b/2)
        digits[j] = d;
        val = (val - d) / b;
    }
    return digits;
}

/* ══════════════════════════════════════════════════
   §7  内积验证: <BitDecomp(x), Powersof(y)> == x·y
   ══════════════════════════════════════════════════ */
inline long inner_prod_mod(const Vec& a, const Vec& b, long q) {
    if (a.size() != b.size())
        throw std::invalid_argument("size mismatch");
    long s = 0;
    for (size_t i = 0; i < a.size(); i++)
        s = mod_q(s + a[i] * b[i], q);
    return s;
}

} // namespace cryptolib
