#pragma once
/**
 * matops.h — 格密码常用矩阵操作(模 q)
 *
 * 实现:
 *   ① mat_add      : 矩阵加法 (A + B) mod q   — 条件减法
 *   ② mat_sub      : 矩阵减法 (A - B) mod q   — 条件加法
 *   ③ mat_mul      : 矩阵乘法 (A · B) mod q   — Barrett 内循环
 *   ④ mat_hcat     : 横向拼接 [A | B]
 *   ⑤ mat_vcat     : 纵向拼接 [A; B]
 *   ⑥ mat_eq       : 相等比较 (验证用)
 *   ⑦ mat_diag_mul : 对角矩阵 × 矩阵 (左)     — NEW
 *   ⑧ vec_diag_mul : 对角矩阵 × 向量          — NEW
 *
 * 设计原则:
 *   - 所有操作严格做维度检查
 *   - mat_mul 内循环使用 Barret 乘-移位替代 % q 硬除法
 *   - mat_add / mat_sub 使用条件加减 (输入范围已知)
 *   - mat_mul 用 i-k-j 顺序 + 提前跳过零元素
 */

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <random>
#include <cassert>

namespace matops {

using Vec = std::vector<long>;
using Mat = std::vector<Vec>;

// ============================================================================
// Barrett 约减 (plain-LWE 侧, 内联零开销)
// ============================================================================

/// Barrett 预计算常数: mu = floor(2^64 / q)
inline unsigned long long barrett_mu(long q) {
    assert(q > 0 && q < (1L << 30));
    return (unsigned long long)(((unsigned __int128)1 << 64) / (unsigned long long)q);
}

/// Barrett 约减: x mod q → [0, q)
/// 正确性: |x| < q * 2^32 (使用 |x| + 符号恢复, 正确处理 INT64_MIN)
inline long barrett_reduce_lwe(long x, long q, unsigned long long mu) {
    bool neg = (x < 0);
    unsigned long long abs_x;
    if (neg) abs_x = (unsigned long long)(-(x + 1)) + 1;  // safe |x|
    else     abs_x = (unsigned long long)x;
    unsigned long long t = (unsigned long long)(
        ((unsigned __int128)abs_x * mu) >> 64);
    unsigned long long r = abs_x - t * (unsigned long long)q;
    if (r >= (unsigned long long)q) r -= q;
    if (neg && r != 0) r = (unsigned long long)q - r;
    return (long)r;
}

// ============================================================================
// 基础约减
// ============================================================================

/// 通用 mod_pos (向后兼容, 非热点路径)
inline long mod_pos(long x, long q) {
    long r = x % q;
    if (r < 0) r += q;
    return r;
}

/// 加法后约减: 输入和 ≤ 2q-2, 至多一次条件减法
inline long mod_add(long x, long q) {
    if (x >= q) x -= q;
    return x;
}

/// 减法后约减: 输入差 ≥ -(q-1), 至多一次条件加法
inline long mod_sub(long x, long q) {
    if (x < 0) x += q;
    return x;
}

// ============================================================================
// 辅助
// ============================================================================

inline Mat make_mat(int rows, int cols, long fill = 0) {
    return Mat(rows, Vec(cols, fill));
}

inline std::pair<int,int> dim(const Mat& A) {
    if (A.empty()) return {0, 0};
    return {(int)A.size(), (int)A[0].size()};
}

inline void check_same_shape(const Mat& A, const Mat& B, const char* op) {
    auto [ra, ca] = dim(A);
    auto [rb, cb] = dim(B);
    if (ra != rb || ca != cb)
        throw std::invalid_argument(std::string(op) +
            ": shape mismatch (" + std::to_string(ra) + "x" + std::to_string(ca) +
            " vs " + std::to_string(rb) + "x" + std::to_string(cb) + ")");
}

// ============================================================================
// §1  矩阵加法 (mod q) — 条件减法
// ============================================================================
inline Mat mat_add(const Mat& A, const Mat& B, long q) {
    check_same_shape(A, B, "mat_add");
    auto [r, c] = dim(A);
    Mat C = make_mat(r, c);
    for (int i = 0; i < r; i++) {
        const Vec& ai = A[i];
        const Vec& bi = B[i];
        Vec& ci = C[i];
        for (int j = 0; j < c; j++) {
            long s = ai[j] + bi[j];
            ci[j] = mod_add(s, q);
        }
    }
    return C;
}

// ============================================================================
// §2  矩阵减法 (mod q) — 条件加法
// ============================================================================
inline Mat mat_sub(const Mat& A, const Mat& B, long q) {
    check_same_shape(A, B, "mat_sub");
    auto [r, c] = dim(A);
    Mat C = make_mat(r, c);
    for (int i = 0; i < r; i++) {
        const Vec& ai = A[i];
        const Vec& bi = B[i];
        Vec& ci = C[i];
        for (int j = 0; j < c; j++) {
            long d = ai[j] - bi[j];
            ci[j] = mod_sub(d, q);
        }
    }
    return C;
}

// ============================================================================
// §3  矩阵乘法 (mod q) — Barrett 内循环
// ============================================================================
inline Mat mat_mul(const Mat& A, const Mat& B, long q) {
    auto [ra, ca] = dim(A);
    auto [rb, cb] = dim(B);
    if (ca != rb)
        throw std::invalid_argument("mat_mul: (" + std::to_string(ra) + "x" +
            std::to_string(ca) + ") · (" + std::to_string(rb) + "x" +
            std::to_string(cb) + ") incompatible");

    Mat C = make_mat(ra, cb);
    unsigned long long mu = barrett_mu(q);
    for (int i = 0; i < ra; i++) {
        Vec& ci = C[i];
        for (int k = 0; k < ca; k++) {
            long aik = A[i][k];
            if (aik == 0) continue;
            const Vec& bk = B[k];
            for (int j = 0; j < cb; j++)
                ci[j] = barrett_reduce_lwe(ci[j] + aik * bk[j], q, mu);
        }
    }
    return C;
}

// ============================================================================
// §4  横向拼接 [A | B]
// ============================================================================
inline Mat mat_hcat(const Mat& A, const Mat& B) {
    auto [rA, cA] = dim(A);
    auto [rB, cB] = dim(B);
    if (rA != rB)
        throw std::invalid_argument("mat_hcat: row count mismatch ("
            + std::to_string(rA) + " vs " + std::to_string(rB) + ")");
    Mat C = make_mat(rA, cA + cB);
    for (int i = 0; i < rA; i++) {
        for (int j = 0; j < cA; j++) C[i][j]      = A[i][j];
        for (int j = 0; j < cB; j++) C[i][cA + j] = B[i][j];
    }
    return C;
}

// ============================================================================
// §5  纵向拼接 [A; B]
// ============================================================================
inline Mat mat_vcat(const Mat& A, const Mat& B) {
    auto [rA, cA] = dim(A);
    auto [rB, cB] = dim(B);
    if (cA != cB)
        throw std::invalid_argument("mat_vcat: col count mismatch ("
            + std::to_string(cA) + " vs " + std::to_string(cB) + ")");
    Mat C = make_mat(rA + rB, cA);
    for (int i = 0; i < rA; i++) C[i] = A[i];
    for (int i = 0; i < rB; i++) C[rA + i] = B[i];
    return C;
}

// ============================================================================
// §6  相等比较 (验证用)
// ============================================================================
inline bool mat_eq(const Mat& A, const Mat& B, long q) {
    auto [ra, ca] = dim(A);
    auto [rb, cb] = dim(B);
    if (ra != rb || ca != cb) return false;
    for (int i = 0; i < ra; i++)
        for (int j = 0; j < ca; j++)
            if (mod_pos(A[i][j], q) != mod_pos(B[i][j], q)) return false;
    return true;
}

// ============================================================================
// §7  对角矩阵 × 矩阵 (左乘): C = diag · B
//     即 C[i][j] = diag[i] * B[i][j] mod q,  O(r·c)
// ============================================================================
inline Mat mat_diag_mul(const Vec& diag, const Mat& B, long q) {
    auto [r, c] = dim(B);
    if ((int)diag.size() != r)
        throw std::invalid_argument("mat_diag_mul: diag size " +
            std::to_string(diag.size()) + " != B rows " + std::to_string(r));
    Mat C = make_mat(r, c);
    unsigned long long mu = barrett_mu(q);
    for (int i = 0; i < r; i++) {
        long d = diag[i];
        if (d == 0) continue;               // 整行为零, C[i] 已初始化为 0
        const Vec& bi = B[i];
        Vec& ci = C[i];
        for (int j = 0; j < c; j++)
            ci[j] = barrett_reduce_lwe(d * bi[j], q, mu);
    }
    return C;
}

// ============================================================================
// §8  对角矩阵 × 向量: result[i] = diag[i] * v[i] mod q,  O(n)
// ============================================================================
inline Vec vec_diag_mul(const Vec& diag, const Vec& v, long q) {
    int n = (int)v.size();
    if ((int)diag.size() != n)
        throw std::invalid_argument("vec_diag_mul: diag size " +
            std::to_string(diag.size()) + " != v size " + std::to_string(n));
    Vec result(n, 0);
    unsigned long long mu = barrett_mu(q);
    for (int i = 0; i < n; i++) {
        result[i] = barrett_reduce_lwe(diag[i] * v[i], q, mu);
    }
    return result;
}

// ============================================================================
// §9  随机矩阵 (基准测试用)
// ============================================================================
inline Mat random_mat(int r, int c, long q, uint64_t seed = 0) {
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::uniform_int_distribution<long> d(0, q - 1);
    Mat M = make_mat(r, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            M[i][j] = d(rng);
    return M;
}

} // namespace matops
