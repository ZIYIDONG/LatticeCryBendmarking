#pragma once
/**
 * matops.h — 格密码常用矩阵操作(模 q)
 *
 * 实现:
 *   ① mat_add    : 矩阵加法 (A + B) mod q
 *   ② mat_sub    : 矩阵减法 (A - B) mod q
 *   ③ mat_mul    : 矩阵乘法 (A · B) mod q  (i-k-j 顺序,缓存友好)
 *   ④ mat_hcat   : 横向拼接 [A | B]
 *   ⑤ mat_vcat   : 纵向拼接 [A; B]      (附赠)
 *   ⑥ mat_eq     : 相等比较             (验证用)
 *
 * 设计原则:
 *   - 所有操作严格做维度检查(避免后期难以追踪的越界)
 *   - 模运算在最内层做,避免 long 溢出(假设 q < 2^31)
 *   - mat_mul 用 i-k-j 顺序 + 提前跳过零元素,对稀疏矩阵特别快
 *   - 与 mp12.h 兼容: 类型 Mat / Vec 一致,可以直接互通
 */

#include <vector>
#include <stdexcept>
#include <cstdint>
#include <string>
#include <random>

namespace matops {

using Vec = std::vector<long>;
using Mat = std::vector<Vec>;

inline long mod_pos(long x, long q) {
    return ((x % q) + q) % q;
}

inline Mat make_mat(int rows, int cols, long fill = 0) {
    return Mat(rows, Vec(cols, fill));
}

/* ══════════════════════════════════════════════════
   §0  辅助: 维度提取与检查
   ══════════════════════════════════════════════════ */
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

/* ══════════════════════════════════════════════════
   §1  矩阵加法 (mod q)
   ══════════════════════════════════════════════════ */
/**
 * C = (A + B) mod q
 * 复杂度: O(r·c)
 */
inline Mat mat_add(const Mat& A, const Mat& B, long q) {
    check_same_shape(A, B, "mat_add");
    auto [r, c] = dim(A);
    Mat C = make_mat(r, c);
    for (int i = 0; i < r; i++) {
        const Vec& ai = A[i];
        const Vec& bi = B[i];
        Vec& ci = C[i];
        for (int j = 0; j < c; j++)
            ci[j] = mod_pos(ai[j] + bi[j], q);
    }
    return C;
}

/* ══════════════════════════════════════════════════
   §2  矩阵减法 (mod q)
   ══════════════════════════════════════════════════ */
/**
 * C = (A - B) mod q
 * 复杂度: O(r·c)
 */
inline Mat mat_sub(const Mat& A, const Mat& B, long q) {
    check_same_shape(A, B, "mat_sub");
    auto [r, c] = dim(A);
    Mat C = make_mat(r, c);
    for (int i = 0; i < r; i++) {
        const Vec& ai = A[i];
        const Vec& bi = B[i];
        Vec& ci = C[i];
        for (int j = 0; j < c; j++)
            ci[j] = mod_pos(ai[j] - bi[j], q);
    }
    return C;
}

/* ══════════════════════════════════════════════════
   §3  矩阵乘法 (mod q)
   ══════════════════════════════════════════════════ */
/**
 * C = (A · B) mod q
 *
 * A: r × n,  B: n × c,  C: r × c
 * 复杂度: O(r·n·c)
 *
 * 优化:
 *   ① 使用 i-k-j 顺序而不是 i-j-k:
 *      内层循环按行连续访问 B[k] 和 C[i],缓存命中率高
 *   ② 当 A[i][k] == 0 时跳过整行,对稀疏矩阵(如 G、陷门 R)有效
 *   ③ 模运算放在最内层,避免长链累加溢出
 */
inline Mat mat_mul(const Mat& A, const Mat& B, long q) {
    auto [ra, ca] = dim(A);
    auto [rb, cb] = dim(B);
    if (ca != rb)
        throw std::invalid_argument("mat_mul: (" + std::to_string(ra) + "x" +
            std::to_string(ca) + ") · (" + std::to_string(rb) + "x" +
            std::to_string(cb) + ") incompatible");

    Mat C = make_mat(ra, cb);
    for (int i = 0; i < ra; i++) {
        Vec& ci = C[i];
        for (int k = 0; k < ca; k++) {
            long aik = A[i][k];
            if (aik == 0) continue;       // 跳过零行
            const Vec& bk = B[k];
            for (int j = 0; j < cb; j++)
                ci[j] = mod_pos(ci[j] + aik * bk[j], q);
        }
    }
    return C;
}

/* ══════════════════════════════════════════════════
   §4  横向拼接 [A | B]
   ══════════════════════════════════════════════════ */
/**
 * 横向拼接: 行数必须相同,列数相加
 * A: r × cA,  B: r × cB  →  [A|B]: r × (cA + cB)
 *
 * 应用: HIBE 中 A_{id_ℓ} = [A_{id_{ℓ-1}} || A_ℓ + H·G]
 * 复杂度: O(r · (cA + cB)),纯拷贝
 */
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

/* ══════════════════════════════════════════════════
   §5  纵向拼接 [A; B] (附赠 — 用于陷门 [R; I])
   ══════════════════════════════════════════════════ */
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

/* ══════════════════════════════════════════════════
   §6  相等比较 (验证用)
   ══════════════════════════════════════════════════ */
inline bool mat_eq(const Mat& A, const Mat& B, long q) {
    auto [ra, ca] = dim(A);
    auto [rb, cb] = dim(B);
    if (ra != rb || ca != cb) return false;
    for (int i = 0; i < ra; i++)
        for (int j = 0; j < ca; j++)
            if (mod_pos(A[i][j], q) != mod_pos(B[i][j], q)) return false;
    return true;
}

/* ══════════════════════════════════════════════════
   §7  随机矩阵 (基准测试用)
   ══════════════════════════════════════════════════ */
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
