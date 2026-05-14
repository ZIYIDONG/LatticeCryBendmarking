#pragma once
/**
 * expand.h — Expand((id_1,…,id_N), i, c) 密文扩展
 *
 *   Ĉ_i 是 N×N 分块矩阵,每个块 R×M (R=(d+1)n+1, M=m):
 *     Ĉ[a][a] = C                          (对角)
 *     Ĉ[i][j] = X_j = Extend(U, b_i, b_j)   (第 i 行非对角)
 *     Ĉ[a][b] = 0                          (其他)
 *
 *   输出尺寸: (N·R) × (N·M)
 */

#include "matops_plain-LWE.h"
#include "extend_plain-LWE.h"
#include <vector>
#include <stdexcept>

namespace cryptolib {

inline Mat expand(const UniEncU& U,
                  const std::vector<Vec>& b_rows,   // 每个身份的 b_k
                  size_t i,
                  const Mat& C,
                  long q)
{
    const size_t N = b_rows.size();
    if (N == 0)         throw std::invalid_argument("expand: empty ids");
    if (i >= N)         throw std::invalid_argument("expand: i out of range");
    if (C.empty() || C[0].empty())
                        throw std::invalid_argument("expand: C empty");

    const size_t R = C.size();        // (d+1)n+1
    const size_t M = C[0].size();     // m

    Mat Chat(N * R, Vec(N * M, 0));

    auto write_block = [&](size_t a, size_t b, const Mat& B) {
        if (B.size() != R || B[0].size() != M)
            throw std::runtime_error("expand: block size mismatch");
        for (size_t r = 0; r < R; ++r)
            for (size_t c = 0; c < M; ++c)
                Chat[a*R + r][b*M + c] = B[r][c];
    };

    // ① 对角线 = C  (注意 (i,i) 也走这个分支)
    for (size_t a = 0; a < N; ++a)
        write_block(a, a, C);

    // ② 第 i 行的非对角 = X_j
    for (size_t j = 0; j < N; ++j) {
        if (j == i) continue;
        Mat Xj = extend(U, b_rows[i], b_rows[j], q);
        write_block(i, j, Xj);
    }

    return Chat;
}

} // namespace cryptolib
