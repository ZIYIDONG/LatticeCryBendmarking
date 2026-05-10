#include "../include/extend.h"
#include "../include/expand.h"
#include "../include/matops.h"
#include "../include/unified_params.h"
#include <iostream>
#include <random>
#include <vector>

using namespace cryptolib;
using matops::Mat;
using matops::Vec;

int main() {
    /* ───── 参数 (来自 unified) ───── */
    const size_t d = 3;
    const size_t N = d;                  // 身份个数
    auto mp = unified::default_midparams_128((int)d, (int)N);
    const long q = mp.q;
    const size_t n = mp.n;
    const size_t R = (d + 1) * n + 1;
    const size_t m = 6;

    std::mt19937 rng(2025);
    std::uniform_int_distribution<long> uni(0, q - 1);
    std::bernoulli_distribution bit(0.5);

    /* ───── 1) 构造一个"明文"掩盖矩阵 R_mat ∈ {0,1}^{m×m} ─────
       并把 U[r][s] 直接设置成 R_mat[r][s] · e1 (长度 N 的单位列向量).
       这样 GSW.LComb 的输出就是一个"明文版本",可以直接代数验证. */
    std::vector<std::vector<long>> R_mat(m, std::vector<long>(m));
    for (size_t r = 0; r < m; ++r)
        for (size_t s = 0; s < m; ++s)
            R_mat[r][s] = bit(rng) ? 1 : 0;

    UniEncU U(m, std::vector<Vec>(m, Vec(R, 0)));
    for (size_t r = 0; r < m; ++r)
        for (size_t s = 0; s < m; ++s)
            U[r][s][0] = R_mat[r][s];     // 把比特放在第 0 个位置

    /* ───── 2) 构造 N 个身份的 b_k ───── */
    std::vector<Vec> b_rows(N, Vec(m));
    for (size_t k = 0; k < N; ++k)
        for (size_t t = 0; t < m; ++t)
            b_rows[k][t] = uni(rng);

    /* ───── 3) 单独验证 extend / GSW.LComb ───── */
    size_t i = 1, j = 2;
    Mat Xj = extend(U, b_rows[i], b_rows[j], q);

    // 期望: Xj[0][s] = Σ_r (b_j[r]-b_i[r]) · R_mat[r][s]   (mod q)
    bool extend_ok = true;
    for (size_t s = 0; s < m; ++s) {
        long expected = 0;
        for (size_t r = 0; r < m; ++r) {
            long diff = (b_rows[j][r] - b_rows[i][r]) % q;
            if (diff < 0) diff += q;
            expected = (expected + diff * R_mat[r][s]) % q;
        }
        if (Xj[0][s] != expected) extend_ok = false;
    }
    std::cout << "Test 1 — extend 代数正确性: "
              << (extend_ok ? "PASS" : "FAIL") << "\n";

    /* ───── 4) 构造一个假密文 C, 测试 expand 的分块结构 ───── */
    Mat C(R, Vec(m));
    for (auto& row : C) for (auto& x : row) x = uni(rng);

    Mat Chat = expand(U, b_rows, i, C, q);

    std::cout << "  Ĉ_" << i << " 尺寸: "
              << Chat.size() << " × " << Chat[0].size()
              << "  (期望 " << N*R << " × " << N*m << ")\n";

    // 检查对角块 = C
    bool diag_ok = true;
    for (size_t a = 0; a < N; ++a)
        for (size_t r = 0; r < R; ++r)
            for (size_t c = 0; c < m; ++c)
                if (Chat[a*R+r][a*m+c] != C[r][c]) diag_ok = false;
    std::cout << "Test 2 — 对角块 = C: "
              << (diag_ok ? "PASS" : "FAIL") << "\n";

    // 检查第 i 行非对角 = extend 的输出
    bool row_ok = true;
    for (size_t jj = 0; jj < N; ++jj) {
        if (jj == i) continue;
        Mat expected = extend(U, b_rows[i], b_rows[jj], q);
        for (size_t r = 0; r < R; ++r)
            for (size_t c = 0; c < m; ++c)
                if (Chat[i*R+r][jj*m+c] != expected[r][c]) row_ok = false;
    }
    std::cout << "Test 3 — 第 " << i << " 行非对角 = X_j: "
              << (row_ok ? "PASS" : "FAIL") << "\n";

    // 检查其他位置 = 0
    bool zero_ok = true;
    for (size_t a = 0; a < N; ++a) {
        if (a == i) continue;
        for (size_t b = 0; b < N; ++b) {
            if (a == b) continue;
            for (size_t r = 0; r < R; ++r)
                for (size_t c = 0; c < m; ++c)
                    if (Chat[a*R+r][b*m+c] != 0) zero_ok = false;
        }
    }
    std::cout << "Test 4 — 其他位置 = 0: "
              << (zero_ok ? "PASS" : "FAIL") << "\n";

    std::cout << "\nDone.\n";
    return 0;
}
