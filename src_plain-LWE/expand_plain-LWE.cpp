#include "extend_plain-LWE.h"
#include "expand_plain-LWE.h"
#include "matops_plain-LWE.h"
#include "unified_params_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <vector>
#include <chrono>
#include <numeric>
#include <fstream>
#include <sstream>
#include "bench_utils_plain-LWE.h"

using namespace cryptolib;
using matops::Mat;
using matops::Vec;

/* ───── 前向声明 ───── */
static void bench_extend();
static void bench_expand_full();

void run_test_expand() {
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

    /* ── 纯基准测试 ── */
    bench_extend();
    bench_expand_full();

    std::cout << "\nDone.\n";
}

/* ────────────────────────────────────────────────────────────────────────────
   Benchmark 1: Extend (GSW.LComb) 纯耗时
   ─────────────────────────────────────────────────── */
static void bench_extend() {
    using Clock = std::chrono::high_resolution_clock;
    constexpr int WARMUP = 3;
    constexpr int ITERS  = 20;

    const size_t d = 3, N = d, m = 6;
    auto mp = unified::default_midparams_128((int)d, (int)N);
    const long q = mp.q;
    const size_t n = mp.n;
    const size_t R = (d + 1) * n + 1;

    std::mt19937 rng(42);
    std::uniform_int_distribution<long> uni(0, q - 1);
    std::bernoulli_distribution bit(0.5);

    /* ── 构造测试数据（不计入时间）── */
    UniEncU U(m, std::vector<Vec>(m, Vec(R, 0)));
    for (size_t r = 0; r < m; ++r)
        for (size_t s = 0; s < m; ++s)
            U[r][s][0] = bit(rng) ? 1 : 0;

    Vec b_i(m), b_j(m);
    for (size_t t = 0; t < m; ++t) { b_i[t] = uni(rng); b_j[t] = uni(rng); }

    /* ── 预热 ── */
    for (int k = 0; k < WARMUP; ++k) (void)extend(U, b_i, b_j, q);

    /* ── 计时 ── */
    std::vector<double> times_us; times_us.reserve(ITERS);
    for (int k = 0; k < ITERS; ++k) {
        b_i[0] = uni(rng); b_j[0] = uni(rng);  // 微调以防止编译器优化
        auto t0 = Clock::now();
        (void)extend(U, b_i, b_j, q);
        auto t1 = Clock::now();
        times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }

    double sum_us = std::accumulate(times_us.begin(), times_us.end(), 0.0);
    double avg_us = sum_us / ITERS;
    double min_us = *std::min_element(times_us.begin(), times_us.end());
    double max_us = *std::max_element(times_us.begin(), times_us.end());
    double var_us = 0;
    for (double t : times_us) { double d = t - avg_us; var_us += d * d; }
    var_us /= (ITERS > 1) ? (ITERS - 1) : 1;

    std::ostringstream oss;
    oss << "\n=== Benchmark: Extend (GSW.LComb) ===\n"
        << "  Parameters: d=" << d << ", N=" << N << ", n=" << n
        << ", q=" << q << ", R=" << R << ", m=" << m << "\n"
        << "  Warmup rounds : " << WARMUP << "\n"
        << "  Timed  rounds : " << ITERS << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << avg_us << " µs\n"
        << "  Min       : " << std::setw(8) << min_us << " µs\n"
        << "  Max       : " << std::setw(8) << max_us << " µs\n"
        << "  StdDev    : " << std::setw(8) << std::sqrt(var_us) << " µs\n"
        << "  Throughput: " << std::setw(8) << (1e6 / avg_us) << " ops/s\n";

    std::cout << "\n--- Benchmark: Extend ---\n" << oss.str();
    bench_write(oss.str());
}

/* ───────────────────────────────────────────────────
   Benchmark 2: Expand 完整操作 纯耗时
   ─────────────────────────────────────────────────── */
static void bench_expand_full() {
    using Clock = std::chrono::high_resolution_clock;
    constexpr int WARMUP = 3;
    constexpr int ITERS  = 20;

    const size_t d = 3, N = d, m = 6;
    auto mp = unified::default_midparams_128((int)d, (int)N);
    const long q = mp.q;
    const size_t n = mp.n;
    const size_t R = (d + 1) * n + 1;

    std::mt19937 rng(42);
    std::uniform_int_distribution<long> uni(0, q - 1);
    std::bernoulli_distribution bit(0.5);

    /* ── 构造测试数据（不计入时间）── */
    UniEncU U(m, std::vector<Vec>(m, Vec(R, 0)));
    for (size_t r = 0; r < m; ++r)
        for (size_t s = 0; s < m; ++s)
            U[r][s][0] = bit(rng) ? 1 : 0;

    std::vector<Vec> b_rows(N, Vec(m));
    for (size_t k = 0; k < N; ++k)
        for (size_t t = 0; t < m; ++t)
            b_rows[k][t] = uni(rng);

    Mat C(R, Vec(m));
    for (size_t r = 0; r < R; ++r)
        for (size_t c = 0; c < m; ++c)
            C[r][c] = uni(rng);

    /* ── 预热 ── */
    for (int k = 0; k < WARMUP; ++k) (void)expand(U, b_rows, 1, C, q);

    /* ── 计时 ── */
    std::vector<double> times_us; times_us.reserve(ITERS);
    for (int k = 0; k < ITERS; ++k) {
        C[0][0] = uni(rng);  // 微调防止优化
        auto t0 = Clock::now();
        (void)expand(U, b_rows, 1, C, q);
        auto t1 = Clock::now();
        times_us.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }

    double sum_us = std::accumulate(times_us.begin(), times_us.end(), 0.0);
    double avg_us = sum_us / ITERS;
    double min_us = *std::min_element(times_us.begin(), times_us.end());
    double max_us = *std::max_element(times_us.begin(), times_us.end());
    double var_us = 0;
    for (double t : times_us) { double d = t - avg_us; var_us += d * d; }
    var_us /= (ITERS > 1) ? (ITERS - 1) : 1;

    std::ostringstream oss;
    oss << "\n=== Benchmark: Expand (N×N block construction) ===\n"
        << "  Parameters: d=" << d << ", N=" << N << ", n=" << n
        << ", q=" << q << ", R=" << R << ", m=" << m << "\n"
        << "  Warmup rounds : " << WARMUP << "\n"
        << "  Timed  rounds : " << ITERS << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << avg_us << " µs\n"
        << "  Min       : " << std::setw(8) << min_us << " µs\n"
        << "  Max       : " << std::setw(8) << max_us << " µs\n"
        << "  StdDev    : " << std::setw(8) << std::sqrt(var_us) << " µs\n"
        << "  Throughput: " << std::setw(8) << (1e6 / avg_us) << " ops/s\n";

    std::cout << "\n--- Benchmark: Expand ---\n" << oss.str();
    bench_write(oss.str());
}

