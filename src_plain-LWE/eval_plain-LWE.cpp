#include "eval_plain-LWE.h"
#include "matops_plain-LWE.h"
#include "unified_params_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <numeric>
#include <fstream>
#include <sstream>
#include "bench_utils_plain-LWE.h"

using namespace cryptolib;
using matops::Mat;
using matops::Vec;

/* ───── 前向声明 ───── */
static void bench_add_eval();
static void bench_mult_eval();

// 用 matops 的 mat_mul 做验证参考 (直接写一个小版本, 避免依赖其命名空间差异)
static Mat naive_mul(const Mat& A, const Mat& B, long q) {
    size_t r = A.size(), mid = A[0].size(), c = B[0].size();
    Mat out(r, Vec(c, 0));
    for (size_t i = 0; i < r; ++i)
        for (size_t l = 0; l < mid; ++l) {
            long a = A[i][l];
            if (!a) continue;
            for (size_t j = 0; j < c; ++j)
                out[i][j] = (out[i][j] + a * B[l][j]) % q;
        }
    for (auto& row : out) for (auto& x : row) if (x < 0) x += q;
    return out;
}

static Mat scalar_mul(long s, const Mat& A, long q) {
    Mat out = A;
    for (auto& row : out)
        for (auto& x : row) {
            long v = (s * x) % q;
            if (v < 0) v += q;
            x = v;
        }
    return out;
}

void run_test_eval() {
    auto mp = unified::default_mp12_params();
    const long q = mp.q;
    const int  b = mp.b;
    const int  k = eval_compute_k(q, b);
    const size_t r = 6;            // "行维度"
    const size_t c = r * k;        // 必须 = r·k 才能乘

    std::cout << "q=" << q << "  b=" << b << "  k=" << k
              << "  r=" << r << "  c=" << c << "\n\n";

    std::mt19937 rng(123);
    std::uniform_int_distribution<long> uni(0, q - 1);

    // 构造 gadget G
    Mat G = build_gadget(r, q, b);

    /* ───── Test 1: G · G⁻¹(X) = X ───── */
    {
        Mat X(r, Vec(c));
        for (auto& row : X) for (auto& x : row) x = uni(rng);

        Mat Ginv = gadget_inverse(X, q, b);
        Mat reconstructed = naive_mul(G, Ginv, q);

        bool ok = (reconstructed == X);
        std::cout << "Test 1 — G · G⁻¹(X) = X: "
                  << (ok ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 2: AddEval ───── */
    {
        Mat C1(r, Vec(c)), C2(r, Vec(c));
        for (auto& row : C1) for (auto& x : row) x = uni(rng);
        for (auto& row : C2) for (auto& x : row) x = uni(rng);

        Mat sum = add_eval(C1, C2, q);
        bool ok = true;
        for (size_t i = 0; i < r; ++i)
            for (size_t j = 0; j < c; ++j) {
                long exp = (C1[i][j] + C2[i][j]) % q;
                if (exp < 0) exp += q;
                if (sum[i][j] != exp) ok = false;
            }
        std::cout << "Test 2 — AddEval 逐元素正确: "
                  << (ok ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 3: MultEval 代数正确性 ─────
       设 C₁ = μ₁·G, C₂ = μ₂·G
       期望 mult_eval(C₁,C₂) = μ₁μ₂·G        */
    {
        long mu1 = 5, mu2 = 7;
        Mat C1 = scalar_mul(mu1, G, q);
        Mat C2 = scalar_mul(mu2, G, q);

        Mat result   = mult_eval(C1, C2, q, b);
        Mat expected = scalar_mul((mu1 * mu2) % q, G, q);

        bool ok = (result == expected);
        std::cout << "Test 3 — MultEval (μ₁G)·G⁻¹(μ₂G) = μ₁μ₂·G: "
                  << (ok ? "PASS" : "FAIL")
                  << "  (μ₁=" << mu1 << ", μ₂=" << mu2 << ")\n";
    }

    /* ───── Test 4: 多次同态乘法链 ───── */
    {
        long mu = 3;
        Mat C = scalar_mul(mu, G, q);

        // 连乘: C² = C·G⁻¹(C),  C³ = C²·G⁻¹(C)
        Mat C2 = mult_eval(C, C, q, b);
        Mat C3 = mult_eval(C2, C, q, b);

        Mat expected2 = scalar_mul((mu * mu) % q, G, q);
        Mat expected3 = scalar_mul((mu * mu * mu) % q, G, q);

        bool ok2 = (C2 == expected2);
        bool ok3 = (C3 == expected3);
        std::cout << "Test 4 — MultEval 连乘 C² / C³: "
                  << (ok2 && ok3 ? "PASS" : "FAIL")
                  << "  (μ=" << mu << ", 期望 μ²=" << (mu*mu) % q
                  << ", μ³=" << (mu*mu*mu) % q << ")\n";
    }

    /* ── 纯基准测试 ── */
    bench_add_eval();
    bench_mult_eval();

    std::cout << "\nDone.\n";
}

/* ────────────────────────────────────────────────────────────────────────────
   Benchmark 1: AddEval 纯耗时
   ─────────────────────────────────────────────────── */
static void bench_add_eval() {
    auto mp = unified::default_mp12_params();
    const long q = mp.q;
    const int  b = mp.b;
    const int  k = eval_compute_k(q, b);
    const size_t r = 6;
    const size_t c = r * k;

    std::mt19937 rng(42);
    std::uniform_int_distribution<long> uni(0, q - 1);

    Mat C1(r, Vec(c)), C2(r, Vec(c));
    for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < c; ++j) {
            C1[i][j] = uni(rng);
            C2[i][j] = uni(rng);
        }

    std::cout << "  Benchmarking AddEval (1000 rounds, 50 warmup)... " << std::flush;
    auto stats = run_benchmark([&]() {
        C1[0][0] = uni(rng);
        (void)add_eval(C1, C2, q);
    }, 1000, 50);
    std::cout << "done\n" << std::flush;

    auto report = [&](std::ostream& os) {
        os << "\n=== Benchmark: AddEval (GSW homomorphic addition) ===\n"
           << "  Parameters: q=" << q << ", b=" << b
           << ", k=" << k << ", r=" << r << ", c=" << c << "\n"
           << "  Total iterations : 1000\n"
           << "  Warmup discarded : 50\n"
           << "  Active samples   : " << stats.active_samples << "\n\n"
           << std::fixed << std::setprecision(1)
           << "  Average   : " << std::setw(8) << stats.avg_us << " us\n"
           << "  Median    : " << std::setw(8) << stats.median_us << " us\n"
           << "  StdDev    : " << std::setw(8) << stats.stddev_us << " us\n"
           << "  Min       : " << std::setw(8) << stats.min_us << " us\n"
           << "  Max       : " << std::setw(8) << stats.max_us << " us\n"
           << "  Throughput: " << std::setw(8)
           << (1e6 / stats.avg_us) << " ops/s\n";
    };

    std::cout << "\n--- Benchmark: AddEval ---\n";
    report(std::cout);

    std::ostringstream oss;
    report(oss);
    bench_write(oss.str());
}

/* ───────────────────────────────────────────────────
   Benchmark 2: MultEval 纯耗时
   ─────────────────────────────────────────────────── */
static void bench_mult_eval() {
    auto mp = unified::default_mp12_params();
    const long q = mp.q;
    const int  b = mp.b;
    const int  k = eval_compute_k(q, b);
    const size_t r = 6;
    const size_t c = r * k;

    std::mt19937 rng(42);
    std::uniform_int_distribution<long> uni(0, q - 1);

    Mat G = build_gadget(r, q, b);
    Mat C1(r, Vec(c)), C2(r, Vec(c));
    for (size_t i = 0; i < r; ++i)
        for (size_t j = 0; j < c; ++j) {
            C1[i][j] = uni(rng);
            C2[i][j] = uni(rng);
        }

    std::cout << "  Benchmarking MultEval (1000 rounds, 50 warmup)... " << std::flush;
    auto stats = run_benchmark([&]() {
        C1[0][0] = uni(rng);
        (void)mult_eval(C1, C2, q, b);
    }, 1000, 50);
    std::cout << "done\n" << std::flush;

    auto report = [&](std::ostream& os) {
        os << "\n=== Benchmark: MultEval (GSW homomorphic multiplication) ===\n"
           << "  Parameters: q=" << q << ", b=" << b
           << ", k=" << k << ", r=" << r << ", c=" << c << "\n"
           << "  Total iterations : 1000\n"
           << "  Warmup discarded : 50\n"
           << "  Active samples   : " << stats.active_samples << "\n\n"
           << std::fixed << std::setprecision(1)
           << "  Average   : " << std::setw(8) << stats.avg_us << " us\n"
           << "  Median    : " << std::setw(8) << stats.median_us << " us\n"
           << "  StdDev    : " << std::setw(8) << stats.stddev_us << " us\n"
           << "  Min       : " << std::setw(8) << stats.min_us << " us\n"
           << "  Max       : " << std::setw(8) << stats.max_us << " us\n"
           << "  Throughput: " << std::setw(8)
           << (1e6 / stats.avg_us) << " ops/s\n";
    };

    std::cout << "\n--- Benchmark: MultEval ---\n";
    report(std::cout);

    std::ostringstream oss;
    report(oss);
    bench_write(oss.str());
}
