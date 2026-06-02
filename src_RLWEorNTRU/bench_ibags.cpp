/**
 * bench_ibags.cpp — IBAGS (RLWE/NTRU) benchmark suite
 *
 * Measures all fundamental cryptographic operations in the IBAGS signature
 * submodule and writes results to benchmarking_output/benchmarking_ibags.txt.
 *
 * Build: 由顶层 CMakeLists.txt 管理，通过 make 构建
 */

#include "params.h"
#include "poly.h"
#include "poly_ntt.h"
#include "mod_reduce.h"
#include "gauss.h"
#include "xof.h"
#include "csprng.h"
#include "encode.h"
#include "hash.h"
#include "reject.h"
#include "ntru_secret.h"
#include "ntru_trapgen.h"
#include "domain.h"
#include "errors.h"
#include "bench_utils.h"

#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <ctime>
#include <vector>
#include <string>
#include <numeric>
#include <cmath>
#include <functional>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace ibags;
using Clock = std::chrono::high_resolution_clock;

/* =========================================================================
   全局参数
   ========================================================================= */
static Params    g_demo_params;
static Params    g_ntt_params;
static NttTable  g_ntt_table;  // 在 init_params() 中赋值

static void init_params() {
    g_demo_params = Params::params_demo_64();

    g_ntt_params.n                   = 64;
    g_ntt_params.q                   = 7681;  // NTT-friendly: 7681 ≡ 1 mod 128
    g_ntt_params.sigma               = 1.0;
    g_ntt_params.eta1                = 4;
    g_ntt_params.eta2                = 2;
    g_ntt_params.eta_ver             = 2;
    g_ntt_params.coefficient_bytes   = 2;
    g_ntt_params.poly_bytes          = 128;
    g_ntt_params.max_signers         = 8;
    g_ntt_params.max_rejection_loops = 10;
    g_ntt_params.kappa               = 8;

    // Precompute NTT table once for all benchmarks
    g_ntt_table = NttTable::create(g_ntt_params);
}

/* =========================================================================
   计时工具
   ========================================================================= */
struct BenchResult {
    double avg_us, min_us, max_us, stddev_us, median_us;
    double p99_us, p999_us;
    double avg_cycles;
    double throughput_ops;
    int warmup, timed_iters;
};

static BenchResult bench(const std::function<void()>& func) {
    auto stats = run_benchmark(func, 1000, 50);
    return { stats.avg_us, stats.min_us, stats.max_us,
             stats.stddev_us, stats.median_us,
             stats.p99_us, stats.p999_us, stats.avg_cycles,
             stats.avg_us > 0 ? 1e6 / stats.avg_us : 0,
             50, 1000 };
}

static BenchResult bench_custom(const std::function<void()>& func,
                                 int iters, int warmup) {
    auto stats = run_benchmark(func, iters, warmup);
    return { stats.avg_us, stats.min_us, stats.max_us,
             stats.stddev_us, stats.median_us,
             stats.p99_us, stats.p999_us, stats.avg_cycles,
             stats.avg_us > 0 ? 1e6 / stats.avg_us : 0,
             warmup, iters };
}

/* =========================================================================
   输出辅助
   ========================================================================= */
static void print_and_save(const std::string& title, const BenchResult& r,
                            const std::string& params_desc) {
    std::ostringstream oss;
    oss << "\n=== " << title << " ===\n"
        << "  " << params_desc << "\n"
        << "  Warmup  : " << std::setw(4) << r.warmup << "\n"
        << "  Timed   : " << std::setw(4) << r.timed_iters << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << r.avg_us  << " us\n"
        << "  Median    : " << std::setw(8) << r.median_us << " us\n"
        << "  Min       : " << std::setw(8) << r.min_us  << " us\n"
        << "  Max       : " << std::setw(8) << r.max_us  << " us\n"
        << "  P99       : " << std::setw(8) << r.p99_us  << " us\n"
        << "  P999      : " << std::setw(8) << r.p999_us << " us\n"
        << "  StdDev    : " << std::setw(8) << r.stddev_us << " us\n"
        << "  Throughput: " << std::setw(8) << r.throughput_ops << " ops/s";
    if (r.avg_cycles > 0) {
        oss << "\n  AvgCycles : " << std::setw(8) << std::setprecision(0) << r.avg_cycles;
    }
    oss << "\n";

    std::cout << "\n--- " << title << " ---\n" << oss.str();

    constexpr const char* OUT_PATH = "benchmarking_output/benchmarking_ibags.txt";
    std::ofstream fout(OUT_PATH, std::ios::app);
    if (fout.is_open()) {
        fout << oss.str();
        fout.close();
    }
}

/* =========================================================================
   辅助: 随机多项式生成 (变化种子，避免固定输入)
   ========================================================================= */
static Poly random_poly(const Params& pp, uint64_t seed = 42) {
    std::vector<int64_t> coeffs(pp.n);
    std::mt19937_64 rng(seed);
    std::uniform_int_distribution<int64_t> dist(0, pp.q - 1);
    for (int i = 0; i < pp.n; ++i) coeffs[i] = dist(rng);
    return Poly::from_canonical(pp, std::move(coeffs));
}

/// 为每个迭代生成不同种子的多项式对 (用于 NTT vs naive 对比)
static std::pair<Poly, Poly> random_poly_pair(const Params& pp, int iter) {
    uint64_t base = static_cast<uint64_t>(iter) * 0x9E3779B97F4A7C15ULL;
    return {random_poly(pp, base), random_poly(pp, base ^ 0xDEADBEEF)};
}

static SeedCSPRNG make_deterministic_csprng() {
    std::vector<uint8_t> seed(32, 0);
    for (int i = 0; i < 32; ++i) seed[i] = (uint8_t)(i * 7 + 1);
    return SeedCSPRNG(seed.data(), seed.size(), "IBAGS-v1/CSPRNG/SEED");
}

/* =========================================================================
   §1  mod_reduce — Barrett Reduction (预计算复用, 与 NTT 路径一致)
   ========================================================================= */
static void bench_barrett_reduce(const Params& pp) {
    const int N = 10000;
    std::vector<int64_t> inputs(N);
    std::mt19937_64 rng(99);
    // 输入范围覆盖 NTT 场景: 系数乘 zeta 可到 ~q²，加减可达 ±2q
    std::uniform_int_distribution<int64_t> dist(-static_cast<int64_t>(pp.q) * 2,
                                                 static_cast<int64_t>(pp.q) * 2);
    for (int i = 0; i < N; ++i) inputs[i] = dist(rng);

    // 预计算 Barrett 常数 (模拟 NTT 中的复用路径)
    BarrettConst bc = make_barrett(pp.q);
    volatile int64_t sink = 0;

    auto r = bench([&]() {
        for (int i = 0; i < N; ++i) sink = barrett_reduce(inputs[i], bc);
    });
    (void)sink;

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " (per " << N << " calls, cached BarrettConst)";
    print_and_save("Barrett Reduction (cached)", r, desc.str());
}

/* =========================================================================
   §1b  mod_reduce — Montgomery Multiplication
   ========================================================================= */
static void bench_montgomery_mul(const Params& pp) {
    MontgomeryConst mc = make_montgomery(pp.q);
    const int N = 10000;
    std::mt19937_64 rng(77);
    std::uniform_int_distribution<int64_t> dist(0, pp.q - 1);

    std::vector<int64_t> a_mont(N), b_mont(N);
    for (int i = 0; i < N; ++i) {
        a_mont[i] = to_montgomery(dist(rng), mc);
        b_mont[i] = to_montgomery(dist(rng), mc);
    }
    volatile int64_t sink = 0;

    auto r = bench([&]() {
        for (int i = 0; i < N; ++i) sink = montgomery_mul(a_mont[i], b_mont[i], mc);
    });
    (void)sink;

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " (per " << N << " calls)";
    print_and_save("Montgomery Mul", r, desc.str());
}

/* =========================================================================
   §1c  mod_reduce — Barrett Reduction (constant-time)
   ========================================================================= */
static void bench_barrett_reduce_ct(const Params& pp) {
    const int N = 10000;
    std::vector<int64_t> inputs(N);
    std::mt19937_64 rng(99);
    std::uniform_int_distribution<int64_t> dist(-static_cast<int64_t>(pp.q) * 2,
                                                 static_cast<int64_t>(pp.q) * 2);
    for (int i = 0; i < N; ++i) inputs[i] = dist(rng);

    BarrettConst bc = make_barrett(pp.q);
    volatile int64_t sink = 0;

    auto r = bench([&]() {
        for (int i = 0; i < N; ++i) sink = barrett_reduce_ct(inputs[i], bc);
    });
    (void)sink;

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " (per " << N << " calls, constant-time)";
    print_and_save("Barrett Reduction (CT)", r, desc.str());
}

/* =========================================================================
   §1d  mod_reduce — 多项式全系数约减 reduce_mod_q
   ========================================================================= */
static void bench_reduce_mod_q(const Params& pp) {
    Poly a = random_poly(pp, 777);
    // 把系数放大到非 canonical，模拟环约减后的输入
    int64_t q = static_cast<int64_t>(pp.q);
    for (int i = 0; i < pp.n; ++i) a[i] = a[i] * 100 + i - q;
    auto r = bench([&]() { a.reduce_mod_q(pp); });
    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly reduce_mod_q", r, desc.str());
}

/* =========================================================================
   §1e  NttTable — 预计算开销
   ========================================================================= */
static void bench_ntt_table_create(const Params& pp) {
    // 只跑少量迭代 — 这是一次性开销，不需要大样本
    const int iters = 100;
    const int warmup = 5;

    auto r = bench_custom([&]() {
        volatile auto tbl = NttTable::create(pp);
        (void)tbl;
    }, iters, warmup);

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NttTable::create (precompute)", r, desc.str());
}

/* =========================================================================
   §2  poly — 多项式加法
   ========================================================================= */
static void bench_poly_add() {
    const auto& pp = g_demo_params;
    Poly a = random_poly(pp);
    Poly b = random_poly(pp);

    auto r = bench([&]() { (void)poly_add(a, b, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Add", r, desc.str());
}

/* =========================================================================
   §3  poly — 多项式减法
   ========================================================================= */
static void bench_poly_sub() {
    const auto& pp = g_demo_params;
    Poly a = random_poly(pp);
    Poly b = random_poly(pp);

    auto r = bench([&]() { (void)poly_sub(a, b, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Sub", r, desc.str());
}

/* =========================================================================
   §3b  poly — 多项式取负 O(n)
   ========================================================================= */
static void bench_poly_neg(const Params& pp) {
    Poly a = random_poly(pp);
    auto r = bench([&]() { (void)poly_neg(a, pp); });
    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Neg", r, desc.str());
}

/* =========================================================================
   §4  poly — 朴素乘法 O(n^2)
   ========================================================================= */
static void bench_poly_mul_naive(const Params& pp) {
    Poly a = random_poly(pp);
    Poly b = random_poly(pp);

    auto r = bench([&]() { (void)poly_mul_naive(a, b, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Mul (naive O(n^2))", r, desc.str());
}

/* =========================================================================
   §4b  poly — 标量乘法 O(n)
   ========================================================================= */
static void bench_poly_mul_scalar(const Params& pp) {
    Poly a = random_poly(pp);
    // 使用 NTT zeta 值作为典型标量: 协议中常乘以挑战系数 ∈ [0, q)
    std::mt19937_64 rng(123);
    int64_t scalar = std::uniform_int_distribution<int64_t>(1, pp.q - 1)(rng);

    auto r = bench([&]() { (void)poly_mul_scalar(a, scalar, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Mul Scalar", r, desc.str());
}

/* =========================================================================
   §4c  poly — 环约减 Y^n = -1
   ========================================================================= */
static void bench_poly_ring_reduce_raw(const Params& pp) {
    // 构造一个 (2n-1) 长度的模拟卷积结果
    int n = pp.n;
    std::vector<int64_t> conv(2 * n - 1);
    std::mt19937_64 rng(555);
    std::uniform_int_distribution<int64_t> dist(0, pp.q - 1);
    for (auto& v : conv) v = dist(rng);
    auto r = bench([&]() { (void)poly_ring_reduce_raw(pp, conv); });
    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Ring Reduce", r, desc.str());
}

/* =========================================================================
   §5  poly — 无穷范数
   ========================================================================= */
static void bench_poly_norm_inf() {
    const auto& pp = g_demo_params;
    Poly a = random_poly(pp);

    auto r = bench([&]() { (void)poly_norm_inf(a, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Norm (inf)", r, desc.str());
}

/* =========================================================================
   §5b  poly — 范数边界检查 (constant-time)
   ========================================================================= */
static void bench_poly_norm_bound_check(const Params& pp) {
    Poly a = random_poly(pp);
    uint64_t bound = static_cast<uint64_t>(pp.eta1);
    auto r = bench([&]() { (void)poly_norm_bound_check(a, bound, pp); });
    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " bound=" << bound;
    print_and_save("Poly Norm Bound Check", r, desc.str());
}

/* =========================================================================
   §5c  poly — constant-time 相等比较
   ========================================================================= */
static void bench_poly_equal_ct(const Params& pp) {
    Poly a = random_poly(pp);
    Poly b = random_poly(pp, 99);
    volatile bool sink = false;
    auto r = bench([&]() { sink = poly_equal_ct(a, b); });
    (void)sink;
    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Equal (CT)", r, desc.str());
}

/* =========================================================================
   §6  poly — 多项式求逆 (扩展欧几里得)
   ========================================================================= */
static void bench_poly_inv(const Params& pp) {
    auto csprng = make_deterministic_csprng();
    SecretPoly f_secret = sample_ntru_secret_f(pp, csprng);

    std::vector<int64_t> coeffs(pp.n);
    for (int i = 0; i < pp.n; ++i) coeffs[i] = f_secret.coeff(i);
    Poly f_poly = Poly::from_coeffs(pp, std::move(coeffs));
    f_poly.reduce_mod_q(pp);

    if (!poly_is_invertible_mod_q(f_poly, pp)) {
        std::cout << "  [SKIP] Poly Inv — sample not invertible\n";
        return;
    }

    auto r = bench([&]() {
        Poly inv(pp.n);
        (void)poly_inv(f_poly, pp, &inv);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Inversion", r, desc.str());
}

/* =========================================================================
   §6b  poly — 可逆性快速检查
   ========================================================================= */
static void bench_poly_is_invertible(const Params& pp) {
    auto csprng = make_deterministic_csprng();
    SecretPoly f_secret = sample_ntru_secret_f(pp, csprng);
    std::vector<int64_t> coeffs(pp.n);
    for (int i = 0; i < pp.n; ++i) coeffs[i] = f_secret.coeff(i);
    Poly f_poly = Poly::from_coeffs(pp, coeffs);
    f_poly.reduce_mod_q(pp);
    volatile bool sink = false;
    auto r = bench([&]() { sink = poly_is_invertible_mod_q(f_poly, pp); });
    (void)sink;
    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Is Invertible", r, desc.str());
}

/* =========================================================================
   §7  poly_ntt — 正向 NTT (Cooley-Tukey)
   ========================================================================= */
static void bench_ntt_forward(const Params& pp, const NttTable& tbl) {
    Poly a = random_poly(pp);
    NttPoly ntt(pp.n);

    auto r = bench([&]() { (void)poly_ntt(a, tbl, &ntt); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Forward", r, desc.str());
}

/* =========================================================================
   §8  poly_ntt — 逆向 NTT (Gentleman-Sande)
   ========================================================================= */
static void bench_ntt_inverse(const Params& pp, const NttTable& tbl) {
    Poly a = random_poly(pp);
    NttPoly ntt(pp.n);
    (void)poly_ntt(a, tbl, &ntt);
    Poly out(pp.n);

    auto r = bench([&]() { (void)poly_invntt(ntt, tbl, &out); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Inverse", r, desc.str());
}

/* =========================================================================
   §9  poly_ntt — NTT 往返 (Forward + Inverse)
   ========================================================================= */
static void bench_ntt_roundtrip(const Params& pp, const NttTable& tbl) {
    Poly a = random_poly(pp);
    NttPoly ntt(pp.n);

    auto r = bench([&]() {
        (void)poly_ntt(a, tbl, &ntt);
        Poly out(pp.n);
        (void)poly_invntt(ntt, tbl, &out);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Roundtrip", r, desc.str());
}

/* =========================================================================
   §10  poly_ntt — NTT 乘法 (NTT→Pointwise→INTT)
   ========================================================================= */
static void bench_ntt_mul(const Params& pp, const NttTable& tbl) {
    Poly a = random_poly(pp);
    Poly b = random_poly(pp);
    NttPoly ntt_a(pp.n), ntt_b(pp.n);
    (void)poly_ntt(a, tbl, &ntt_a);
    (void)poly_ntt(b, tbl, &ntt_b);

    NttPoly ntt_c(pp.n);
    Poly c(pp.n);

    auto r = bench([&]() {
        (void)poly_pointwise_mul_ntt(ntt_a, ntt_b, tbl, &ntt_c);
        (void)poly_invntt(ntt_c, tbl, &c);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Mul (Pointwise+INTT)", r, desc.str());
}

/* =========================================================================
   §10b  poly_ntt — 完整 NTT 乘法管线 (2×NTT + PW + INTT)
   ========================================================================= */
static void bench_poly_mul_ntt_full(const Params& pp, const NttTable& tbl) {
    Poly a = random_poly(pp);
    Poly b = random_poly(pp);
    Poly c(pp.n);

    auto r = bench([&]() { (void)poly_mul_ntt(a, b, tbl, &c); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Mul Full (2xNTT+PW+INTT)", r, desc.str());
}

/* =========================================================================
   §10c  poly_ntt — NTT Mul vs Naive Mul 对比 (同输入)
   ========================================================================= */
static void bench_ntt_vs_naive(const Params& pp, const NttTable& tbl, int trials = 5) {
    // 对同一对输入分别测 NTT 乘法和 Naive 乘法 — 完整管线公平对比
    for (int t = 0; t < trials; ++t) {
        auto [a, b] = random_poly_pair(pp, t);

        // Naive (O(n²))
        auto r_naive = bench_custom([&]() { (void)poly_mul_naive(a, b, pp); }, 100, 10);

        // NTT 完整管线 (2×NTT + Pointwise + INTT)
        Poly c(pp.n);
        auto r_ntt = bench_custom([&]() { (void)poly_mul_ntt(a, b, tbl, &c); }, 100, 10);

        double speedup = (r_naive.avg_us > 0) ? r_naive.avg_us / r_ntt.avg_us : 0.0;

        std::ostringstream desc;
        desc << "n=" << pp.n << " q=" << pp.q
             << " trial=" << (t + 1)
             << " | naive=" << std::fixed << std::setprecision(1) << r_naive.avg_us << " us"
             << " ntt(full)=" << r_ntt.avg_us << " us"
             << " speedup=" << std::setprecision(2) << speedup << "x";
        print_and_save("NTT vs Naive Mul", r_ntt, desc.str());
    }
}

/* =========================================================================
   §11  gauss — 单系数离散高斯采样
   ========================================================================= */
static void bench_gauss_sample_coeff() {
    const auto& pp = g_demo_params;
    Xof xof(std::string(domain::GAUSS_FUNCTION));
    std::vector<uint8_t> seed(32, 0x42);
    xof.absorb(ByteSpan(seed.data(), seed.size()));
    xof.finalize();

    auto r = bench([&]() { (void)gauss_sample_coeff(xof, pp.sigma, pp.q); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " sigma=" << pp.sigma;
    print_and_save("Gauss Sample Coeff", r, desc.str());
}

/* =========================================================================
   §12  gauss — 多项式高斯采样 (n 个系数)
   ========================================================================= */
static void bench_gauss_sample_poly() {
    const auto& pp = g_demo_params;
    Xof xof(std::string(domain::GAUSS_FUNCTION));
    std::vector<uint8_t> seed(32, 0x42);
    xof.absorb(ByteSpan(seed.data(), seed.size()));
    xof.finalize();

    auto r = bench([&]() { (void)gauss_sample_poly(xof, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " sigma=" << pp.sigma;
    print_and_save("Gauss Sample Poly", r, desc.str());
}

/* =========================================================================
   §13  xof — SHAKE256 squeeze 吞吐量
   ========================================================================= */
static void bench_xof_squeeze() {
    Xof xof(std::string(domain::H1_TO_RING));
    std::vector<uint8_t> data(64, 0xAB);
    xof.absorb(ByteSpan(data.data(), data.size()));
    xof.finalize();

    const int N = 1000;
    auto r = bench([&]() {
        for (int i = 0; i < N; ++i) (void)xof.squeeze_u64();
    });

    std::ostringstream desc;
    desc << "SHAKE256 (per " << N << " calls)";
    print_and_save("XOF squeeze_u64", r, desc.str());
}

/* =========================================================================
   §14  csprng — SystemCSPRNG (OS entropy)
   ========================================================================= */
static void bench_csprng_system() {
    SystemCSPRNG rng;

    const int N = 5000;
    auto r = bench([&]() {
        for (int i = 0; i < N; ++i) (void)rng.rand_u64();
    });

    std::ostringstream desc;
    desc << "System entropy (per " << N << " calls)";
    print_and_save("SystemCSPRNG rand_u64", r, desc.str());
}

/* =========================================================================
   §15  csprng — SeedCSPRNG (deterministic)
   ========================================================================= */
static void bench_csprng_seed() {
    auto rng = make_deterministic_csprng();

    const int N = 10000;
    auto r = bench([&]() {
        for (int i = 0; i < N; ++i) (void)rng.rand_u64();
    });

    std::ostringstream desc;
    desc << "Seed (per " << N << " calls)";
    print_and_save("SeedCSPRNG rand_u64", r, desc.str());
}

/* =========================================================================
   §16  encode — 多项式编码 (n × 16 字节)
   ========================================================================= */
static void bench_poly_encode() {
    const auto& pp = g_demo_params;
    Poly a = random_poly(pp);

    std::vector<uint8_t> buf(pp.poly_bytes);
    auto r = bench([&]() { (void)poly_encode(a, pp, buf.data(), buf.size()); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " bytes=" << pp.poly_bytes;
    print_and_save("Poly Encode", r, desc.str());
}

/* =========================================================================
   §17  encode — 多项式解码
   ========================================================================= */
static void bench_poly_decode() {
    const auto& pp = g_demo_params;
    Poly a = random_poly(pp);
    std::vector<uint8_t> buf(pp.poly_bytes);
    (void)poly_encode(a, pp, buf.data(), buf.size());

    auto r = bench([&]() {
        Poly dec(pp.n);
        (void)poly_decode(buf.data(), buf.size(), pp, &dec);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " bytes=" << pp.poly_bytes;
    print_and_save("Poly Decode", r, desc.str());
}

/* =========================================================================
   §18  hash — sample_ring_element (H1 core)
   ========================================================================= */
static void bench_sample_ring_element() {
    const auto& pp = g_demo_params;

    Xof xof(std::string(domain::H1_TO_RING));
    std::vector<uint8_t> data(64, 0xAB);
    xof.absorb(ByteSpan(data.data(), data.size()));
    xof.finalize();

    Poly out(pp.n);
    auto r = bench([&]() { sample_ring_element(xof, pp, &out); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("H1 sample_ring_element", r, desc.str());
}

/* =========================================================================
   §19  hash — sample_challenge_polynomial (H2 core)
   ========================================================================= */
static void bench_sample_challenge_polynomial() {
    const auto& pp = g_demo_params;

    Xof xof(std::string(domain::H2_CHALLENGE));
    std::vector<uint8_t> data(64, 0xCD);
    xof.absorb(ByteSpan(data.data(), data.size()));
    xof.finalize();

    Poly out(pp.n);
    auto r = bench([&]() { sample_challenge_polynomial(xof, pp, pp.kappa, &out); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " kappa=" << pp.kappa;
    print_and_save("H2 sample_challenge_poly", r, desc.str());
}

/* =========================================================================
   §20  hash — coeff_reject_small (H3 core)
   ========================================================================= */
static void bench_coeff_reject_small() {
    const auto& pp = g_demo_params;

    Xof xof(std::string(domain::H3_AGG_COEFF));
    std::vector<uint8_t> data(64, 0xEF);
    xof.absorb(ByteSpan(data.data(), data.size()));
    xof.finalize();

    auto r = bench([&]() { (void)coeff_reject_small(xof, pp.eta2, 24); });

    std::ostringstream desc;
    desc << "B=" << pp.eta2;
    print_and_save("coeff_reject_small", r, desc.str());
}


/* =========================================================================
   §21  reject — hash_to_uniform (n 系数拒绝采样)
   ========================================================================= */
static void bench_hash_to_uniform() {
    const auto& pp = g_demo_params;

    Xof xof(std::string(domain::H1_TO_RING));
    std::vector<uint8_t> data(64, 0xAB);
    xof.absorb(ByteSpan(data.data(), data.size()));
    xof.finalize();

    auto r = bench([&]() { (void)hash_to_uniform(xof, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Rejection Sampling (hash_to_uniform)", r, desc.str());
}

/* =========================================================================
   §22  ntru_secret — 采样 f+g 密钥对
   ========================================================================= */
static void bench_ntru_secret_pair() {
    const auto& pp = g_demo_params;
    auto csprng = make_deterministic_csprng();

    SecretPoly f(pp.n), g(pp.n);

    auto r = bench([&]() {
        sample_ntru_secret_pair(pp, csprng, &f, &g);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " sigma=" << pp.sigma;
    print_and_save("NTRU Secret Pair (f,g)", r, desc.str());
}

/* =========================================================================
   §23  ntru_trapgen — 完整 TrapGen 流水线
   ========================================================================= */
static void bench_ntru_trapgen() {
    const auto& pp = g_ntt_params;
    auto csprng = make_deterministic_csprng();
    auto config = TrapGenConfig::default_config();

    auto r = bench_custom([&]() {
        PublicTrapdoorParams pub(pp);
        MasterTrapdoorSecret sec(pp.n);
        (void)ntru_trapgen(pp, csprng, config, &pub, &sec);
    }, 5, 1);

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " sigma=" << pp.sigma;
    print_and_save("NTRU TrapGen (full pipeline)", r, desc.str());
}

/* =========================================================================
   §23b  ntru_trapgen — 陷门验证
   ========================================================================= */
static void bench_check_trapdoor_basis(const Params& pp) {
    auto csprng = make_deterministic_csprng();
    auto config = TrapGenConfig::default_config();
    config.max_invertibility_attempts = 10;
    PublicTrapdoorParams pub(pp);
    MasterTrapdoorSecret sec(pp.n);
    auto st = ntru_trapgen(pp, csprng, config, &pub, &sec);
    if (!st.ok()) { std::cout << "  [SKIP] Trapdoor check — TrapGen failed\n"; return; }
    volatile bool sink = false;
    auto r = bench_custom([&]() { sink = check_trapdoor_basis(pp, sec.basis, pub.h); }, 100, 10);
    (void)sink;
    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Check Trapdoor Basis", r, desc.str());
}

/* =========================================================================
   §24  ntru_trapgen — sample_f_until_invertible
   ========================================================================= */
static void bench_sample_f_until_invertible() {
    const auto& pp = g_ntt_params;
    auto csprng = make_deterministic_csprng();

    auto r = bench([&]() {
        Poly f_out(pp.n);
        (void)sample_f_until_invertible(pp, csprng, 100, &f_out);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTRU sample_f_until_invertible", r, desc.str());
}

/* =========================================================================
   多级参数表 — NTT 相关基准在所有安全级别上运行
   ========================================================================= */
struct ParamEntry {
    const char* name;
    Params      params;
};

static std::vector<ParamEntry> make_ntt_param_entries() {
    return {
        {"Demo (n=64)",       Params::params_demo_64()},
        {"Level 1 (n=512)",   Params::params_level2_512()},
        {"Level 3 (n=1024)",  Params::params_level3_1024()},
        {"Level 5 (n=1024)",  Params::params_level5_1024()},
    };
}

static std::string make_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

/* =========================================================================
   main — 统一基准测试入口 (多级 NTT + 单级非 NTT)
   ========================================================================= */
int main() {
    init_params();

    constexpr const char* OUT_PATH = "benchmarking_output/benchmarking_ibags.txt";

    // ── 文件头 (覆盖 + 时间戳) ──
    {
        std::ofstream fout(OUT_PATH, std::ios::trunc);
        if (fout.is_open()) {
            fout << "==========================================================\n"
                 << "  IBAGS (RLWE/NTRU) — Benchmark Run\n"
                 << "  Timestamp: " << make_timestamp() << "\n"
                 << "==========================================================\n\n";
            fout.close();
        }
    }

    std::cout << "============================================================\n"
              << "  IBAGS (RLWE/NTRU) — Cryptographic Operation Benchmarks\n"
              << "  Timestamp: " << make_timestamp() << "\n"
              << "============================================================\n\n";

    // ════════════════════════════════════════════════════════════════
    // §A  NTT 专项基准 — 遍历所有安全级别
    // ════════════════════════════════════════════════════════════════
    auto entries = make_ntt_param_entries();
    for (const auto& entry : entries) {
        const Params& pp = entry.params;
        if (!is_ntt_friendly(pp)) {
            std::cout << "\n--- SKIP " << entry.name
                      << " (q not NTT-friendly) ---\n";
            continue;
        }
        NttTable tbl = NttTable::create(pp);

        // 级别分隔标题
        {
            std::ostringstream header;
            header << "\n========== " << entry.name
                   << " n=" << pp.n << " q=" << pp.q << " ==========\n";
            std::cout << header.str();
            std::ofstream fout(OUT_PATH, std::ios::app);
            if (fout.is_open()) { fout << header.str(); fout.close(); }
        }

        // 模约减 (每个级别 q 不同)
        bench_barrett_reduce(pp);
        bench_barrett_reduce_ct(pp);
        bench_montgomery_mul(pp);
        bench_reduce_mod_q(pp);

        // NTT 预计算开销
        bench_ntt_table_create(pp);

        // NTT 核心操作
        bench_ntt_forward(pp, tbl);
        bench_ntt_inverse(pp, tbl);
        bench_ntt_roundtrip(pp, tbl);
        bench_ntt_mul(pp, tbl);
        bench_poly_mul_ntt_full(pp, tbl);

        // 基础多项式操作 (所有级别 O(n))
        bench_poly_mul_scalar(pp);
        bench_poly_neg(pp);
        bench_poly_norm_bound_check(pp);
        bench_poly_equal_ct(pp);
        bench_poly_ring_reduce_raw(pp);

        // O(n²) 操作 (仅 Demo 和 L1)
        if (pp.n <= 512) {
            bench_poly_mul_naive(pp);
            bench_poly_inv(pp);
            bench_poly_is_invertible(pp);
            bench_check_trapdoor_basis(pp);
        }

        // NTT vs Naive 对比
        int vs_trials = (pp.n <= 512) ? 5 : 1;
        bench_ntt_vs_naive(pp, tbl, vs_trials);
    }

    // ════════════════════════════════════════════════════════════════
    // §B  非 NTT 基准 — 仅在 Demo 参数运行 (与多项式/签名逻辑相关)
    // ════════════════════════════════════════════════════════════════
    std::cout << "\n========== Non-NTT Operations (Demo) ==========\n";
    {
        std::ofstream fout(OUT_PATH, std::ios::app);
        if (fout.is_open()) {
            fout << "\n========== Non-NTT Operations (Demo) ==========\n";
            fout.close();
        }
    }

    bench_poly_add();
    bench_poly_sub();
    bench_poly_mul_scalar(g_demo_params);
    bench_poly_mul_naive(g_demo_params);
    bench_poly_norm_inf();
    bench_poly_inv(g_demo_params);

    bench_gauss_sample_coeff();
    bench_gauss_sample_poly();
    bench_xof_squeeze();
    bench_csprng_system();
    bench_csprng_seed();
    bench_poly_encode();
    bench_poly_decode();
    bench_sample_ring_element();
    bench_sample_challenge_polynomial();
    bench_coeff_reject_small();
    bench_hash_to_uniform();
    bench_ntru_secret_pair();
    bench_sample_f_until_invertible();
    bench_ntru_trapgen();

    std::cout << "\n============================================================\n"
              << "  IBAGS benchmark complete.\n"
              << "  Results written to " << OUT_PATH << "\n"
              << "============================================================\n";

    return 0;
}
