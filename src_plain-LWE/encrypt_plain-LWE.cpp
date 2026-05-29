/**
 * encrypt_plain-LWE.cpp — Encryption primitives test and benchmark
 *
 * Tests and benchmarks the three encryption layers of plain-LWE FHE:
 *   ① LWE Encrypt (single bit)       — standard LWE encryption
 *   ② GSW Encrypt (full ciphertext)  — GSW matrix encryption: C = A·S + μ·G
 *   ③ UniEnc (mask scheme)           — C = A·R + μ·G  with LWE-encrypted R bits
 */

#include "unienc_plain-LWE.h"
#include "matops_plain-LWE.h"
#include "eval_plain-LWE.h"
#include "unified_params_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace matops;
using namespace unienc;

/* ───── 文件输出 ───── */
#include "bench_utils_plain-LWE.h"
static std::ostringstream enc_oss;

/* ═══════════════════════════════════════════════════
   §1  LWE Encrypt 单比特加密正确性测试
   ═══════════════════════════════════════════════════ */
static void test_lwe_encrypt_correctness() {
    std::cout << "\n--- Test 1: LWE Encrypt (single bit) correctness ---\n";
    auto mp = unified::default_mp12_params();
    const long q = mp.q;
    const int  N = mp.n + 1;
    const int  m = 20;

    std::mt19937_64 rng(42);
    std::uniform_int_distribution<long> uni(0, q - 1);

    Mat A = make_mat(N, m);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < m; ++j) A[i][j] = uni(rng);

    int trials = 10, pass = 0;
    for (int t = 0; t < trials; ++t) {
        long mu = t % 2;
        Vec V = lwe_encrypt_bit(A, mu, q, 3.2, rng);
        if ((int)V.size() == N) pass++;
        if (t < 2) {
            std::cout << "  trial " << t << ": mu=" << mu
                      << "  |V|=" << V.size() << "\n";
        }
    }
    std::cout << "  Dimension check: " << pass << "/" << trials
              << (pass == trials ? " PASS" : " FAIL") << "\n";
}

/* ═══════════════════════════════════════════════════
   §2  GSW Encrypt 正确性验证
   ═══════════════════════════════════════════════════ */
static Mat gsw_encrypt(const Mat& A, const Mat& G, int mu, long q, std::mt19937_64& rng) {
    int R = (int)A.size(), M = (int)A[0].size();
    std::uniform_int_distribution<int> bit(0, 1);
    Mat S = make_mat(M, M);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j) S[i][j] = bit(rng);
    Mat AS = mat_mul(A, S, q);
    Mat muG = make_mat(R, M, 0);
    if (mu) for (int i = 0; i < R; ++i)
                for (int j = 0; j < M; ++j)
                    muG[i][j] = mod_pos((long)mu * G[i][j], q);
    return mat_add(AS, muG, q);
}

static int gsw_trial_decrypt(const Vec& t, const Mat& C, const Vec& w_hat, long q, int b) {
    int R = (int)C.size(), M = (int)C[0].size();
    Vec tC(M, 0);
    for (int i = 0; i < M; ++i) {
        long acc = 0;
        for (int r = 0; r < R; ++r) acc = mod_pos(acc + t[r] * C[r][i], q);
        tC[i] = acc;
    }
    Mat w_col(R, Vec(1));
    for (int i = 0; i < R; ++i) w_col[i][0] = w_hat[i];
    Mat u_col = cryptolib::gadget_inverse(w_col, q, b);
    long inner = 0;
    for (int i = 0; i < M; ++i) inner = mod_pos(inner + tC[i] * u_col[i][0], q);
    long c = inner; if (c > q/2) c -= q;
    return (c < 0 ? -c : c) > (q/4) ? 1 : 0;
}

static void test_gsw_encrypt_correctness() {
    std::cout << "\n--- Test 2: GSW Encrypt correctness ---\n";
    auto mp = unified::default_mp12_params();
    const long q = mp.q; const int b = mp.b;
    const int k = cryptolib::eval_compute_k(q, b);
    const int N = 8, M = N * k;

    std::mt19937_64 rng(99);
    std::uniform_int_distribution<long> uni(0, q-1);

    Mat G = cryptolib::build_gadget(N, q, b);
    Mat A = make_mat(N, M);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < M; ++j) A[i][j] = uni(rng);

    Vec t(N, 0); t[N-1] = 1;
    Vec w_hat(N, 0); w_hat[N-1] = (q+1)/2;

    int pass = 0, trials = 10;
    for (int tr = 0; tr < trials; ++tr) {
        int mu = tr % 2;
        Mat C = gsw_encrypt(A, G, mu, q, rng);
        if (gsw_trial_decrypt(t, C, w_hat, q, b) == mu) pass++;
    }
    std::cout << "  Roundtrip: " << pass << "/" << trials
              << "  (note: uses random A without LWE structure; "
              << "full LWE-based test in decrypt_plain-LWE.cpp)\n";
}

/* ═══════════════════════════════════════════════════
   §3  UniEnc 正确性验证
   ═══════════════════════════════════════════════════ */
static void test_uni_enc_correctness() {
    std::cout << "\n--- Test 3: UniEnc correctness ---\n";
    auto mp = unified::default_mp12_params();
    const long q = mp.q; const int n = mp.n; const int m = 4;
    auto params = unienc::Params::make(n, 0, m, q, 3.2);

    std::mt19937_64 rng(123);
    std::uniform_int_distribution<long> uni(0, q-1);
    Mat A = make_mat(params.N, params.m);
    for (int i = 0; i < params.N; ++i)
        for (int j = 0; j < params.m; ++j) A[i][j] = uni(rng);
    Mat Gg = unienc::make_gadget(params.N, params.m, mp.b, q);

    for (int mu = 0; mu <= 1; ++mu) {
        auto out = uni_enc(A, Gg, mu, params, (uint64_t)(mu+1)*777);
        bool ok = (int)out.C.size() == params.N && (int)out.C[0].size() == params.m
                  && (int)out.U.size() == params.m * params.m;
        std::cout << "  mu=" << mu << " C:" << out.C.size() << "x" << out.C[0].size()
                  << " |U|=" << out.U.size() << (ok ? " PASS" : " FAIL") << "\n";
    }
}

static void bench_lwe_encrypt() {
    auto mp = unified::default_mp12_params();
    long q = mp.q; int Ndim = mp.n + 1, mdim = 20;
    std::mt19937_64 rng(42);
    Mat A(Ndim, Vec(mdim));
    for (int i = 0; i < Ndim; ++i)
        for (int j = 0; j < mdim; ++j)
            A[i][j] = std::uniform_int_distribution<long>(0,q-1)(rng);

    auto stats = run_benchmark([&]{ (void)lwe_encrypt_bit(A, 0, q, 3.2, rng); }, 1000, 50);
    std::ostringstream oss;
    oss << "\n=== Benchmark: LWE Encrypt (single bit) ===\n"
        << "  Parameters: q=" << q << ", N=" << Ndim << ", m=" << mdim << ", sigma=3.2\n"
        << "  Total iterations: 1000  Warmup: 50  Active samples: " << stats.active_samples << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << stats.avg_us    << " us\n"
        << "  Median    : " << std::setw(8) << stats.median_us << " us\n"
        << "  Min       : " << std::setw(8) << stats.min_us    << " us\n"
        << "  Max       : " << std::setw(8) << stats.max_us    << " us\n"
        << "  StdDev    : " << std::setw(8) << stats.stddev_us << " us\n"
        << "  Throughput: " << std::setw(8) << (1e6/stats.avg_us) << " ops/s\n";
    std::cout << "\n--- Benchmark: LWE Encrypt ---\n" << oss.str();
    enc_oss << oss.str();
}

static void bench_gsw_encrypt() {
    auto mp = unified::default_mp12_params();
    long q = mp.q; int b = mp.b;
    int k = cryptolib::eval_compute_k(q, b), Ndim = 8, Mdim = Ndim * k;
    std::mt19937_64 rng(99);
    Mat G = cryptolib::build_gadget(Ndim, q, b);
    Mat A = make_mat(Ndim, Mdim);
    for (int i = 0; i < Ndim; ++i)
        for (int j = 0; j < Mdim; ++j)
            A[i][j] = std::uniform_int_distribution<long>(0,q-1)(rng);

    auto stats = run_benchmark([&]{ (void)gsw_encrypt(A, G, 1, q, rng); }, 1000, 50);
    std::ostringstream oss;
    oss << "\n=== Benchmark: GSW Encrypt (full ciphertext) ===\n"
        << "  Parameters: q=" << q << ", N=" << Ndim << ", M=" << Mdim << "\n"
        << "  Total iterations: 1000  Warmup: 50  Active samples: " << stats.active_samples << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << stats.avg_us    << " us\n"
        << "  Median    : " << std::setw(8) << stats.median_us << " us\n"
        << "  Min       : " << std::setw(8) << stats.min_us    << " us\n"
        << "  Max       : " << std::setw(8) << stats.max_us    << " us\n"
        << "  StdDev    : " << std::setw(8) << stats.stddev_us << " us\n"
        << "  Throughput: " << std::setw(8) << (1e6/stats.avg_us) << " ops/s\n";
    std::cout << "\n--- Benchmark: GSW Encrypt ---\n" << oss.str();
    enc_oss << oss.str();
}

static void bench_uni_enc() {
    auto mp = unified::default_mp12_params();
    long q = mp.q; int n = mp.n, mdim = 4;
    auto params = unienc::Params::make(n, 0, mdim, q, 3.2);
    std::mt19937_64 rng(99);
    Mat A = make_mat(params.N, params.m);
    for (int i = 0; i < params.N; ++i)
        for (int j = 0; j < params.m; ++j)
            A[i][j] = std::uniform_int_distribution<long>(0,q-1)(rng);
    Mat Gg = unienc::make_gadget(params.N, params.m, mp.b, q);

    auto stats = run_benchmark([&]{ (void)uni_enc(A, Gg, 1, params, 777); }, 1000, 50);
    std::ostringstream oss;
    oss << "\n=== Benchmark: UniEnc (A*R + mu*G with LWE-encrypted R) ===\n"
        << "  Parameters: q=" << q << ", N=" << params.N << ", m=" << mdim << "\n"
        << "  Total iterations: 1000  Warmup: 50  Active samples: " << stats.active_samples << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << stats.avg_us    << " us\n"
        << "  Median    : " << std::setw(8) << stats.median_us << " us\n"
        << "  Min       : " << std::setw(8) << stats.min_us    << " us\n"
        << "  Max       : " << std::setw(8) << stats.max_us    << " us\n"
        << "  StdDev    : " << std::setw(8) << stats.stddev_us << " us\n"
        << "  Throughput: " << std::setw(8) << (1e6/stats.avg_us) << " ops/s\n";
    std::cout << "\n--- Benchmark: UniEnc ---\n" << oss.str();
    enc_oss << oss.str();
}

/* ═══════════════════════════════════════════════════
   UniEnc sub-operation benchmarks (1-4, 5=lwe_encrypt already exists)
   ═══════════════════════════════════════════════════ */

// 1 — sample_binary_mat
static void bench_uni_sub01_binary_mat() {
    int m = 20;
    std::mt19937_64 rng(42);
    auto stats = run_benchmark([&]{ (void)unienc::sample_binary_mat(m, m, rng); }, 1000, 50);
    std::ostringstream oss;
    oss << "\n=== Benchmark: UniEnc.1 — sample_binary_mat (R in {0,1}^(mxm)) ===\n"
        << "  Parameters: m=" << m << "\n"
        << "  Total iterations: 1000  Warmup: 50  Active samples: " << stats.active_samples << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average: " << std::setw(8) << stats.avg_us    << " us  Min: " << std::setw(8) << stats.min_us
        << "  Max: " << std::setw(8) << stats.max_us << "  Median: " << std::setw(8) << stats.median_us
        << "  StdDev: " << std::setw(8) << stats.stddev_us
        << "  Ops/s: " << std::setw(8) << (1e6/stats.avg_us) << "\n";
    std::cout << "\n--- Sub-bench: UniEnc.1 sample_binary_mat ---\n" << oss.str();
    enc_oss << oss.str();
}

// 2 — mat_mul(A, R)
static void bench_uni_sub02_mat_mul_AR() {
    auto mp = unified::default_mp12_params();
    long q = mp.q; int Ndim = mp.n + 1, mdim = 20;
    std::mt19937_64 rng(99);
    Mat A = make_mat(Ndim, mdim);
    Mat R = unienc::sample_binary_mat(mdim, mdim, rng);
    for (int i = 0; i < Ndim; ++i)
        for (int j = 0; j < mdim; ++j) A[i][j] = std::uniform_int_distribution<long>(0,q-1)(rng);
    auto stats = run_benchmark([&]{ (void)mat_mul(A, R, q); }, 1000, 50);
    std::ostringstream oss;
    oss << "\n=== Benchmark: UniEnc.2 — mat_mul(A*R) ===\n"
        << "  Parameters: q=" << q << ", N=" << Ndim << ", m=" << mdim << "\n"
        << "  Total iterations: 1000  Warmup: 50  Active samples: " << stats.active_samples << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average: " << std::setw(8) << stats.avg_us    << " us  Min: " << std::setw(8) << stats.min_us
        << "  Max: " << std::setw(8) << stats.max_us << "  Median: " << std::setw(8) << stats.median_us
        << "  StdDev: " << std::setw(8) << stats.stddev_us
        << "  Ops/s: " << std::setw(8) << (1e6/stats.avg_us) << "\n";
    std::cout << "\n--- Sub-bench: UniEnc.2 mat_mul ---\n" << oss.str();
    enc_oss << oss.str();
}

// 3 — scalar_mat_mul(mu, G)
static void bench_uni_sub03_scalar_mul_G() {
    auto mp = unified::default_mp12_params();
    long q = mp.q; int Ndim = mp.n + 1, mdim = 20;
    Mat Gg = unienc::make_gadget(Ndim, mdim, mp.b, q);
    auto stats = run_benchmark([&]{ (void)unienc::scalar_mat_mul(1, Gg, q); }, 1000, 50);
    std::ostringstream oss;
    oss << "\n=== Benchmark: UniEnc.3 — scalar_mat_mul(mu*G) ===\n"
        << "  Parameters: q=" << q << ", N=" << Ndim << ", m=" << mdim << "\n"
        << "  Total iterations: 1000  Warmup: 50  Active samples: " << stats.active_samples << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average: " << std::setw(8) << stats.avg_us    << " us  Min: " << std::setw(8) << stats.min_us
        << "  Max: " << std::setw(8) << stats.max_us << "  Median: " << std::setw(8) << stats.median_us
        << "  StdDev: " << std::setw(8) << stats.stddev_us
        << "  Ops/s: " << std::setw(8) << (1e6/stats.avg_us) << "\n";
    std::cout << "\n--- Sub-bench: UniEnc.3 scalar_mul ---\n" << oss.str();
    enc_oss << oss.str();
}

// 4 — mat_add(AR, mu*G)
static void bench_uni_sub04_mat_add_C() {
    auto mp = unified::default_mp12_params();
    long q = mp.q; int Ndim = mp.n + 1, mdim = 20;
    std::mt19937_64 rng(99);
    Mat A = make_mat(Ndim, mdim), B = make_mat(Ndim, mdim);
    for (int i = 0; i < Ndim; ++i)
        for (int j = 0; j < mdim; ++j)
            { A[i][j] = std::uniform_int_distribution<long>(0,q-1)(rng);
              B[i][j] = std::uniform_int_distribution<long>(0,q-1)(rng); }
    auto stats = run_benchmark([&]{ (void)mat_add(A, B, q); }, 1000, 50);
    std::ostringstream oss;
    oss << "\n=== Benchmark: UniEnc.4 — mat_add(C = AR + mu*G) ===\n"
        << "  Parameters: q=" << q << ", N=" << Ndim << ", m=" << mdim << "\n"
        << "  Total iterations: 1000  Warmup: 50  Active samples: " << stats.active_samples << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average: " << std::setw(8) << stats.avg_us    << " us  Min: " << std::setw(8) << stats.min_us
        << "  Max: " << std::setw(8) << stats.max_us << "  Median: " << std::setw(8) << stats.median_us
        << "  StdDev: " << std::setw(8) << stats.stddev_us
        << "  Ops/s: " << std::setw(8) << (1e6/stats.avg_us) << "\n";
    std::cout << "\n--- Sub-bench: UniEnc.4 mat_add ---\n" << oss.str();
    enc_oss << oss.str();
}

/* ═══════════════════════════════════════════════════
   Main entry
   ═══════════════════════════════════════════════════ */
void run_test_encrypt() {
    enc_oss.str(""); enc_oss.clear();
    enc_oss << "=== Benchmark: Encryption Primitives ===\n";
    std::cout << "==========================================================\n"
              << "  Encryption Primitives — LWE / GSW / UniEnc\n"
              << "==========================================================\n";
    test_lwe_encrypt_correctness();
    test_gsw_encrypt_correctness();
    test_uni_enc_correctness();
    bench_lwe_encrypt();           // 5
    bench_gsw_encrypt();
    bench_uni_enc();               // 1-5 combined
    bench_uni_sub01_binary_mat();  // 1
    bench_uni_sub02_mat_mul_AR();  // 2
    bench_uni_sub03_scalar_mul_G();// 3
    bench_uni_sub04_mat_add_C();   // 4
    bench_write(enc_oss.str());
    std::cout << "\nDone.\n";
}
