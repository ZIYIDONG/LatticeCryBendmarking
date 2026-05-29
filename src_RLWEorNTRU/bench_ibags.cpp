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
static Params g_demo_params;
static Params g_ntt_params;

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
}

/* =========================================================================
   计时工具
   ========================================================================= */
struct BenchResult {
    double avg_us;
    double min_us;
    double max_us;
    double stddev_us;
    double median_us;
    double throughput_ops;
    int warmup;
    int timed_iters;
};

static BenchResult bench(const std::function<void()>& func) {
    auto stats = run_benchmark(func, 1000, 50);
    return { stats.avg_us, stats.min_us, stats.max_us,
             stats.stddev_us, stats.median_us,
             stats.avg_us > 0 ? 1e6 / stats.avg_us : 0,
             50, 1000 };
}

static BenchResult bench_custom(const std::function<void()>& func,
                                 int iters, int warmup) {
    auto stats = run_benchmark(func, iters, warmup);
    return { stats.avg_us, stats.min_us, stats.max_us,
             stats.stddev_us, stats.median_us,
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
        << "  StdDev    : " << std::setw(8) << r.stddev_us << " us\n"
        << "  Throughput: " << std::setw(8) << r.throughput_ops << " ops/s\n";

    std::cout << "\n--- " << title << " ---\n" << oss.str();

    constexpr const char* OUT_PATH = "benchmarking_output/benchmarking_ibags.txt";
    std::ofstream fout(OUT_PATH, std::ios::app);
    if (fout.is_open()) {
        fout << oss.str();
        fout.close();
    }
}

/* =========================================================================
   辅助: 随机多项式生成
   ========================================================================= */
static Poly random_poly(const Params& pp) {
    std::vector<int64_t> coeffs(pp.n);
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(0, pp.q - 1);
    for (int i = 0; i < pp.n; ++i) coeffs[i] = dist(rng);
    return Poly::from_canonical(pp, std::move(coeffs));
}

static SeedCSPRNG make_deterministic_csprng() {
    std::vector<uint8_t> seed(32, 0);
    for (int i = 0; i < 32; ++i) seed[i] = (uint8_t)(i * 7 + 1);
    return SeedCSPRNG(seed.data(), seed.size(), "IBAGS-v1/CSPRNG/SEED");
}

/* =========================================================================
   §1  mod_reduce — Barrett Reduction
   ========================================================================= */
static void bench_barrett_reduce() {
    const auto& pp = g_demo_params;
    const int N = 10000;
    std::vector<int64_t> inputs(N);
    std::mt19937_64 rng(99);
    std::uniform_int_distribution<int64_t> dist(0, pp.q * 3);
    for (int i = 0; i < N; ++i) inputs[i] = dist(rng);

    auto r = bench([&]() {
        for (int i = 0; i < N; ++i) (void)barrett_reduce(inputs[i], pp.q);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q << " (per " << N << " calls)";
    print_and_save("Barrett Reduction", r, desc.str());
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
   §4  poly — 朴素乘法 O(n^2)
   ========================================================================= */
static void bench_poly_mul_naive() {
    const auto& pp = g_demo_params;
    Poly a = random_poly(pp);
    Poly b = random_poly(pp);

    auto r = bench([&]() { (void)poly_mul_naive(a, b, pp); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("Poly Mul (naive O(n^2))", r, desc.str());
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
   §6  poly — 多项式求逆 (扩展欧几里得)
   ========================================================================= */
static void bench_poly_inv() {
    const auto& pp = g_ntt_params;

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
   §7  poly_ntt — 正向 NTT (Cooley-Tukey)
   ========================================================================= */
static void bench_ntt_forward() {
    const auto& pp = g_ntt_params;
    if (!is_ntt_friendly(pp)) {
        std::cout << "  [SKIP] NTT Forward — q not NTT-friendly\n";
        return;
    }
    Poly a = random_poly(pp);
    NttPoly ntt(pp.n);

    auto r = bench([&]() { (void)poly_ntt(a, pp, &ntt); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Forward", r, desc.str());
}

/* =========================================================================
   §8  poly_ntt — 逆向 NTT (Gentleman-Sande)
   ========================================================================= */
static void bench_ntt_inverse() {
    const auto& pp = g_ntt_params;
    if (!is_ntt_friendly(pp)) {
        std::cout << "  [SKIP] NTT Inverse — q not NTT-friendly\n";
        return;
    }
    Poly a = random_poly(pp);
    NttPoly ntt(pp.n);
    (void)poly_ntt(a, pp, &ntt);
    Poly out(pp.n);

    auto r = bench([&]() { (void)poly_invntt(ntt, pp, &out); });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Inverse", r, desc.str());
}

/* =========================================================================
   §9  poly_ntt — NTT 往返 (Forward + Inverse)
   ========================================================================= */
static void bench_ntt_roundtrip() {
    const auto& pp = g_ntt_params;
    if (!is_ntt_friendly(pp)) {
        std::cout << "  [SKIP] NTT Roundtrip — q not NTT-friendly\n";
        return;
    }
    Poly a = random_poly(pp);
    NttPoly ntt(pp.n);

    auto r = bench([&]() {
        (void)poly_ntt(a, pp, &ntt);
        Poly out(pp.n);
        (void)poly_invntt(ntt, pp, &out);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Roundtrip", r, desc.str());
}

/* =========================================================================
   §10  poly_ntt — NTT 乘法 (NTT→Pointwise→INTT)
   ========================================================================= */
static void bench_ntt_mul() {
    const auto& pp = g_ntt_params;
    if (!is_ntt_friendly(pp)) {
        std::cout << "  [SKIP] NTT Mul — q not NTT-friendly\n";
        return;
    }
    Poly a = random_poly(pp);
    Poly b = random_poly(pp);
    NttPoly ntt_a(pp.n), ntt_b(pp.n);
    (void)poly_ntt(a, pp, &ntt_a);
    (void)poly_ntt(b, pp, &ntt_b);

    NttPoly ntt_c(pp.n);
    Poly c(pp.n);

    auto r = bench([&]() {
        (void)poly_pointwise_mul_ntt(ntt_a, ntt_b, pp, &ntt_c);
        (void)poly_invntt(ntt_c, pp, &c);
    });

    std::ostringstream desc;
    desc << "n=" << pp.n << " q=" << pp.q;
    print_and_save("NTT Mul (Pointwise+INTT)", r, desc.str());
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
   main — 统一基准测试入口
   ========================================================================= */
int main() {
    init_params();

    std::cout << "============================================================\n"
              << "  IBAGS (RLWE/NTRU) — Cryptographic Operation Benchmarks\n"
              << "============================================================\n"
              << "  Demo  Params: n=" << g_demo_params.n
              << " q=" << g_demo_params.q
              << " sigma=" << g_demo_params.sigma << "\n"
              << "  NTT   Params: n=" << g_ntt_params.n
              << " q=" << g_ntt_params.q << "\n"
              << "============================================================\n\n";

    {
        constexpr const char* OUT_PATH = "benchmarking_output/benchmarking_ibags.txt";
        std::ofstream fout(OUT_PATH, std::ios::trunc);
        if (fout.is_open()) {
            fout << "==========================================================\n"
                 << "  IBAGS (RLWE/NTRU) — Benchmark Run\n"
                 << "==========================================================\n"
                 << "  Demo  Params: n=" << g_demo_params.n
                 << " q=" << g_demo_params.q
                 << " sigma=" << (int)g_demo_params.sigma << "\n"
                 << "  NTT   Params: n=" << g_ntt_params.n
                 << " q=" << g_ntt_params.q << "\n";
            fout.close();
        }
    }

    /* §1  mod_reduce */
    bench_barrett_reduce();

    /* §2-6  poly */
    bench_poly_add();
    bench_poly_sub();
    bench_poly_mul_naive();
    bench_poly_norm_inf();
    bench_poly_inv();

    /* §7-10  poly_ntt */
    bench_ntt_forward();
    bench_ntt_inverse();
    bench_ntt_roundtrip();
    bench_ntt_mul();

    /* §11-12  gauss */
    bench_gauss_sample_coeff();
    bench_gauss_sample_poly();

    /* §13  xof */
    bench_xof_squeeze();

    /* §14-15  csprng */
    bench_csprng_system();
    bench_csprng_seed();

    /* §16-17  encode */
    bench_poly_encode();
    bench_poly_decode();

    /* §18-20  hash */
    bench_sample_ring_element();
    bench_sample_challenge_polynomial();
    bench_coeff_reject_small();

    /* §21  reject */
    bench_hash_to_uniform();

    /* §22-24  ntru */
    bench_ntru_secret_pair();
    bench_sample_f_until_invertible();
    bench_ntru_trapgen();

    std::cout << "\n============================================================\n"
              << "  IBAGS benchmark complete.\n"
              << "  Results written to benchmarking_output/benchmarking_ibags.txt\n"
              << "============================================================\n";

    return 0;
}
