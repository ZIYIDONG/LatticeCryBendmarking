/**
 * security_plain-LWE.cpp — SIS/LWE security parameter analysis
 *
 * Computes conservative estimates for:
 *   ① SIS hardness (Hermite factor, BKZ block size)
 *   ② LWE hardness (dual attack cost)
 *   ③ Correctness threshold (q/4 vs noise budget)
 */

#include "../include_plain-LWE/mp12_plain-LWE.h"
#include "../include_plain-LWE/unified_params_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <sstream>

/* ───── 文件输出 ───── */
static std::ostringstream sec_oss;
static void write_to_bench_file() {
    constexpr const char* OUT_PATH = "../bendmarking_output/bendmarking_plain-LWE.txt";
    std::ofstream fout(OUT_PATH, std::ios::app);
    if (fout.is_open()) { fout << sec_oss.str(); fout.close(); }
}

/* ═══════════════════════════════════════════════════
   §1  Hermite factor → BKZ block size lookup
   ═══════════════════════════════════════════════════ */
static double hermite_factor(int n, double norm, long q) {
    // Hermite delta^n ≈ ||v|| / q^(n/m)
    // For SIS we approximate: delta = (norm / q^(n/m))^(1/n)
    (void)n; (void)norm; (void)q;
    return 0.0; // placeholder
}

/* ═══════════════════════════════════════════════════
   §2  SIS hardness estimate
   ═══════════════════════════════════════════════════ */
struct SISEstimate {
    int n; long q; int m;
    double norm_bound;        // expected preimage norm (from sample_pre)
    double root_hermite;      // root-Hermite factor δ
    int bkz_blocksize;        // estimated BKZ block size
    int bits_security;        // estimated classical bit-security
};

static SISEstimate estimate_sis(const mp12::Params& p) {
    SISEstimate e;
    e.n = p.n; e.q = p.q; e.m = p.m;
    e.norm_bound = p.s * std::sqrt((double)p.m);  // approximate ||x|| ≤ s√m

    // Root-Hermite factor: δ = (||v|| / q^(n/m))^(1/n)
    double vol = std::pow((double)p.q, (double)p.n / p.m);
    e.root_hermite = std::pow(e.norm_bound / vol, 1.0 / p.n);

    // Rough BKZ block-size from δ ≈ ((β·π)^(1/β) · β/(2πe))^(1/(2β−2))
    // Simplified: β ≈ 0.009/ln(δ)^2  (experimental fit for δ ∈ [1.005, 1.02])
    if (e.root_hermite > 1.0001) {
        double ln_delta = std::log(e.root_hermite);
        e.bkz_blocksize = (int)(0.009 / (ln_delta * ln_delta) + 0.5);
        if (e.bkz_blocksize < 3) e.bkz_blocksize = 3;
        if (e.bkz_blocksize > 1000) e.bkz_blocksize = 1000;
    } else {
        e.bkz_blocksize = 1;
    }

    // Classical bit-security: ≈ 0.292·β  (BKZ core-SVP cost)
    e.bits_security = (int)(0.292 * e.bkz_blocksize + 0.5);

    return e;
}

/* ═══════════════════════════════════════════════════
   §3  LWE dual attack estimate
   ═══════════════════════════════════════════════════ */
struct LWEEstimate {
    int n; long q; double sigma;
    double alpha;             // noise rate σ/q
    double root_hermite;
    int bkz_blocksize;
    int bits_security;
};

static LWEEstimate estimate_lwe(const mp12::Params& p) {
    LWEEstimate e;
    e.n = p.n; e.q = p.q; e.sigma = 3.2; // LWE noise (unienc)
    e.alpha = e.sigma / e.q;

    // Dual attack: δ ≈ exp((α·q)^2 / (2·n))  (simplified)
    double log_delta = (e.alpha * e.q) * (e.alpha * e.q) / (2.0 * e.n);
    e.root_hermite = std::exp(log_delta);

    double ln_delta2 = log_delta; // ≈ ln(δ)
    if (ln_delta2 > 1e-6) {
        e.bkz_blocksize = (int)(0.009 / (ln_delta2 * ln_delta2) + 0.5);
        if (e.bkz_blocksize < 3) e.bkz_blocksize = 3;
        if (e.bkz_blocksize > 1000) e.bkz_blocksize = 1000;
    } else {
        e.root_hermite = 1.0;
        e.bkz_blocksize = 1;
    }
    e.bits_security = (int)(0.292 * e.bkz_blocksize + 0.5);

    return e;
}

/* ═══════════════════════════════════════════════════
   §4  Noise budget analysis
   ═══════════════════════════════════════════════════ */
void run_test_security() {
    sec_oss.str(""); sec_oss.clear();
    sec_oss << "=== Benchmark: Security Estimates ===\n";

    std::cout << "==========================================================\n";
    std::cout << "  Security Estimates — SIS / LWE / Noise Budget\n";
    std::cout << "==========================================================\n";

    auto p = unified::default_mp12_params();

    /* ── SIS ── */
    auto sis = estimate_sis(p);
    std::cout << "\n--- SIS Hardness Estimate ---\n"
              << "  n=" << sis.n << ", q=" << sis.q << ", m=" << sis.m << "\n"
              << "  Preimage norm bound: " << std::fixed << std::setprecision(1)
              << sis.norm_bound << "\n"
              << "  Root-Hermite factor δ: " << std::setprecision(6) << sis.root_hermite << "\n"
              << "  Est. BKZ block-size: " << sis.bkz_blocksize << "\n"
              << "  Classical bit-security: " << sis.bits_security << " bits\n";
    sec_oss << "\n--- SIS ---\n"
            << "  n=" << sis.n << " q=" << sis.q << " m=" << sis.m
            << "  ||x||=" << sis.norm_bound
            << "  δ=" << std::fixed << std::setprecision(6) << sis.root_hermite
            << "  β=" << sis.bkz_blocksize
            << "  λ=" << sis.bits_security << " bits\n";

    /* ── LWE ── */
    auto lwe = estimate_lwe(p);
    std::cout << "\n--- LWE Hardness Estimate (Dual Attack) ---\n"
              << "  n=" << lwe.n << ", q=" << lwe.q
              << ", σ=" << lwe.sigma << ", α=" << std::setprecision(4) << lwe.alpha << "\n"
              << "  Root-Hermite factor δ: " << std::setprecision(6) << lwe.root_hermite << "\n"
              << "  Est. BKZ block-size: " << lwe.bkz_blocksize << "\n"
              << "  Classical bit-security: " << lwe.bits_security << " bits\n";
    sec_oss << "\n--- LWE ---\n"
            << "  n=" << lwe.n << " q=" << lwe.q
            << "  σ=" << lwe.sigma << " α=" << std::fixed << std::setprecision(4) << lwe.alpha
            << "  δ=" << std::setprecision(6) << lwe.root_hermite
            << "  β=" << lwe.bkz_blocksize
            << "  λ=" << lwe.bits_security << " bits\n";

    /* ── Correctness ── */
    double noise_bound = p.s * std::sqrt((double)p.m);  // ‖x‖ after SamplePre
    double threshold = p.q / 4.0;
    std::cout << "\n--- Correctness Threshold ---\n"
              << "  Noise bound ‖x‖: " << std::fixed << std::setprecision(1) << noise_bound << "\n"
              << "  Decryption threshold q/4: " << threshold << "\n"
              << "  Safety margin: " << std::setprecision(1) << (threshold / noise_bound) << "×\n"
              << "  Verdict: " << (noise_bound < threshold ? "NOISE < q/4 (correct)" : "WARN: noise > q/4") << "\n";
    sec_oss << "\n--- Correctness ---\n"
            << "  ‖x‖=" << noise_bound << " q/4=" << threshold
            << "  margin=" << std::setprecision(1) << (threshold / noise_bound) << "x\n";

    write_to_bench_file();
    std::cout << "\nDone.\n";
}
