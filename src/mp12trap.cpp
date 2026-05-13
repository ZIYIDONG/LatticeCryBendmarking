/* mp12trap.cpp
 * Extracted trapgen-related tests and helpers from main.cpp
 */

#include "../include/mp12trap.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <algorithm>

using namespace mp12;

/* ─────────────── pretty-print helpers ──────────────── */
static void print_vec(const char* label, const Vec& v, int maxshow = 8) {
    std::cout << label << " [";
    int show = std::min((int)v.size(), maxshow);
    for (int i = 0; i < show; i++)
        std::cout << (i ? ", " : "") << v[i];
    if ((int)v.size() > show) std::cout << " ...";
    std::cout << "]  (dim=" << v.size() << ")\n";
}

static void print_mat(const char* label, const Mat& M,
                      int maxr = 4, int maxc = 8) {
    std::cout << label << "  (" << M.size() << "×" << M[0].size() << ")\n";
    int r = std::min((int)M.size(), maxr);
    int c = std::min((int)M[0].size(), maxc);
    for (int i = 0; i < r; i++) {
        std::cout << "  [";
        for (int j = 0; j < c; j++)
            std::cout << std::setw(6) << M[i][j];
        if ((int)M[0].size() > c) std::cout << " ...";
        std::cout << " ]\n";
    }
    if ((int)M.size() > r) std::cout << "  ...\n";
}

double mp12::vec_norm(const Vec& v) {
    double s = 0;
    for (auto x : v) s += (double)x * x;
    return std::sqrt(s);
}

/* ─────────────────────── tests ─────────────────────── */
void mp12::test_gadget(const Params& p) {
    std::cout << "\n=== Test 1: Gadget Matrix G ===\n";
    Mat G = gadget_matrix(p);
    print_mat("G", G);

    // Verify: G[i][i*k+j] = b^j
    long bpow = 1;
    bool ok = true;
    for (int j = 0; j < p.k; j++) {
        for (int i = 0; i < p.n; i++) {
            long expected = bpow % p.q;
            if (G[i][i * p.k + j] != expected) { ok = false; break; }
        }
        bpow *= p.b;
    }
    std::cout << "Gadget structure check: " << (ok ? "PASS" : "FAIL") << "\n";
}

void mp12::test_gadget_basis(const Params& p) {
    std::cout << "\n=== Test 2: Gadget Basis S_g ===\n";
    Mat Sg = gadget_basis(p);
    print_mat("S_g (first block)", Sg, p.k, p.k);

    // Verify: G · S_g = 0 (mod q)
    Mat G = gadget_matrix(p);
    Mat GS = mat_mul_mod(G, Sg, p.q);
    bool all_zero = true;
    for (int i = 0; i < p.n; i++)
        for (int j = 0; j < p.n * p.k; j++)
            if (GS[i][j] != 0) { all_zero = false; }
    std::cout << "G · S_g = 0 (mod q): " << (all_zero ? "PASS" : "FAIL") << "\n";
}

void mp12::test_sample_g(const Params& p) {
    std::cout << "\n=== Test 3: SampleG (G-lattice short preimage) ===\n";
    UniformSampler usampler(p.q, 42);
    int trials = 20, pass = 0;
    double avg_norm = 0;
    Mat G = gadget_matrix(p);
    for (int t = 0; t < trials; t++) {
        // random target
        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();
        Vec z = sample_g(p, u);
        Vec Gz = mat_vec_mod(G, z, p.q);
        bool ok = true;
        for (int i = 0; i < p.n; i++)
            if (Gz[i] != mod(u[i], p.q)) { ok = false; break; }
        if (ok) pass++;
        avg_norm += vec_norm(z);
    }
    avg_norm /= trials;
    std::cout << "SampleG correctness: " << pass << "/" << trials << " PASS\n";
    std::cout << "Average preimage norm: " << std::fixed << std::setprecision(2)
              << avg_norm << "\n";
}

void mp12::test_gen_trap(const Params& p) {
    std::cout << "\n=== Test 4: GenTrap ===\n";
    auto t0 = std::chrono::high_resolution_clock::now();
    Trapdoor td = gen_trap(p, 123);
    auto t1 = std::chrono::high_resolution_clock::now();
    double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    print_mat("A", td.A, 3, 10);
    print_mat("R (trapdoor)", td.R, 3, 8);
    std::cout << "GenTrap time: " << std::fixed << std::setprecision(2) << ms << " ms\n";

    // Verify: A · [R; I] = G (mod q)
    // Build T = [R; I_{nk}]
    int nk = p.n * p.k;
    Mat T = make_mat(p.m, nk);
    for (int i = 0; i < p.m_bar; i++)
        for (int j = 0; j < nk; j++)
            T[i][j] = td.R[i][j];
    for (int j = 0; j < nk; j++)
        T[p.m_bar + j][j] = 1;

    Mat AT = mat_mul_mod(td.A, T, p.q);
    Mat G  = gadget_matrix(p);
    bool ok = true;
    for (int i = 0; i < p.n && ok; i++)
        for (int j = 0; j < nk && ok; j++)
            if (AT[i][j] != G[i][j]) ok = false;
    std::cout << "A · [R; I] = G (mod q): " << (ok ? "PASS" : "FAIL") << "\n";
}

void mp12::test_sample_pre(const Params& p) {
    std::cout << "\n=== Test 5: SamplePre (full round-trip) ===\n";
    Trapdoor td = gen_trap(p, 999);
    UniformSampler usampler(p.q, 77);

    int trials = 10, pass = 0;
    double avg_norm = 0;
    for (int t = 0; t < trials; t++) {
        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();

        auto t0 = std::chrono::high_resolution_clock::now();
        Vec x = sample_pre(p, td, u, (uint64_t)(t * 100 + 1));
        auto t1 = std::chrono::high_resolution_clock::now();

        bool ok = verify(p, td.A, x, u);
        if (ok) pass++;
        double nrm = vec_norm(x);
        avg_norm += nrm;

        if (t < 3) {
            std::cout << "  Trial " << t << ": ||x||=" << std::setprecision(1)
                      << std::fixed << nrm
                      << "  A·x=u? " << (ok ? "YES" : "NO")
                      << "  time=" << std::chrono::duration<double, std::milli>(t1-t0).count()
                      << "ms\n";
        }
    }
    avg_norm /= trials;
    std::cout << "SamplePre correctness: " << pass << "/" << trials << " PASS\n";
    std::cout << "Average preimage norm: " << std::fixed << std::setprecision(2)
              << avg_norm << "\n";
    std::cout << "Gaussian width s = " << p.s << "\n";
}

void mp12::test_uniformity(const Params& p) {
    std::cout << "\n=== Test 6: A distribution (uniformity check) ===\n";
    // Generate several As and check column averages ≈ q/2
    int num = 5;
    double avg = 0;
    for (int i = 0; i < num; i++) {
        Trapdoor td = gen_trap(p, (uint64_t)i * 31337);
        for (int r = 0; r < p.n; r++)
            for (int c = 0; c < p.m; c++)
                avg += td.A[r][c];
    }
    avg /= (double)(num * p.n * p.m);
    double expected = (p.q - 1.0) / 2.0;
    double rel_err = std::abs(avg - expected) / expected * 100.0;
    std::cout << "Mean entry of A: " << std::fixed << std::setprecision(1) << avg
              << "  (expected ≈ " << expected << ")  relative error: "
              << std::setprecision(2) << rel_err << "%\n";
    std::cout << "Uniformity check: " << (rel_err < 5.0 ? "PASS" : "WARN") << "\n";
}

void mp12::test_full_roundtrip_large(const Params& p2) {
    std::cout << "\n\n=== Larger params: n=16, q=8209 ===\n";
    std::cout << "  m = " << p2.m << ",  k = " << p2.k << "\n";
    Trapdoor td2 = gen_trap(p2, 42);
    UniformSampler us2(p2.q, 13);
    Vec u2(p2.n);
    for (int i = 0; i < p2.n; i++) u2[i] = us2.sample();
    Vec x2 = sample_pre(p2, td2, u2, 17);
    bool ok2 = verify(p2, td2.A, x2, u2);
    std::cout << "  SamplePre correctness: " << (ok2 ? "PASS" : "FAIL") << "\n";
    std::cout << "  ||x|| = " << std::fixed << std::setprecision(1) << vec_norm(x2) << "\n";
}

void mp12::run_mp12_trap_tests(const Params& p) {
    test_gadget(p);
    test_gadget_basis(p);
    test_sample_g(p);
    test_gen_trap(p);
    test_sample_pre(p);
    test_uniformity(p);
    test_full_roundtrip_large(p);
}
