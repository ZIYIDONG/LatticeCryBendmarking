/**
 * @file test_ntru_trapgen.cpp
 * @brief NTRU TrapGen 单元测试 — 四级参数全覆盖
 */

#include "../include_RLWEorNTRU/ntru_trapgen.h"
#include "../include_RLWEorNTRU/params.h"
#include "../include_RLWEorNTRU/poly.h"
#include "../include_RLWEorNTRU/csprng.h"
#include "../include_RLWEorNTRU/test_colors.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>
#include <array>

extern bool g_ibags_quiet;
extern bool g_ibags_nocolor;

static int g_passed = 0;
static int g_failed = 0;

#define C(suffix) (g_ibags_nocolor ? "" : COLOR_##suffix)

#define TEST(name) std::cout << "\n=== " << name << " ===\n"
#define CHECK(cond, msg) do { \
    if (!(cond)) { std::cerr << C(RED) << "  FAIL: " << msg << COLOR_RESET << std::endl; ++g_failed; } \
    else { if (!g_ibags_quiet) std::cout << C(GREEN) << "  PASS: " << msg << COLOR_RESET << std::endl; ++g_passed; } \
} while(0)

#define CHECK_OK(status, msg) do { \
    if (!(status).ok()) { std::cerr << C(RED) << "  FAIL: " << msg << " — " << (status).message() << COLOR_RESET << std::endl; ++g_failed; } \
    else { if (!g_ibags_quiet) std::cout << C(GREEN) << "  PASS: " << msg << COLOR_RESET << std::endl; ++g_passed; } \
} while(0)

#define CHECK_FAIL(status, msg) do { \
    if ((status).ok()) { std::cerr << C(RED) << "  FAIL: " << msg << " — expected error" << COLOR_RESET << std::endl; ++g_failed; } \
    else { if (!g_ibags_quiet) std::cout << C(GREEN) << "  PASS: " << msg << COLOR_RESET << std::endl; ++g_passed; } \
} while(0)

// ── 四级参数 ──
static std::array<ibags::Params, 4> all_params() {
    return {{
        ibags::Params::params_demo_64(),
        ibags::Params::params_level2_512(),
        ibags::Params::params_level3_1024(),
        ibags::Params::params_level5_1024(),
    }};
}

static const char* param_label(const ibags::Params& pp) {
    if (pp.n == 64)   return "Demo(n=64)";
    if (pp.n == 512)  return "L1(n=512)";
    if (pp.n == 1024 && pp.q == 16908289) return "L3(n=1024)";
    if (pp.n == 1024 && pp.q == 4206593)  return "L5(n=1024)";
    return "?";
}

static std::vector<uint8_t> make_seed(int id) {
    std::vector<uint8_t> seed(32, 0);
    seed[0] = static_cast<uint8_t>(id);
    std::snprintf(reinterpret_cast<char*>(seed.data() + 1), 31, "trapgen_%d", id);
    return seed;
}

// ============================================================================
// Test 1
// ============================================================================
void test_sample_f_until_invertible(const ibags::Params& pp) {
    TEST("sample_f_until_invertible – basic flow");
    auto seed = make_seed(1);
    ibags::SeedCSPRNG csprng(seed, "test_sample_f");
    ibags::Poly f_out(pp.n);
    CHECK_OK(ibags::sample_f_until_invertible(pp, csprng, 50, &f_out), "succeeds");
    CHECK(f_out.n() == pp.n, "output dim == pp.n");
    CHECK(ibags::poly_is_invertible_mod_q(f_out, pp), "f invertible");
    CHECK_FAIL(ibags::sample_f_until_invertible(pp, csprng, 10, nullptr), "null → error");
}

// ============================================================================
// Test 2
// ============================================================================
void test_solve_ntru_equation(const ibags::Params& pp) {
    TEST("solve_ntru_equation – stub behavior");
    ibags::NtruEquationResult result(ibags::Poly::zero(pp), ibags::Poly::zero(pp));
    auto st = ibags::solve_ntru_equation(pp, ibags::Poly::zero(pp), ibags::Poly::zero(pp), &result);
    CHECK_FAIL(st, "stub returns error");
    CHECK(st.message().find("STUB") != std::string::npos, "contains STUB");
    CHECK(result.F.n() == pp.n, "F dim correct");
    CHECK(result.G.n() == pp.n, "G dim correct");
}

// ============================================================================
// Test 3
// ============================================================================
void test_ntru_trapgen_default(const ibags::Params& pp) {
    TEST("ntru_trapgen – default config");
    auto seed = make_seed(2);
    ibags::SeedCSPRNG csprng(seed, "test_trapgen_d");
    auto config = ibags::TrapGenConfig::default_config();
    config.max_invertibility_attempts = 25;
    ibags::PublicTrapdoorParams pub(pp);
    ibags::MasterTrapdoorSecret sec(pp.n);
    CHECK_OK(ibags::ntru_trapgen(pp, csprng, config, &pub, &sec), "succeeds");
    CHECK(pub.h.n() == pp.n, "pub.h dim");
    if (pub.h_is_invertible) CHECK(pub.h_inv.n() == pp.n, "pub.h_inv dim");
    CHECK(sec.basis.dimensions_consistent(), "basis consistent");
    CHECK(ibags::check_h_from_f_g(pp, sec.basis.f, sec.basis.g, pub.h), "h==g*f^{-1}");
    CHECK(ibags::check_trapdoor_basis(pp, sec.basis, pub.h), "full check passes");
}

// ============================================================================
// Test 4
// ============================================================================
void test_ntru_trapgen_strict(const ibags::Params& pp) {
    TEST("ntru_trapgen – strict config");
    auto seed = make_seed(3);
    ibags::SeedCSPRNG csprng(seed, "test_trapgen_s");
    auto config = ibags::TrapGenConfig::strict_config();
    config.max_invertibility_attempts = 10;
    ibags::PublicTrapdoorParams pub(pp);
    ibags::MasterTrapdoorSecret sec(pp.n);
    CHECK_FAIL(ibags::ntru_trapgen(pp, csprng, config, &pub, &sec), "strict fails (stub)");
}

// ============================================================================
// Test 5
// ============================================================================
void test_gauss_sample_preimage(const ibags::Params& pp) {
    TEST("gauss_sample_preimage – stub");
    auto seed = make_seed(4);
    ibags::SeedCSPRNG csprng(seed, "test_gauss");
    ibags::TrapdoorBasis basis(pp.n);
    ibags::Poly s1(pp.n), s2(pp.n);
    auto st = ibags::gauss_sample_preimage(pp, basis, ibags::Poly::zero(pp), pp.sigma, csprng, &s1, &s2);
    CHECK_FAIL(st, "stub returns error");
    CHECK(st.message().find("STUB") != std::string::npos, "contains STUB");
    CHECK(s1.n() == pp.n, "s1 dim");
    CHECK(s2.n() == pp.n, "s2 dim");
}

// ============================================================================
// Test 6
// ============================================================================
void test_check_trapdoor_basis_edge_cases(const ibags::Params& pp) {
    TEST("check_trapdoor_basis – edge cases");
    ibags::TrapdoorBasis zero_basis(pp.n);
    ibags::Poly h_one(pp.n);
    for (int i = 0; i < pp.n; ++i) h_one[i] = 0;
    h_one[0] = 1; h_one.reduce_mod_q(pp);
    CHECK(!ibags::check_trapdoor_basis(pp, zero_basis, h_one), "zero basis + h=1 fails");

    ibags::Poly f2(pp.n), g2(pp.n), h2(pp.n);
    f2[0] = 1; g2[0] = 1; h2[0] = 1;
    CHECK(ibags::check_h_from_f_g(pp, f2, g2, h2), "f=1,g=1,h=1 passes");
}

// ============================================================================
// Test 7
// ============================================================================
void test_ntru_trapgen_reproducibility(const ibags::Params& pp) {
    TEST("ntru_trapgen – reproducibility");
    auto config = ibags::TrapGenConfig::default_config();
    config.max_invertibility_attempts = 25;
    auto s1 = make_seed(99), s2 = make_seed(99);
    ibags::SeedCSPRNG rng1(s1, "test_repro"), rng2(s2, "test_repro");
    ibags::PublicTrapdoorParams pub1(pp), pub2(pp);
    ibags::MasterTrapdoorSecret sec1(pp.n), sec2(pp.n);
    CHECK_OK(ibags::ntru_trapgen(pp, rng1, config, &pub1, &sec1), "first succeeds");
    CHECK_OK(ibags::ntru_trapgen(pp, rng2, config, &pub2, &sec2), "second succeeds");
    CHECK(ibags::poly_equal_ct(pub1.h, pub2.h), "h reproducible");
    CHECK(ibags::poly_equal_ct(sec1.basis.f, sec2.basis.f), "f reproducible");
}

// ============================================================================
// Test 8
// ============================================================================
void test_sample_f_max_attempts_boundary(const ibags::Params& pp) {
    TEST("sample_f_until_invertible – max_attempts boundary");
    auto seed = make_seed(5);
    ibags::SeedCSPRNG csprng(seed, "test_boundary");
    ibags::Poly f_out(pp.n);
    CHECK_FAIL(ibags::sample_f_until_invertible(pp, csprng, 0, &f_out), "max=0 → error");
    ibags::SeedCSPRNG csprng2(make_seed(6), "test_boundary2");
    auto st = ibags::sample_f_until_invertible(pp, csprng2, 1, &f_out);
    std::cout << "  INFO: max_attempts=1: " << (st.ok() ? "Ok" : st.message()) << std::endl;
    CHECK(true, "max=1 doesn't crash");
}

// ============================================================================
// Unified runner — 四级参数迭代
// ============================================================================
int run_ntru_trapgen_tests() {
    auto params = all_params();
    std::cout << "============================================\n"
              << " test_ntru_trapgen — NTRU TrapGen Unit Tests (4 levels)\n"
              << "============================================\n";

    for (const auto& pp : params) {
        std::cout << "\n========== " << param_label(pp)
                  << " n=" << pp.n << " q=" << pp.q << " ==========\n";
        test_sample_f_until_invertible(pp);
        test_solve_ntru_equation(pp);
        test_ntru_trapgen_default(pp);
        test_ntru_trapgen_strict(pp);
        test_gauss_sample_preimage(pp);
        test_check_trapdoor_basis_edge_cases(pp);
        test_ntru_trapgen_reproducibility(pp);
        test_sample_f_max_attempts_boundary(pp);
    }

    std::cout << "\n============================================\n";
    std::cout << " Results: " << g_passed << " passed, " << g_failed << " failed\n";
    std::cout << "============================================\n";
    return (g_failed == 0) ? 0 : 1;
}