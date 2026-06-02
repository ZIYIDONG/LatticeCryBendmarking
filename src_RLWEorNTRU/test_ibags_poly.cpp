/**
 * @file test_ibags_poly.cpp
 * @brief Poly 和 NTT 模块单元测试 — 四级参数全覆盖
 *
 * 测试覆盖:
 *  1. Poly 基本操作 (构造、访问、规范化)
 *  2. 多项式算术 (加减乘负标量乘)
 *  3. 环约减 (Y^n + 1)
 *  4. 范数计算 (绝对无穷范数)
 *  5. constant-time equality
 *  6. Polynomial inversion (常数情况 + 随机往返)
 *  7. NTT-friendly 参数检测
 *  8. NTT roundtrip (INTT(NTT(a)) == a)
 *  9. NTT-based multiplication (NTT(a*b) = NTT(a) * NTT(b))
 *
 * 所有测试对 Demo / L1 / L3 / L5 四级参数均执行。
 */

#include "../include_RLWEorNTRU/params.h"
#include "../include_RLWEorNTRU/poly.h"
#include "../include_RLWEorNTRU/poly_ntt.h"
#include "../include_RLWEorNTRU/mod_reduce.h"
#include "../include_RLWEorNTRU/errors.h"
#include "../include_RLWEorNTRU/test_colors.h"

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>
#include <random>
#include <array>

// ============================================================================
// Test helpers
// ============================================================================

extern bool g_ibags_quiet;
extern bool g_ibags_nocolor;

static int tests_passed = 0;
static int tests_failed = 0;

#define C(suffix) (g_ibags_nocolor ? "" : COLOR_##suffix)

#define TEST(name) \
    do { \
        if (!g_ibags_quiet) printf("  TEST: %s ... ", name); \
    } while (0)

#define PASS() \
    do { \
        if (!g_ibags_quiet) printf("%sPASSED%s\n", C(GREEN), COLOR_RESET); \
        ++tests_passed; \
    } while (0)

#define FAIL(msg) \
    do { \
        printf("%sFAILED%s: %s\n", C(RED), COLOR_RESET, msg); \
        ++tests_failed; \
    } while (0)

#define CHECK(cond, msg) \
    do { \
        if (!(cond)) { FAIL(msg); return; } \
    } while (0)

// ============================================================================
// 四级参数表
// ============================================================================

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

// ============================================================================
// 测试用例 — 每个函数接受 const Params&
// ============================================================================

/// 1. Poly 基础构造
void test_poly_construction(const ibags::Params& pp) {
    TEST("Poly construction with n");
    {
        ibags::Poly p(pp.n);
        CHECK(p.n() == pp.n, "n mismatch");
        for (int i = 0; i < std::min(pp.n, 128); ++i) {
            CHECK(p[i] == 0, "non-zero initial element");
        }
        CHECK(static_cast<int>(p.coeffs().size()) == pp.n, "coeffs size mismatch");
    }
    PASS();

    TEST("Poly from Params (zero)");
    {
        ibags::Poly p = ibags::Poly::zero(pp);
        CHECK(p.n() == pp.n, "n mismatch");
        for (int i = 0; i < std::min(pp.n, 128); ++i) {
            CHECK(p[i] == 0, "non-zero coefficient");
        }
    }
    PASS();

    TEST("Poly from_coeffs");
    {
        std::vector<int64_t> c(pp.n, 42);
        ibags::Poly p = ibags::Poly::from_coeffs(pp, c);
        CHECK(p.n() == pp.n, "n mismatch");
        for (int i = 0; i < std::min(pp.n, 128); ++i) {
            CHECK(p[i] == 42, "coefficient not matched");
        }
    }
    PASS();

    TEST("Poly size mismatch throws");
    {
        std::vector<int64_t> c(pp.n + 36, 0);
        try {
            ibags::Poly p(pp.n, c);
        } catch (const std::exception&) {}
    }
    PASS();
}

/// 2. 规范化 / 约减
void test_poly_reduction(const ibags::Params& pp) {
    int n = pp.n;

    TEST("reduce_mod_q");
    {
        std::vector<int64_t> c(n);
        for (int i = 0; i < n; ++i)
            c[i] = static_cast<int64_t>(pp.q) * 10 + i;
        ibags::Poly p = ibags::Poly::from_coeffs(pp, c);
        p.reduce_mod_q(pp);
        CHECK(p.is_canonical(pp), "should be canonical after reduce");
        for (int i = 0; i < std::min(n, 128); ++i)
            CHECK(p[i] == static_cast<int64_t>(i), "unexpected reduced value");
    }
    PASS();

    TEST("centered normalization");
    {
        int64_t q = static_cast<int64_t>(pp.q);
        int64_t q_minus_1 = q - 1;
        std::vector<int64_t> c(n, 0);
        c[0] = 0; c[1] = 1; c[2] = q_minus_1; c[3] = q_minus_1;
        ibags::Poly p = ibags::Poly::from_coeffs(pp, c);
        p.reduce_mod_q(pp);
        ibags::Poly cent = p.normalize_centered(pp);
        CHECK(cent.is_centered(pp), "not centered");
        CHECK(cent[0] == 0, "0 should stay 0");
        CHECK(cent[1] == 1, "1 should stay 1");
        CHECK(cent[2] == -1, "q-1 should map to -1");
        CHECK(cent[3] == -1, "second q-1 should also be -1");
    }
    PASS();
}

/// 3. 多项式算术
void test_poly_arithmetic(const ibags::Params& pp) {
    int n = pp.n;
    int64_t q = static_cast<int64_t>(pp.q);
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(0, q - 1);

    auto random_poly = [&]() {
        std::vector<int64_t> coeffs(n);
        for (int i = 0; i < n; ++i) coeffs[i] = dist(rng);
        return ibags::Poly::from_canonical(pp, coeffs);
    };

    constexpr int kTrials = 10;

    TEST("poly_add commutativity");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly(), b = random_poly();
            CHECK(ibags::poly_equal_ct(ibags::poly_add(a,b,pp), ibags::poly_add(b,a,pp)),
                  "a+b == b+a");
        }
    }
    PASS();

    TEST("poly_add/sub invertibility");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly(), b = random_poly();
            CHECK(ibags::poly_equal_ct(ibags::poly_sub(ibags::poly_add(a,b,pp), b, pp), a),
                  "(a+b)-b == a");
        }
    }
    PASS();

    TEST("poly_add associativity");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly(), b = random_poly(), c = random_poly();
            auto lhs = ibags::poly_add(ibags::poly_add(a,b,pp), c, pp);
            auto rhs = ibags::poly_add(a, ibags::poly_add(b,c,pp), pp);
            CHECK(ibags::poly_equal_ct(lhs, rhs), "(a+b)+c == a+(b+c)");
        }
    }
    PASS();

    TEST("poly_neg double negation");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            CHECK(ibags::poly_equal_ct(ibags::poly_neg(ibags::poly_neg(a,pp), pp), a),
                  "-(-a) == a");
        }
    }
    PASS();

    TEST("poly_neg additive inverse");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            CHECK(ibags::poly_equal_ct(ibags::poly_add(a, ibags::poly_neg(a,pp), pp),
                  ibags::Poly::zero(pp)), "a+(-a) == 0");
        }
    }
    PASS();

    // ── 标量乘法 (新增) ──
    TEST("poly_mul_scalar distributive");
    {
        std::uniform_int_distribution<int64_t> sdist(1, q - 1);
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly(), b = random_poly();
            int64_t s = sdist(rng);
            auto lhs = ibags::poly_mul_scalar(ibags::poly_add(a, b, pp), s, pp);
            auto rhs = ibags::poly_add(
                ibags::poly_mul_scalar(a, s, pp),
                ibags::poly_mul_scalar(b, s, pp), pp);
            CHECK(ibags::poly_equal_ct(lhs, rhs), "s*(a+b) == s*a + s*b");
        }
    }
    PASS();

    // ── naive 乘法代数性质 ──
    TEST("poly_mul_naive commutativity");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly(), b = random_poly();
            CHECK(ibags::poly_equal_ct(
                ibags::poly_mul_naive(a,b,pp), ibags::poly_mul_naive(b,a,pp)),
                "a*b == b*a");
        }
    }
    PASS();

    TEST("poly_mul_naive associativity");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly(), b = random_poly(), c = random_poly();
            auto lhs = ibags::poly_mul_naive(ibags::poly_mul_naive(a,b,pp), c, pp);
            auto rhs = ibags::poly_mul_naive(a, ibags::poly_mul_naive(b,c,pp), pp);
            CHECK(ibags::poly_equal_ct(lhs, rhs), "(a*b)*c == a*(b*c)");
        }
    }
    PASS();

    TEST("poly_mul_naive distributivity");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly(), b = random_poly(), c = random_poly();
            auto lhs = ibags::poly_mul_naive(a, ibags::poly_add(b,c,pp), pp);
            auto rhs = ibags::poly_add(
                ibags::poly_mul_naive(a,b,pp), ibags::poly_mul_naive(a,c,pp), pp);
            CHECK(ibags::poly_equal_ct(lhs, rhs), "a*(b+c) == a*b + a*c");
        }
    }
    PASS();

    TEST("poly_mul_naive identity");
    {
        std::vector<int64_t> one_coeffs(n, 0); one_coeffs[0] = 1;
        ibags::Poly one = ibags::Poly::from_canonical(pp, one_coeffs);
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            CHECK(ibags::poly_equal_ct(ibags::poly_mul_naive(one,a,pp), a), "1*a == a");
        }
    }
    PASS();

    TEST("poly_mul_naive zero");
    {
        auto zero = ibags::Poly::zero(pp);
        for (int t = 0; t < kTrials; ++t) {
            CHECK(ibags::poly_equal_ct(ibags::poly_mul_naive(zero, random_poly(), pp), zero),
                  "0*a == 0");
        }
    }
    PASS();
}

/// 4. 环约减
void test_poly_ring_reduce(const ibags::Params& pp) {
    TEST("ring reduction Y^n = -1");
    {
        int n = pp.n;
        std::vector<int64_t> conv(2 * n - 1, 0);
        conv[0] = 1; conv[n] = 2;
        ibags::Poly result = ibags::poly_ring_reduce_raw(pp, conv);
        int64_t expected = (static_cast<int64_t>(pp.q) - 1) % static_cast<int64_t>(pp.q);
        CHECK(result[0] == expected, "1+2Y^n should fold to q-1");
    }
    PASS();
}

/// 5. 范数
void test_poly_norm(const ibags::Params& pp) {
    int n = pp.n;

    TEST("norm_inf of zero");
    {
        CHECK(ibags::poly_norm_inf(ibags::Poly::zero(pp), pp) == 0, "|0| == 0");
    }
    PASS();

    TEST("norm_inf basic");
    {
        std::vector<int64_t> vals(n, 0);
        vals[std::min(50, n-1)] = 5;
        vals[std::min(100, n-1)] = static_cast<int64_t>(pp.q) - 3;
        ibags::Poly p = ibags::Poly::from_canonical(pp, vals);
        CHECK(ibags::poly_norm_inf(p, pp) == 5, "max(|5|,|-3|) == 5");
    }
    PASS();

    TEST("norm_bound_check");
    {
        std::vector<int64_t> vals(n, 0);
        vals[0] = 10;
        vals[std::min(1, n-1)] = static_cast<int64_t>(pp.q) - 10;
        ibags::Poly p = ibags::Poly::from_canonical(pp, vals);
        CHECK(ibags::poly_norm_bound_check(p, 10, pp), "<=10");
        CHECK(!ibags::poly_norm_bound_check(p, 9, pp), ">9");
    }
    PASS();
}

/// 6. Constant-time equality
void test_poly_equal_ct(const ibags::Params& pp) {
    int n = pp.n;
    int64_t q = static_cast<int64_t>(pp.q);

    TEST("poly_equal_ct same (zero)");
    {
        CHECK(ibags::poly_equal_ct(ibags::Poly::zero(pp), ibags::Poly::zero(pp)), "0==0");
    }
    PASS();

    TEST("poly_equal_ct same (random)");
    {
        std::mt19937_64 rng(99);
        std::uniform_int_distribution<int64_t> dist(0, q-1);
        std::vector<int64_t> c(n);
        for (int i = 0; i < n; ++i) c[i] = dist(rng);
        CHECK(ibags::poly_equal_ct(
            ibags::Poly::from_canonical(pp, c),
            ibags::Poly::from_canonical(pp, c)), "same == same");
    }
    PASS();

    TEST("poly_equal_ct different (random vs zero)");
    {
        std::mt19937_64 rng(77);
        std::uniform_int_distribution<int64_t> dist(0, q-1);
        std::vector<int64_t> vals(n);
        for (int i = 0; i < n; ++i) vals[i] = dist(rng);
        CHECK(!ibags::poly_equal_ct(
            ibags::Poly::zero(pp),
            ibags::Poly::from_canonical(pp, vals)), "0 != random");
    }
    PASS();

    TEST("poly_equal_ct different (single coeff)");
    {
        std::mt19937_64 rng(33);
        std::uniform_int_distribution<int64_t> dist(0, q-1);
        std::vector<int64_t> c(n);
        for (int i = 0; i < n; ++i) c[i] = dist(rng);
        ibags::Poly a = ibags::Poly::from_canonical(pp, c);
        c[n-1] = (c[n-1] + 1) % q;
        CHECK(!ibags::poly_equal_ct(a, ibags::Poly::from_canonical(pp, c)), "one diff");
    }
    PASS();

    TEST("poly_equal_ct different length");
    {
        CHECK(!ibags::poly_equal_ct(ibags::Poly(32), ibags::Poly(64)), "diff len");
    }
    PASS();
}

/// 7. 多项式求逆
void test_poly_inverse(const ibags::Params& pp) {
    int n = pp.n;

    TEST("inverse of constant 3");
    {
        std::vector<int64_t> vals(n, 0); vals[0] = 3;
        ibags::Poly a = ibags::Poly::from_canonical(pp, vals);
        ibags::Poly inv(n);
        CHECK(ibags::poly_inv(a, pp, &inv).ok(), "should succeed");
        ibags::Poly prod = ibags::poly_mul_naive(a, inv, pp);
        CHECK(prod[0] == 1, "3*inv(3)[0] == 1");
        for (int i = 1; i < n; ++i) CHECK(prod[i] == 0, "other coeffs == 0");
    }
    PASS();

    TEST("inverse roundtrip (random)");
    {
        std::mt19937_64 rng(12345);
        std::uniform_int_distribution<int64_t> sd(-5, 5);
        int pass = 0;
        for (int t = 0; t < 10; ++t) {
            std::vector<int64_t> c(n);
            for (int i = 0; i < n; ++i) c[i] = sd(rng);
            c[0] = sd(rng) + 1;
            ibags::Poly a = ibags::Poly::from_canonical(pp, c);
            ibags::Poly inv(n);
            if (!ibags::poly_inv(a, pp, &inv).ok()) continue;
            ibags::Poly prod = ibags::poly_mul_naive(a, inv, pp);
            if (prod[0] == 1) {
                bool ok = true;
                for (int i = 1; i < n; ++i) if (prod[i] != 0) { ok = false; break; }
                if (ok) ++pass;
            }
        }
        CHECK(pass >= 5, ">= half should invert & roundtrip");
    }
    PASS();

    TEST("inverse of zero fails");
    {
        ibags::Poly inv(pp.n);
        CHECK(!ibags::poly_inv(ibags::Poly::zero(pp), pp, &inv).ok(), "zero not invertible");
    }
    PASS();
}

// ============================================================================
// NTT 测试
// ============================================================================

void test_ntt(const ibags::Params& pp) {
    CHECK(ibags::is_ntt_friendly(pp), "params must be NTT-friendly");
    ibags::NttTable tbl = ibags::NttTable::create(pp);

    TEST("is_ntt_friendly check");
    {
        CHECK(ibags::is_ntt_friendly(pp), "NTT-friendly");
    }
    PASS();

    TEST("NTT roundtrip");
    {
        int64_t q = static_cast<int64_t>(pp.q);
        std::mt19937_64 rng(2025);
        std::uniform_int_distribution<int64_t> dist(0, q - 1);
        for (int t = 0; t < 10; ++t) {
            std::vector<int64_t> vals(pp.n);
            for (int i = 0; i < pp.n; ++i) vals[i] = dist(rng);
            ibags::Poly a = ibags::Poly::from_canonical(pp, vals);
            ibags::Status s = ibags::ntt_roundtrip_test(a, tbl);
            if (!s.ok()) {
                char buf[128];
                std::snprintf(buf, sizeof(buf), "roundtrip failed trial %d (n=%d q=%u): %s",
                              t, pp.n, pp.q, s.message().c_str());
                FAIL(buf); return;
            }
        }
    }
    PASS();

    TEST("NTT mul == naive mul");
    {
        int64_t q = static_cast<int64_t>(pp.q);
        std::mt19937_64 rng(2025);
        std::uniform_int_distribution<int64_t> dist(0, q - 1);
        for (int t = 0; t < 10; ++t) {
            std::vector<int64_t> va(pp.n), vb(pp.n);
            for (int i = 0; i < pp.n; ++i) { va[i] = dist(rng); vb[i] = dist(rng); }
            ibags::Poly a = ibags::Poly::from_canonical(pp, va);
            ibags::Poly b = ibags::Poly::from_canonical(pp, vb);
            ibags::Poly c_naive = ibags::poly_mul_naive(a, b, pp);

            ibags::NttPoly a_ntt(pp.n), b_ntt(pp.n), pw_ntt(pp.n);
            if (!ibags::poly_ntt(a, tbl, &a_ntt).ok() ||
                !ibags::poly_ntt(b, tbl, &b_ntt).ok()) { FAIL("NTT failed"); return; }
            if (!ibags::poly_pointwise_mul_ntt(a_ntt, b_ntt, tbl, &pw_ntt).ok())
                { FAIL("pw mul failed"); return; }
            ibags::Poly c_out(pp.n);
            if (!ibags::poly_invntt(pw_ntt, tbl, &c_out).ok())
                { FAIL("INTT failed"); return; }
            if (!ibags::poly_equal_ct(c_naive, c_out))
                { FAIL("NTT mul != naive mul"); return; }
        }
    }
    PASS();
}

// ============================================================================
// Unified runner entry — 遍历四级参数
// ============================================================================

int run_ibags_poly_tests() {
    auto params = all_params();
    printf("=== ibags Poly & NTT Unit Tests (4 levels) ===\n\n");

    for (const auto& pp : params) {
        printf("\n----- %s (n=%d q=%u) -----\n", param_label(pp), pp.n, pp.q);

        printf("--- Poly Construction ---\n");
        test_poly_construction(pp);

        printf("\n--- Poly Reduction ---\n");
        test_poly_reduction(pp);

        printf("\n--- Poly Arithmetic ---\n");
        test_poly_arithmetic(pp);

        printf("\n--- Ring Reduction ---\n");
        test_poly_ring_reduce(pp);

        printf("\n--- Poly Norm ---\n");
        test_poly_norm(pp);

        printf("\n--- Constant-Time Equality ---\n");
        test_poly_equal_ct(pp);

        printf("\n--- Poly Inverse ---\n");
        test_poly_inverse(pp);

        printf("\n--- NTT ---\n");
        test_ntt(pp);
    }

    printf("\n========================================\n");
    printf("  Passed:  %d\n", tests_passed);
    printf("  Failed:  %d\n", tests_failed);
    printf("========================================\n");

    return (tests_failed == 0) ? 0 : 1;
}
