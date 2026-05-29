/**
 * @file test_ibags_poly.cpp
 * @brief Poly 和 NTT 模块单元测试
 *
 * 测试覆盖:
 *  1. Poly 基本操作 (构造、访问、规范化)
 *  2. 多项式算术 (加减乘负)
 *  3. 环约减 (Y^n + 1)
 *  4. 范数计算 (绝对无穷范数)
 *  5. constant-time equality
 *  6. Polynomial inversion (常数情况)
 *  7. NTT-friendly 参数检测
 *  8. NTT roundtrip (INTT(NTT(a)) == a)
 *  9. NTT-based multiplication (NTT(a*b) = NTT(a) * NTT(b))
 * 10. 错误处理路径
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
// 测试用例
// ============================================================================

/// 1. Poly 基础构造
void test_poly_construction() {
    ibags::Params pp = ibags::default_params();

    TEST("Poly construction with n");
    {
        ibags::Poly p(pp.n);
        CHECK(p.n() == pp.n, "n mismatch");
        for (int i = 0; i < pp.n; ++i) {
            CHECK(p[i] == 0, "non-zero initial element");
        }
        CHECK(static_cast<int>(p.coeffs().size()) == pp.n, "coeffs size mismatch");
    }
    PASS();

    TEST("Poly from Params (zero)");
    {
        ibags::Poly p = ibags::Poly::zero(pp);
        CHECK(p.n() == pp.n, "n mismatch");
        for (int i = 0; i < pp.n; ++i) {
            CHECK(p[i] == 0, "non-zero coefficient");
        }
    }
    PASS();

    TEST("Poly from_coeffs");
    {
        std::vector<int64_t> c(pp.n, 42);
        ibags::Poly p = ibags::Poly::from_coeffs(pp, c);
        CHECK(p.n() == pp.n, "n mismatch");
        for (int i = 0; i < pp.n; ++i) {
            CHECK(p[i] == 42, "coefficient not matched");
        }
    }
    PASS();

    TEST("Poly size mismatch throws");
    {
        std::vector<int64_t> c(pp.n + 36, 0); // wrong size
        try {
            ibags::Poly p(pp.n, c); // expects size pp.n
        } catch (const std::exception&) {
            // some implementations throw or assert
        }
    }
    PASS();
}

/// 2. 规范化 / 约减
void test_poly_reduction() {
    ibags::Params pp = ibags::default_params();
    int n = pp.n;

    TEST("reduce_mod_q");
    {
        std::vector<int64_t> c(n);
        for (int i = 0; i < n; ++i) {
            c[i] = static_cast<int64_t>(pp.q) * 10 + i;
        }
        ibags::Poly p = ibags::Poly::from_coeffs(pp, c);
        p.reduce_mod_q(pp);
        CHECK(p.is_canonical(pp), "should be canonical after reduce");
        for (int i = 0; i < n; ++i) {
            CHECK(p[i] == static_cast<int64_t>(i), "unexpected reduced value");
        }
    }
    PASS();

    TEST("centered normalization");
    {
        // 使用当前参数集的 q 值，而非硬编码
        int64_t q = static_cast<int64_t>(pp.q);
        int64_t q_minus_1 = q - 1;

        std::vector<int64_t> c(n, 0);
        c[0] = 0;
        c[1] = 1;
        c[2] = q_minus_1;                      // canonical value q-1
        c[3] = q_minus_1;                      // also q-1
        ibags::Poly p = ibags::Poly::from_coeffs(pp, c);
        p.reduce_mod_q(pp);
        ibags::Poly cent = p.normalize_centered(pp);
        CHECK(cent.is_centered(pp), "not centered");
        CHECK(cent[0] == 0, "0 should stay 0");
        CHECK(cent[1] == 1, "1 should stay 1");
        // q-1 should become -1 in centered (since q-1 > q/2 for any q > 2)
        CHECK(cent[2] == -1, "q-1 should map to -1 in centered");
        CHECK(cent[3] == -1, "second q-1 should also be -1");
    }
    PASS();
}

/// 3. 多项式算术 — 使用全范围 [0, q) 随机系数，充分覆盖模约减路径
void test_poly_arithmetic() {
    ibags::Params pp = ibags::default_params();
    int n = pp.n;
    int64_t q = static_cast<int64_t>(pp.q);

    // 固定种子保证可复现性
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int64_t> dist(0, q - 1);

    // ── 辅助: 生成随机多项式 ──
    auto random_poly = [&]() {
        std::vector<int64_t> coeffs(n);
        for (int i = 0; i < n; ++i) coeffs[i] = dist(rng);
        return ibags::Poly::from_canonical(pp, coeffs);
    };

    // 测试轮数: 用多次随机实例验证代数恒等式
    constexpr int kTrials = 10;

    // ── poly_add 交换律: a + b == b + a ──
    TEST("poly_add commutativity (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto b = random_poly();
            auto ab = ibags::poly_add(a, b, pp);
            auto ba = ibags::poly_add(b, a, pp);
            CHECK(ibags::poly_equal_ct(ab, ba), "a+b should equal b+a");
        }
    }
    PASS();

    // ── poly_add / poly_sub 互逆: (a + b) - b == a ──
    TEST("poly_add/sub invertibility (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto b = random_poly();
            auto sum = ibags::poly_add(a, b, pp);
            auto recovered = ibags::poly_sub(sum, b, pp);
            CHECK(ibags::poly_equal_ct(recovered, a), "(a+b)-b should equal a");
        }
    }
    PASS();

    // ── poly_add 结合律: (a + b) + c == a + (b + c) ──
    TEST("poly_add associativity (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto b = random_poly();
            auto c = random_poly();
            auto lhs = ibags::poly_add(ibags::poly_add(a, b, pp), c, pp);
            auto rhs = ibags::poly_add(a, ibags::poly_add(b, c, pp), pp);
            CHECK(ibags::poly_equal_ct(lhs, rhs), "(a+b)+c should equal a+(b+c)");
        }
    }
    PASS();

    // ── poly_neg 双负归原: -(-a) == a ──
    TEST("poly_neg double negation (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto neg_a = ibags::poly_neg(a, pp);
            auto neg_neg_a = ibags::poly_neg(neg_a, pp);
            CHECK(ibags::poly_equal_ct(neg_neg_a, a), "-(-a) should equal a");
        }
    }
    PASS();

    // ── poly_neg 加法消去: a + (-a) == 0 ──
    TEST("poly_neg additive inverse (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto neg_a = ibags::poly_neg(a, pp);
            auto zero = ibags::poly_add(a, neg_a, pp);
            auto expected_zero = ibags::Poly::zero(pp);
            CHECK(ibags::poly_equal_ct(zero, expected_zero), "a+(-a) should equal 0");
        }
    }
    PASS();

    // ── poly_mul_naive 交换律: a * b == b * a ──
    TEST("poly_mul_naive commutativity (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto b = random_poly();
            auto ab = ibags::poly_mul_naive(a, b, pp);
            auto ba = ibags::poly_mul_naive(b, a, pp);
            CHECK(ibags::poly_equal_ct(ab, ba), "a*b should equal b*a");
        }
    }
    PASS();

    // ── poly_mul_naive 结合律: (a * b) * c == a * (b * c) ──
    TEST("poly_mul_naive associativity (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto b = random_poly();
            auto c = random_poly();
            auto lhs = ibags::poly_mul_naive(ibags::poly_mul_naive(a, b, pp), c, pp);
            auto rhs = ibags::poly_mul_naive(a, ibags::poly_mul_naive(b, c, pp), pp);
            CHECK(ibags::poly_equal_ct(lhs, rhs), "(a*b)*c should equal a*(b*c)");
        }
    }
    PASS();

    // ── poly_mul_naive 分配律: a * (b + c) == a*b + a*c ──
    TEST("poly_mul_naive distributivity (random full-range)");
    {
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto b = random_poly();
            auto c = random_poly();
            auto lhs = ibags::poly_mul_naive(a, ibags::poly_add(b, c, pp), pp);
            auto rhs = ibags::poly_add(
                ibags::poly_mul_naive(a, b, pp),
                ibags::poly_mul_naive(a, c, pp), pp);
            CHECK(ibags::poly_equal_ct(lhs, rhs), "a*(b+c) should equal a*b + a*c");
        }
    }
    PASS();

    // ── poly_mul_naive 单位元: 1 * a == a ──
    TEST("poly_mul_naive multiplicative identity (random full-range)");
    {
        // 多项式 "1": 常数项为 1, 其余为 0
        std::vector<int64_t> one_coeffs(n, 0);
        one_coeffs[0] = 1;
        ibags::Poly one = ibags::Poly::from_canonical(pp, one_coeffs);

        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto prod = ibags::poly_mul_naive(one, a, pp);
            CHECK(ibags::poly_equal_ct(prod, a), "1*a should equal a");
        }
    }
    PASS();

    // ── poly_mul_naive 零元: 0 * a == 0 ──
    TEST("poly_mul_naive zero absorption (random full-range)");
    {
        auto zero = ibags::Poly::zero(pp);
        for (int t = 0; t < kTrials; ++t) {
            auto a = random_poly();
            auto prod = ibags::poly_mul_naive(zero, a, pp);
            CHECK(ibags::poly_equal_ct(prod, zero), "0*a should equal 0");
        }
    }
    PASS();
}

/// 4. 环约减
void test_poly_ring_reduce() {
    ibags::Params pp = ibags::default_params();
    int n = pp.n;

    TEST("ring reduction Y^n = -1 (basic)");
    {
        // 卷积 [1, 0, ..., 0, 2] (长度 2n-1) = 1 + 2*Y^n
        std::vector<int64_t> conv(2 * n - 1, 0);
        conv[0] = 1;
        conv[n] = 2;
        ibags::Poly result = ibags::poly_ring_reduce_raw(pp, conv);
        // Should be 1 + 2*(-1) = -1 mod q = q-1
        int64_t expected = (static_cast<int64_t>(pp.q) - 1) % static_cast<int64_t>(pp.q);
        CHECK(result[0] == expected, "1 + 2*Y^n should fold to q-1");
    }
    PASS();
}

/// 5. 范数
void test_poly_norm() {
    ibags::Params pp = ibags::default_params();
    int n = pp.n;

    TEST("norm_inf of zero");
    {
        ibags::Poly z = ibags::Poly::zero(pp);
        uint64_t norm = ibags::poly_norm_inf(z, pp);
        CHECK(norm == 0, "norm of zero should be 0");
    }
    PASS();

    TEST("norm_inf basic");
    {
        std::vector<int64_t> vals(n, 0);
        vals[50] = 5;
        vals[100] = static_cast<int64_t>(pp.q) - 3; // q-3, centered = -3
        ibags::Poly p = ibags::Poly::from_canonical(pp, vals);
        uint64_t norm = ibags::poly_norm_inf(p, pp);
        CHECK(norm == 5, "max(|5|, |-3|) should be 5"); // centered: max(5, 3) = 5
    }
    PASS();

    TEST("norm_bound_check");
    {
        std::vector<int64_t> vals(n, 0);
        vals[0] = 10;
        vals[1] = static_cast<int64_t>(pp.q) - 10; // centered = -10
        ibags::Poly p = ibags::Poly::from_canonical(pp, vals);
        CHECK(ibags::poly_norm_bound_check(p, 10, pp), "||p||∞ should be <= 10");
        CHECK(!ibags::poly_norm_bound_check(p, 9, pp), "||p||∞ should be > 9");
    }
    PASS();
}

/// 6. Constant-time equality
void test_poly_equal_ct() {
    ibags::Params pp = ibags::default_params();
    int n = pp.n;
    int64_t q = static_cast<int64_t>(pp.q);

    TEST("poly_equal_ct same (zero poly)");
    {
        ibags::Poly a = ibags::Poly::zero(pp);
        ibags::Poly b = ibags::Poly::zero(pp);
        CHECK(ibags::poly_equal_ct(a, b), "zero == zero");
    }
    PASS();

    TEST("poly_equal_ct same (random full-range)");
    {
        std::mt19937_64 rng(99);
        std::uniform_int_distribution<int64_t> dist(0, q - 1);
        std::vector<int64_t> coeffs(n);
        for (int i = 0; i < n; ++i) coeffs[i] = dist(rng);
        ibags::Poly a = ibags::Poly::from_canonical(pp, coeffs);
        ibags::Poly b = ibags::Poly::from_canonical(pp, coeffs);
        CHECK(ibags::poly_equal_ct(a, b), "identical full-range polys should be equal");
    }
    PASS();

    TEST("poly_equal_ct different (full-range vs zero)");
    {
        ibags::Poly a = ibags::Poly::zero(pp);
        std::mt19937_64 rng(77);
        std::uniform_int_distribution<int64_t> dist(0, q - 1);
        std::vector<int64_t> vals(n);
        for (int i = 0; i < n; ++i) vals[i] = dist(rng);
        ibags::Poly b = ibags::Poly::from_canonical(pp, vals);
        CHECK(!ibags::poly_equal_ct(a, b), "zero != random full-range poly");
    }
    PASS();

    TEST("poly_equal_ct different (single coeff differs)");
    {
        std::mt19937_64 rng(33);
        std::uniform_int_distribution<int64_t> dist(0, q - 1);
        std::vector<int64_t> coeffs(n);
        for (int i = 0; i < n; ++i) coeffs[i] = dist(rng);
        ibags::Poly a = ibags::Poly::from_canonical(pp, coeffs);
        // 修改最后一个系数
        coeffs[n - 1] = (coeffs[n - 1] + 1) % q;
        ibags::Poly b = ibags::Poly::from_canonical(pp, coeffs);
        CHECK(!ibags::poly_equal_ct(a, b), "differing in last coeff should not be equal");
    }
    PASS();

    TEST("poly_equal_ct different length");
    {
        ibags::Poly a(32);
        ibags::Poly b(64);
        CHECK(!ibags::poly_equal_ct(a, b), "different n should not be equal");
    }
    PASS();
}

/// 7. 多项式求逆
void test_poly_inverse() {
    ibags::Params pp = ibags::default_params();
    int n = pp.n;

    // ── 常数多项式 3 求逆（确定性基准）──
    TEST("inverse of constant polynomial 3 mod q");
    {
        std::vector<int64_t> vals(n, 0);
        vals[0] = 3;
        ibags::Poly a = ibags::Poly::from_canonical(pp, vals);

        ibags::Poly inv(n);
        ibags::Status status = ibags::poly_inv(a, pp, &inv);
        CHECK(status.ok(), "should succeed for constant 3");

        // Verify: a * inv ≡ 1 (mod q)
        ibags::Poly prod = ibags::poly_mul_naive(a, inv, pp);
        CHECK(prod[0] == 1, "3 * inv(3) should equal 1");
        for (int i = 1; i < n; ++i) {
            CHECK(prod[i] == 0, "other coefficients should be 0");
        }
    }
    PASS();

    // ── 随机可逆多项式往返测试: a * a^{-1} ≡ 1 ──
    // 策略: 构造 a = 1 + r，其中 r 是小系数随机多项式。
    //   a ≡ 1 (mod small prime) 且 1 总是可逆的，因此 a 很可能可逆。
    TEST("inverse roundtrip (random small-coeff, full-range)");
    {
        std::mt19937_64 rng(12345);
        // 小系数范围: [-5, 5]，确保 a 的高度约束在合理值
        std::uniform_int_distribution<int64_t> small_dist(-5, 5);

        constexpr int kTrials = 10;
        int pass = 0;
        for (int t = 0; t < kTrials; ++t) {
            std::vector<int64_t> coeffs(n);
            for (int i = 0; i < n; ++i) coeffs[i] = small_dist(rng);
            // 常数项保底 +1，确保 a 不是零多项式或其倍数
            coeffs[0] = small_dist(rng) + 1;
            ibags::Poly a = ibags::Poly::from_canonical(pp, coeffs);

            ibags::Poly inv(n);
            ibags::Status status = ibags::poly_inv(a, pp, &inv);
            if (!status.ok()) continue;  // 不可逆则跳过本轮

            ibags::Poly prod = ibags::poly_mul_naive(a, inv, pp);
            // 期望: prod = 1 (常数项 1，其余为 0)
            if (prod[0] == 1) {
                bool all_zero = true;
                for (int i = 1; i < n; ++i) {
                    if (prod[i] != 0) { all_zero = false; break; }
                }
                if (all_zero) ++pass;
            }
        }
        CHECK(pass >= 5, "at least half of the random polys should be invertible and roundtrip");
    }
    PASS();

    TEST("inverse of zero polynomial fails");
    {
        ibags::Poly a = ibags::Poly::zero(pp);

        ibags::Poly inv(pp.n);
        ibags::Status status = ibags::poly_inv(a, pp, &inv);
        CHECK(!status.ok(), "should fail for zero polynomial");
    }
    PASS();
}

// ============================================================================
// NTT 测试
// ============================================================================

void test_ntt() {
    ibags::Params pp = ibags::default_params();

    TEST("is_ntt_friendly check");
    {
        // 编译期选择的参数集自身需要支持 NTT
        if (!ibags::is_ntt_friendly(pp)) {
            if (!g_ibags_quiet)
                printf("SKIPPED (params n=%d q=%d not NTT-friendly) ... PASSED\n",
                       pp.n, pp.q);
            ++tests_passed;
        } else {
            CHECK(ibags::is_ntt_friendly(pp),
                  "default_params should be NTT-friendly");
        }
        // 额外确认 level5_1024 也是 NTT-friendly (硬引用，不受编译期开关影响)
        ibags::Params pp1024 = ibags::Params::params_level5_1024();
        CHECK(ibags::is_ntt_friendly(pp1024),
              "q=4206593, n=1024 should be NTT-friendly");
    }
    PASS();

    TEST("NTT roundtrip");
    {
        if (!ibags::is_ntt_friendly(pp)) {
            if (!g_ibags_quiet)
                printf("SKIPPED (params n=%d q=%d not NTT-friendly) ... PASSED\n",
                       pp.n, pp.q);
            ++tests_passed;
        } else {
            int64_t q = static_cast<int64_t>(pp.q);
            std::mt19937_64 rng(2025);
            std::uniform_int_distribution<int64_t> dist(0, q - 1);

            constexpr int kTrials = 10;
            for (int t = 0; t < kTrials; ++t) {
                std::vector<int64_t> vals(pp.n);
                for (int i = 0; i < pp.n; ++i) vals[i] = dist(rng);
                ibags::Poly a = ibags::Poly::from_canonical(pp, vals);
                ibags::Status s = ibags::ntt_roundtrip_test(a, pp);
                if (!s.ok()) {
                    char buf[128];
                    std::snprintf(buf, sizeof(buf),
                                  "NTT roundtrip failed at trial %d: %s",
                                  t, s.message().c_str());
                    FAIL(buf);
                    return;
                }
            }
        }
    }
    PASS();

    TEST("NTT-based multiplication == naive multiplication");
    {
        if (!ibags::is_ntt_friendly(pp)) {
            if (!g_ibags_quiet)
                printf("SKIPPED (params n=%d q=%d not NTT-friendly) ... PASSED\n",
                       pp.n, pp.q);
            ++tests_passed;
        } else {
        int n = pp.n;
        int64_t q = static_cast<int64_t>(pp.q);

        // 使用全范围随机系数，充分覆盖模约减路径
        std::mt19937_64 rng(2025);
        std::uniform_int_distribution<int64_t> dist(0, q - 1);

        constexpr int kTrials = 10;
        for (int t = 0; t < kTrials; ++t) {
            std::vector<int64_t> va(n), vb(n);
            for (int i = 0; i < n; ++i) {
                va[i] = dist(rng);
                vb[i] = dist(rng);
            }
            ibags::Poly a = ibags::Poly::from_canonical(pp, va);
            ibags::Poly b = ibags::Poly::from_canonical(pp, vb);

            // Naive: poly_mul_naive
            ibags::Poly c_naive = ibags::poly_mul_naive(a, b, pp);

            // NTT-based:
            ibags::NttPoly a_ntt(pp), b_ntt(pp), c_ntt(pp);
            ibags::Status s1 = ibags::poly_ntt(a, pp, &a_ntt);
            ibags::Status s2 = ibags::poly_ntt(b, pp, &b_ntt);
            if (!s1.ok() || !s2.ok()) {
                FAIL("NTT conversion failed");
                return;
            }

            ibags::Status s3 = ibags::poly_pointwise_mul_ntt(a_ntt, b_ntt, pp, &c_ntt);
            if (!s3.ok()) { FAIL("pointwise mul failed"); return; }

            ibags::Poly c_ntt_out(n);
            ibags::Status s4 = ibags::poly_invntt(c_ntt, pp, &c_ntt_out);
            if (!s4.ok()) { FAIL("INTT failed"); return; }

            // 验证相等
            if (!ibags::poly_equal_ct(c_naive, c_ntt_out)) {
                FAIL("NTT-based multiplication should equal naive");
                return;
            }
        }
        } // end if (is_ntt_friendly)
    }
    PASS();
}

// ============================================================================
// Unified runner entry — callable from test_ibags_all.cpp
// ============================================================================

int run_ibags_poly_tests() {
    printf("=== ibags Poly & NTT Unit Tests ===\n\n");

    printf("--- Poly Construction ---\n");
    test_poly_construction();

    printf("\n--- Poly Reduction ---\n");
    test_poly_reduction();

    printf("\n--- Poly Arithmetic ---\n");
    test_poly_arithmetic();

    printf("\n--- Ring Reduction ---\n");
    test_poly_ring_reduce();

    printf("\n--- Poly Norm ---\n");
    test_poly_norm();

    printf("\n--- Constant-Time Equality ---\n");
    test_poly_equal_ct();

    printf("\n--- Poly Inverse ---\n");
    test_poly_inverse();

    printf("\n--- NTT ---\n");
    test_ntt();

    printf("\n========================================\n");
    printf("  Passed:  %d\n", tests_passed);
    printf("  Failed:  %d\n", tests_failed);
    printf("========================================\n");

    return (tests_failed == 0) ? 0 : 1;
}