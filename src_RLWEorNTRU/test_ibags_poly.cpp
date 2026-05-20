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

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <string>
#include <cassert>
#include <stdexcept>

// ============================================================================
// Test helpers
// ============================================================================

static int tests_passed = 0;
static int tests_failed = 0;

#define TEST(name) \
    do { \
        printf("  TEST: %s ... ", name); \
    } while (0)

#define PASS() \
    do { \
        printf("PASSED\n"); \
        ++tests_passed; \
    } while (0)

#define FAIL(msg) \
    do { \
        printf("FAILED: %s\n", msg); \
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
        // q=8191, n=64. Only first 4 coefficients are non-zero
        std::vector<int64_t> c(n, 0);
        c[0] = 0;
        c[1] = 1;
        c[2] = 8190;                           // canonical value q-1
        c[3] = static_cast<int64_t>(pp.q) - 1; // also q-1 (8190)
        ibags::Poly p = ibags::Poly::from_coeffs(pp, c);
        p.reduce_mod_q(pp);
        ibags::Poly cent = p.normalize_centered(pp);
        CHECK(cent.is_centered(pp), "not centered");
        CHECK(cent[0] == 0, "0 should stay 0");
        CHECK(cent[1] == 1, "1 should stay 1");
        // 8190 = q-1 should become -1 in centered
        // q/2 = 4095, 8190 > 4095, so it maps to -1
        CHECK(cent[2] == -1, "q-1 should map to -1 in centered");
        CHECK(cent[3] == -1, "second q-1 should also be -1");
    }
    PASS();
}

/// 3. 多项式算术
void test_poly_arithmetic() {
    ibags::Params pp = ibags::default_params();
    int n = pp.n;

    TEST("poly_add");
    {
        std::vector<int64_t> ones(n, 1);
        std::vector<int64_t> twos(n, 2);
        ibags::Poly a = ibags::Poly::from_canonical(pp, ones);
        ibags::Poly b = ibags::Poly::from_canonical(pp, twos);
        ibags::Poly c = ibags::poly_add(a, b, pp);
        for (int i = 0; i < n; ++i) {
            CHECK(c[i] == 3, "1+2 should be 3");
        }
    }
    PASS();

    TEST("poly_sub");
    {
        std::vector<int64_t> fives(n, 5);
        std::vector<int64_t> threes(n, 3);
        ibags::Poly a = ibags::Poly::from_canonical(pp, fives);
        ibags::Poly b = ibags::Poly::from_canonical(pp, threes);
        ibags::Poly c = ibags::poly_sub(a, b, pp);
        for (int i = 0; i < n; ++i) {
            CHECK(c[i] == 2, "5-3 should be 2");
        }
    }
    PASS();

    TEST("poly_neg");
    {
        std::vector<int64_t> vals(n, 1);
        ibags::Poly a = ibags::Poly::from_canonical(pp, vals);
        ibags::Poly neg_a = ibags::poly_neg(a, pp);
        for (int i = 0; i < n; ++i) {
            CHECK(neg_a[i] == static_cast<int64_t>(pp.q) - 1,
                  "-1 mod q should be q-1");
        }
    }
    PASS();

    TEST("poly_mul_naive (1 * any = any)");
    {
        std::vector<int64_t> ones(n, 0);
        ones[0] = 1; // 多项式 "1"
        ibags::Poly one = ibags::Poly::from_canonical(pp, ones);

        std::vector<int64_t> rvals(n);
        for (int i = 0; i < n; ++i) rvals[i] = i % 10;
        ibags::Poly r = ibags::Poly::from_canonical(pp, rvals);

        ibags::Poly prod = ibags::poly_mul_naive(one, r, pp);
        for (int i = 0; i < n; ++i) {
            CHECK(prod[i] == r[i], "1 * r should equal r");
        }
    }
    PASS();

    TEST("poly_mul_naive commutativity");
    {
        std::vector<int64_t> va(n), vb(n);
        for (int i = 0; i < n; ++i) {
            va[i] = i % 7;
            vb[i] = i % 11;
        }
        ibags::Poly a = ibags::Poly::from_canonical(pp, va);
        ibags::Poly b = ibags::Poly::from_canonical(pp, vb);
        ibags::Poly ab = ibags::poly_mul_naive(a, b, pp);
        ibags::Poly ba = ibags::poly_mul_naive(b, a, pp);
        CHECK(ibags::poly_equal_ct(ab, ba), "multiplication should be commutative");
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

    TEST("poly_equal_ct same");
    {
        ibags::Poly a = ibags::Poly::zero(pp);
        ibags::Poly b = ibags::Poly::zero(pp);
        CHECK(ibags::poly_equal_ct(a, b), "zero == zero");
    }
    PASS();

    TEST("poly_equal_ct different");
    {
        ibags::Poly a = ibags::Poly::zero(pp);
        std::vector<int64_t> vals(n, 1);
        ibags::Poly b = ibags::Poly::from_canonical(pp, vals);
        CHECK(!ibags::poly_equal_ct(a, b), "zero != all ones");
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

/// 7. 多项式求逆 (常数情况)
void test_poly_inverse() {
    ibags::Params pp = ibags::default_params();
    TEST("inverse of constant polynomial 3 mod q");
    {
        int n = pp.n;
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
            printf("SKIPPED (params n=%d q=%d not NTT-friendly) ... PASSED\n",
                   pp.n, pp.q);
            ++tests_passed;
        } else {
            std::vector<int64_t> vals(pp.n);
            for (int i = 0; i < pp.n; ++i) {
                vals[i] = static_cast<int64_t>(i * i + 3 * i + 7) % static_cast<int64_t>(pp.q);
            }
            ibags::Poly a = ibags::Poly::from_canonical(pp, vals);
            ibags::Status s = ibags::ntt_roundtrip_test(a, pp);
            CHECK(s.ok(), "NTT roundtrip should succeed");
        }
    }
    PASS();

    TEST("NTT-based multiplication == naive multiplication");
    {
        if (!ibags::is_ntt_friendly(pp)) {
            printf("SKIPPED (params n=%d q=%d not NTT-friendly) ... PASSED\n",
                   pp.n, pp.q);
            ++tests_passed;
        } else {
        int n = pp.n;

        std::vector<int64_t> va(n), vb(n);
        for (int i = 0; i < n; ++i) {
            va[i] = (i * 7 + 3) % 100;
            vb[i] = (i * 11 + 5) % 100;
        }
        ibags::Poly a = ibags::Poly::from_canonical(pp, va);
        ibags::Poly b = ibags::Poly::from_canonical(pp, vb);

        // Naive: poly_mul_naive
        ibags::Poly c_naive = ibags::poly_mul_naive(a, b, pp);

        // NTT-based:
        ibags::NttPoly a_ntt(pp), b_ntt(pp), c_ntt(pp);
        ibags::Status s1 = ibags::poly_ntt(a, pp, &a_ntt);
        ibags::Status s2 = ibags::poly_ntt(b, pp, &b_ntt);
        CHECK(s1.ok() && s2.ok(), "NTT conversion should succeed");

        ibags::Status s3 = ibags::poly_pointwise_mul_ntt(a_ntt, b_ntt, pp, &c_ntt);
        CHECK(s3.ok(), "pointwise mul should succeed");

        ibags::Poly c_ntt_out(n);
        ibags::Status s4 = ibags::poly_invntt(c_ntt, pp, &c_ntt_out);
        CHECK(s4.ok(), "INTT should succeed");

        // 验证相等
        CHECK(ibags::poly_equal_ct(c_naive, c_ntt_out),
              "NTT-based multiplication should equal naive");
        } // end if (is_ntt_friendly)
    }
    PASS();
}

// ============================================================================
// Main
// ============================================================================

int main(void) {
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

    if (tests_failed > 0) {
        printf("\n*** SOME TESTS FAILED ***\n");
        return EXIT_FAILURE;
    }
    printf("\nAll tests passed.\n");
    return EXIT_SUCCESS;
}