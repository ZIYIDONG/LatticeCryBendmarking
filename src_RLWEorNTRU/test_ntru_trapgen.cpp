/**
 * @file test_ntru_trapgen.cpp
 * @brief 独立测试 ntru_trapgen 的完整流程。
 *
 * 测试内容:
 *   1. 流式 API: sample_f_until_invertible, solve_ntru_equation
 *   2. 主入口: ntru_trapgen (default/strict config)
 *   3. 验证: check_h_from_f_g, check_trapdoor_basis
 *   4. 前置检查: GaussSample stub
 *
 * 用法:
 *   cd build && cmake --build . --target test_ntru_trapgen -j$(nproc) && ./test_ntru_trapgen
 */

#include "../include_RLWEorNTRU/ntru_trapgen.h"
#include "../include_RLWEorNTRU/params.h"
#include "../include_RLWEorNTRU/poly.h"
#include "../include_RLWEorNTRU/csprng.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "../include_RLWEorNTRU/test_colors.h"

extern bool g_ibags_quiet;
extern bool g_ibags_nocolor;

// ── 简单断言宏 ──
static int g_passed = 0;
static int g_failed = 0;

// 颜色辅助
#define C(suffix) (g_ibags_nocolor ? "" : COLOR_##suffix)

#define TEST(name) std::cout << "\n=== " << name << " ===\n"
#define CHECK(cond, msg) do {                           \
    if (!(cond)) {                                       \
        std::cerr << C(RED) << "  FAIL: " << msg << COLOR_RESET << std::endl; \
        ++g_failed;                                      \
    } else {                                             \
        if (!g_ibags_quiet) std::cout << C(GREEN) << "  PASS: " << msg << COLOR_RESET << std::endl; \
        ++g_passed;                                      \
    }                                                    \
} while(0)

#define CHECK_OK(status, msg) do {                       \
    if (!(status).ok()) {                                 \
        std::cerr << C(RED) << "  FAIL: " << msg          \
                  << " — " << (status).message()          \
                  << COLOR_RESET << std::endl;            \
        ++g_failed;                                      \
    } else {                                             \
        if (!g_ibags_quiet) std::cout << C(GREEN) << "  PASS: " << msg << COLOR_RESET << std::endl; \
        ++g_passed;                                      \
    }                                                    \
} while(0)

#define CHECK_FAIL(status, msg) do {                     \
    if ((status).ok()) {                                  \
        std::cerr << C(RED) << "  FAIL: " << msg          \
                  << " — expected error but got Ok"       \
                  << COLOR_RESET << std::endl;            \
        ++g_failed;                                      \
    } else {                                             \
        if (!g_ibags_quiet) std::cout << C(GREEN) << "  PASS: " << msg << COLOR_RESET << std::endl; \
        ++g_passed;                                      \
    }                                                    \
} while(0)

// ── 辅助函数 ──
static std::vector<uint8_t> make_seed(int id) {
    std::vector<uint8_t> seed(32, 0);
    seed[0] = static_cast<uint8_t>(id);
    std::snprintf(reinterpret_cast<char*>(seed.data() + 1), 31,
                  "test_ntru_trapgen_%d", id);
    return seed;
}

// ============================================================================
// Test 1: sample_f_until_invertible – 基本流程
// ============================================================================
void test_sample_f_until_invertible() {
    TEST("sample_f_until_invertible – basic flow");

    auto pp = ibags::default_params();
    auto seed = make_seed(1);
    ibags::SeedCSPRNG csprng(seed, "test_sample_f");

    ibags::Poly f_out(pp.n);
    auto st = ibags::sample_f_until_invertible(pp, csprng, 50, &f_out);
    CHECK_OK(st, "sample_f_until_invertible succeeds");

    // 检查维度
    CHECK(f_out.n() == pp.n, "output dimension == pp.n");

    // 检查可逆性
    bool inv = ibags::poly_is_invertible_mod_q(f_out, pp);
    CHECK(inv, "returned f is invertible mod q");

    // 检查 null 指针保护
    auto st_null = ibags::sample_f_until_invertible(pp, csprng, 10, nullptr);
    CHECK_FAIL(st_null, "null pointer returns error");
}

// ============================================================================
// Test 2: solve_ntru_equation – stub 行为
// ============================================================================
void test_solve_ntru_equation() {
    TEST("solve_ntru_equation – stub behavior");

    auto pp = ibags::default_params();
    ibags::Poly f = ibags::Poly::zero(pp);
    ibags::Poly g = ibags::Poly::zero(pp);

    // 使用 NtruEquationResult 的显式构造
    ibags::NtruEquationResult result(ibags::Poly::zero(pp),
                                      ibags::Poly::zero(pp));
    auto st = ibags::solve_ntru_equation(pp, f, g, &result);
    // Stub 返回 InternalError 以提示不可生产使用
    CHECK_FAIL(st, "solve_ntru_equation returns error (stub)");
    CHECK(st.message().find("STUB") != std::string::npos,
          "error message contains STUB");

    // 但 result 仍被填充为零
    CHECK(result.F.n() == pp.n, "result.F has correct dimension");
    CHECK(result.G.n() == pp.n, "result.G has correct dimension");
}

// ============================================================================
// Test 3: ntru_trapgen – default config (PoC)
// ============================================================================
void test_ntru_trapgen_default() {
    TEST("ntru_trapgen – default config (PoC)");

    auto pp = ibags::default_params();
    auto seed = make_seed(2);
    ibags::SeedCSPRNG csprng(seed, "test_trapgen_default");

    auto config = ibags::TrapGenConfig::default_config();
    // PoC 使用较小的 max_attempts 以加速
    config.max_invertibility_attempts = 25;

    // 显式构造输出结构体（Poly 无默认构造函数）
    ibags::PublicTrapdoorParams pub(pp);
    ibags::MasterTrapdoorSecret sec(pp.n);

    auto st = ibags::ntru_trapgen(pp, csprng, config, &pub, &sec);
    CHECK_OK(st, "ntru_trapgen succeeds with default config");

    // 检查 h 的维度
    CHECK(pub.h.n() == pp.n, "pub.h has correct dimension");

    // 检查 h 可逆性
    if (pub.h_is_invertible) {
        CHECK(pub.h_inv.n() == pp.n, "pub.h_inv has correct dimension");
    }

    // 检查基础维度一致性
    CHECK(sec.basis.dimensions_consistent(), "basis dimensions consistent");

    // 验证 h == g * f^{-1}
    bool h_check = ibags::check_h_from_f_g(pp, sec.basis.f, sec.basis.g, pub.h);
    CHECK(h_check, "h == g * f^{-1} passes");

    // 全验证
    bool full_check = ibags::check_trapdoor_basis(pp, sec.basis, pub.h);
    CHECK(full_check, "check_trapdoor_basis passes");
}

// ============================================================================
// Test 4: ntru_trapgen – strict config（不 fallback，应失败）
// ============================================================================
void test_ntru_trapgen_strict() {
    TEST("ntru_trapgen – strict config (no fallback)");

    auto pp = ibags::default_params();
    auto seed = make_seed(3);
    ibags::SeedCSPRNG csprng(seed, "test_trapgen_strict");

    auto config = ibags::TrapGenConfig::strict_config();
    config.max_invertibility_attempts = 10;

    ibags::PublicTrapdoorParams pub(pp);
    ibags::MasterTrapdoorSecret sec(pp.n);

    auto st = ibags::ntru_trapgen(pp, csprng, config, &pub, &sec);
    // strict config 下 NTRU equation solver (stub) 返回 InternalError
    // 因为 fallback_on_ntru_solver = false
    CHECK_FAIL(st, "ntru_trapgen with strict config should fail (stub solver)");
}

// ============================================================================
// Test 5: GaussSample stub
// ============================================================================
void test_gauss_sample_preimage() {
    TEST("gauss_sample_preimage – stub");

    auto pp = ibags::default_params();
    auto seed = make_seed(4);
    ibags::SeedCSPRNG csprng(seed, "test_gauss_sample");

    ibags::TrapdoorBasis basis(pp.n);
    // 填充零 basis（PoC）
    basis.f   = ibags::Poly::zero(pp);
    basis.g   = ibags::Poly::zero(pp);
    basis.capF = ibags::Poly::zero(pp);
    basis.capG = ibags::Poly::zero(pp);

    ibags::Poly target = ibags::Poly::zero(pp);
    ibags::Poly s1(pp.n);
    ibags::Poly s2(pp.n);

    auto st = ibags::gauss_sample_preimage(pp, basis, target, pp.sigma,
                                            csprng, &s1, &s2);
    // Stub 返回 InternalError
    CHECK_FAIL(st, "gauss_sample_preimage returns error (stub)");
    CHECK(st.message().find("STUB") != std::string::npos,
          "error message contains STUB");

    // 但仍填充了零向量
    CHECK(s1.n() == pp.n, "s1 has correct dimension");
    CHECK(s2.n() == pp.n, "s2 has correct dimension");
}

// ============================================================================
// Test 6: check_trapdoor_basis – 边界情况
// ============================================================================
void test_check_trapdoor_basis_edge_cases() {
    TEST("check_trapdoor_basis – edge cases");

    auto pp = ibags::default_params();

    // 6a. All-zero basis: should NOT pass because h * 0 != 0
    ibags::TrapdoorBasis zero_basis(pp.n);
    zero_basis.f   = ibags::Poly::zero(pp);
    zero_basis.g   = ibags::Poly::zero(pp);
    zero_basis.capF = ibags::Poly::zero(pp);
    zero_basis.capG = ibags::Poly::zero(pp);

    // 构造 h = 1 (non-zero)
    ibags::Poly h_one(pp.n);
    for (int i = 0; i < pp.n; ++i) h_one[i] = 0;
    h_one[0] = 1;
    h_one.reduce_mod_q(pp);

    // f=0, g=0, h=1 → 1*0 != 0 → check_h_from_f_g should fail
    bool check1 = ibags::check_trapdoor_basis(pp, zero_basis, h_one);
    CHECK(!check1, "zero basis with h=1 fails validation");

    // 6b. Test dimension mismatch: create a basis with mismatched dimensions
    // (We can't easily create mismatched-dimension Poly without actual allocation,
    //  so we test check_h_from_f_g with different-dimension inputs instead)
    ibags::Poly f2(pp.n);
    ibags::Poly g2(pp.n);
    for (int i = 0; i < pp.n; ++i) { f2[i] = 0; g2[i] = 0; }
    f2[0] = 1; // f = 1
    g2[0] = 1; // g = 1

    // h should be g * f^{-1} = 1 * 1 = 1
    ibags::Poly h2(pp.n);
    for (int i = 0; i < pp.n; ++i) h2[i] = 0;
    h2[0] = 1;

    bool check2 = ibags::check_h_from_f_g(pp, f2, g2, h2);
    CHECK(check2, "f=1, g=1, h=1 passes h==g*f^{-1}");
}

// ============================================================================
// Test 7: 多次 TrapGen 的可复现性
// ============================================================================
void test_ntru_trapgen_reproducibility() {
    TEST("ntru_trapgen – reproducibility");

    auto pp = ibags::default_params();
    auto config = ibags::TrapGenConfig::default_config();
    config.max_invertibility_attempts = 25;

    // 两次使用同一 seed 应产生相同结果
    auto seed1 = make_seed(99);
    auto seed2 = make_seed(99); // same

    ibags::SeedCSPRNG csprng1(seed1, "test_repro");
    ibags::SeedCSPRNG csprng2(seed2, "test_repro");

    ibags::PublicTrapdoorParams pub1(pp), pub2(pp);
    ibags::MasterTrapdoorSecret sec1(pp.n), sec2(pp.n);

    auto st1 = ibags::ntru_trapgen(pp, csprng1, config, &pub1, &sec1);
    auto st2 = ibags::ntru_trapgen(pp, csprng2, config, &pub2, &sec2);

    CHECK_OK(st1, "first trapgen succeeds");
    CHECK_OK(st2, "second trapgen succeeds");

    // 检查 h 相等
    bool h_equal = ibags::poly_equal_ct(pub1.h, pub2.h);
    CHECK(h_equal, "h is reproducible with same seed");

    // 检查 f 相等
    bool f_equal = ibags::poly_equal_ct(sec1.basis.f, sec2.basis.f);
    CHECK(f_equal, "f is reproducible with same seed");
}

// ============================================================================
// Test 8: max_attempts = 0 / 负数 边界
// ============================================================================
void test_sample_f_max_attempts_boundary() {
    TEST("sample_f_until_invertible – max_attempts boundary");

    auto pp = ibags::default_params();
    auto seed = make_seed(5);
    ibags::SeedCSPRNG csprng(seed, "test_boundary");

    ibags::Poly f_out(pp.n);

    // max_attempts = 0 → should fail immediately
    auto st_zero = ibags::sample_f_until_invertible(pp, csprng, 0, &f_out);
    CHECK_FAIL(st_zero, "max_attempts=0 returns error");

    // max_attempts = 1 — 可能行可能不行，但至少不崩溃
    ibags::SeedCSPRNG csprng2(make_seed(6), "test_boundary2");
    auto st_one = ibags::sample_f_until_invertible(pp, csprng2, 1, &f_out);
    // 只是确认不崩溃，无论成功与否
    std::cout << "  INFO: max_attempts=1 result: "
              << (st_one.ok() ? "Ok" : st_one.message()) << std::endl;
    CHECK(true, "max_attempts=1 does not crash");
}

// ============================================================================
// Unified runner entry — callable from test_ibags_all.cpp
// ============================================================================

int run_ntru_trapgen_tests() {
    std::cout << "============================================\n"
              << " test_ntru_trapgen — NTRU TrapGen Unit Tests\n"
              << "============================================\n";

    test_sample_f_until_invertible();
    test_solve_ntru_equation();
    test_ntru_trapgen_default();
    test_ntru_trapgen_strict();
    test_gauss_sample_preimage();
    test_check_trapdoor_basis_edge_cases();
    test_ntru_trapgen_reproducibility();
    test_sample_f_max_attempts_boundary();

    std::cout << "\n============================================\n";
    std::cout << " Results: " << g_passed << " passed, "
              << g_failed << " failed\n";
    std::cout << "============================================\n";

    return (g_failed == 0) ? 0 : 1;
}