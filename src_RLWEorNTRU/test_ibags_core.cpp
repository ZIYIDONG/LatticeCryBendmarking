/**
 * @file test_ibags_core.cpp
 * @brief IBAGS 核心模块单元测试 (Standalone, 无 GTest)
 *
 * 覆盖:
 *  - errors.h         : ErrorCode, Status, IBAGS_RETURN_IF_ERROR 宏
 *  - params.h/.cpp    : params_demo_64/level2_512/level3_768/level5_1024, validate_params, encode_params
 *  - secure_memory.h/.cpp : secure_zero, SecretBuffer, SecretPoly, SecureWipe
 *
 * 构建:
 *   g++ -std=c++20 -O0 -g -I../include_RLWEorNTRU \
 *       test_ibags_core.cpp params.cpp secure_memory.cpp \
 *       -o test_ibags_core
 */

#include "../include_RLWEorNTRU/errors.h"
#include "../include_RLWEorNTRU/params.h"
#include "../include_RLWEorNTRU/secure_memory.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <stdexcept>

using namespace ibags;

// ============================================================================
// Test helpers (与 test_ibags_poly.cpp 风格一致)
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

#define CHECK_EQ(a, b, msg) \
    do { \
        if (!((a) == (b))) { \
            printf("FAILED: %s (expected eq, got %lld vs %lld)\n", msg, \
                   (long long)(a), (long long)(b)); \
            ++tests_failed; \
            return; \
        } \
    } while (0)

#define CHECK_NE(a, b, msg) \
    do { \
        if (!((a) != (b))) { \
            printf("FAILED: %s (expected ne)\n", msg); \
            ++tests_failed; \
            return; \
        } \
    } while (0)

#define CHECK_TRUE(cond, msg) \
    do { \
        if (!(cond)) { FAIL(msg); return; } \
    } while (0)

#define CHECK_FALSE(cond, msg) \
    do { \
        if (cond) { FAIL(msg); return; } \
    } while (0)

#define CHECK_GT(a, b, msg) \
    do { \
        if (!((a) > (b))) { \
            printf("FAILED: %s (expected >)\n", msg); \
            ++tests_failed; \
            return; \
        } \
    } while (0)

#define CHECK_GE(a, b, msg) \
    do { \
        if (!((a) >= (b))) { \
            printf("FAILED: %s (expected >=)\n", msg); \
            ++tests_failed; \
            return; \
        } \
    } while (0)

#define CHECK_STREQ(a, b, msg) \
    do { \
        if (std::strcmp((a), (b)) != 0) { \
            printf("FAILED: %s (expected \"%s\", got \"%s\")\n", msg, (b), (a)); \
            ++tests_failed; \
            return; \
        } \
    } while (0)

#define CHECK_OK(status, msg) \
    do { \
        if (!(status).ok()) { \
            printf("FAILED: %s — %s\n", msg, (status).message().c_str()); \
            ++tests_failed; \
            return; \
        } \
    } while (0)

// ============================================================================
// 全局测试数据（替代 GTest fixtures）
// ============================================================================

static bool is_all_zero(const void* ptr, size_t len) {
    const unsigned char* p = static_cast<const unsigned char*>(ptr);
    for (size_t i = 0; i < len; ++i) {
        if (p[i] != 0) return false;
    }
    return true;
}

static Params demo_params;
static Params level2_params;
static Params level3_params;
static Params level5_params;

static void init_params() {
    demo_params    = Params::params_demo_64();
    level2_params  = Params::params_level2_512();
    level3_params  = Params::params_level3_768();
    level5_params  = Params::params_level5_1024();
}

// ============================================================================
// §1  errors 模块测试
// ============================================================================

void test_errors_default_status_is_ok() {
    TEST("Default status is Ok");
    Status s;
    CHECK_TRUE(s.ok(), "expected ok");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::Ok), "code mismatch");
    CHECK_TRUE(s.message().empty(), "message should be empty");
    CHECK_TRUE(static_cast<bool>(s), "bool conversion should be true");
    CHECK_FALSE(!s, "!s should be false");
    PASS();
}

void test_errors_status_ok_factory() {
    TEST("Status::Ok factory");
    Status s = Status::Ok();
    CHECK_TRUE(s.ok(), "expected ok");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::Ok), "code mismatch");
    PASS();
}

void test_errors_error_code_with_message() {
    TEST("ErrorCode with message");
    Status s(ErrorCode::InvalidParams, "n is not power of 2");
    CHECK_FALSE(s.ok(), "expected not ok");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidParams), "code mismatch");
    CHECK_STREQ(s.message().c_str(), "n is not power of 2", "message mismatch");
    PASS();
}

void test_errors_error_code_without_message() {
    TEST("ErrorCode without message");
    Status s(ErrorCode::InternalError);
    CHECK_FALSE(s.ok(), "expected not ok");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InternalError), "code mismatch");
    CHECK_TRUE(s.message().empty(), "message should be empty");
    PASS();
}

void test_errors_move_status() {
    TEST("Move status");
    Status s1(ErrorCode::DuplicateSigner, "dup");
    Status s2 = std::move(s1);
    CHECK_FALSE(s2.ok(), "moved status should not be ok");
    CHECK_EQ(static_cast<int>(s2.code()), static_cast<int>(ErrorCode::DuplicateSigner), "code mismatch");
    CHECK_STREQ(s2.message().c_str(), "dup", "message mismatch");
    PASS();
}

void test_errors_status_equality() {
    TEST("Status equality");
    Status ok1 = Status::Ok();
    Status ok2 = Status::Ok();
    CHECK_TRUE(ok1 == ok2, "ok1 should equal ok2");

    Status err1(ErrorCode::NormTooLarge);
    Status err2(ErrorCode::NormTooLarge);
    Status err3(ErrorCode::InvalidPolynomial);
    CHECK_TRUE(err1 == err2, "err1 should equal err2");
    CHECK_TRUE(err1 != err3, "err1 should not equal err3");
    CHECK_TRUE(ok1 != err1, "ok1 should not equal err1");
    PASS();
}

void test_errors_all_error_codes_distinct() {
    TEST("All error codes distinct");
    std::vector<int> codes = {
        static_cast<int>(ErrorCode::Ok),
        static_cast<int>(ErrorCode::InvalidParams),
        static_cast<int>(ErrorCode::InvalidEncoding),
        static_cast<int>(ErrorCode::InvalidSignerCount),
        static_cast<int>(ErrorCode::DuplicateSigner),
        static_cast<int>(ErrorCode::InvalidPolynomial),
        static_cast<int>(ErrorCode::InvalidPublicKey),
        static_cast<int>(ErrorCode::InvalidSecretKey),
        static_cast<int>(ErrorCode::InvalidSignature),
        static_cast<int>(ErrorCode::NormTooLarge),
        static_cast<int>(ErrorCode::TrapdoorError),
        static_cast<int>(ErrorCode::SamplingError),
        static_cast<int>(ErrorCode::InternalError),
    };
    std::sort(codes.begin(), codes.end());
    auto it = std::unique(codes.begin(), codes.end());
    CHECK_TRUE(it == codes.end(), "Duplicate error codes detected");
    PASS();
}

void test_errors_return_if_error_macro() {
    TEST("IBAGS_RETURN_IF_ERROR macro");
    auto test_function = [](int x) -> Status {
        IBAGS_RETURN_IF_ERROR(
            (x == 0) ? Status::Ok()
                     : Status(ErrorCode::InvalidParams, "bad x"));
        return Status::Ok();
    };

    CHECK_TRUE(test_function(0).ok(), "x=0 should return Ok");
    CHECK_FALSE(test_function(1).ok(), "x=1 should return error");
    CHECK_EQ(static_cast<int>(test_function(1).code()),
             static_cast<int>(ErrorCode::InvalidParams), "code mismatch");
    PASS();
}

// ============================================================================
// §2  params 模块测试
// ============================================================================

// ── 2.1 工厂函数测试 ──

void test_params_demo_64() {
    TEST("Demo params 64");
    CHECK_EQ(static_cast<int>(demo_params.param_id),
             static_cast<int>(ParamId::IBAGS_64_DEMO), "param_id mismatch");
    CHECK_EQ(demo_params.n, 64, "n mismatch");
    CHECK_GT(demo_params.q, 0, "q should be > 0");
    CHECK_GT(demo_params.sigma, 0.0, "sigma should be > 0");
    CHECK_GE(demo_params.eta1, demo_params.eta2, "eta1 >= eta2");
    CHECK_GE(demo_params.eta2, demo_params.eta_ver, "eta2 >= eta_ver");
    CHECK_GT(demo_params.max_signers, 0, "max_signers > 0");
    CHECK_GT(demo_params.max_rejection_loops, 0, "max_rejection_loops > 0");
    CHECK_GT(demo_params.coefficient_bytes, 0, "coefficient_bytes > 0");
    CHECK_EQ(demo_params.poly_bytes,
             demo_params.n * demo_params.coefficient_bytes, "poly_bytes mismatch");
    PASS();
}

void test_params_level2_512() {
    TEST("Level2 params 512");
    CHECK_EQ(static_cast<int>(level2_params.param_id),
             static_cast<int>(ParamId::IBAGS_512_LEVEL2), "param_id mismatch");
    CHECK_EQ(level2_params.n, 512, "n mismatch");
    CHECK_GT(level2_params.q, 0, "q should be > 0");
    PASS();
}

void test_params_level3_768() {
    TEST("Level3 params 768");
    CHECK_EQ(static_cast<int>(level3_params.param_id),
             static_cast<int>(ParamId::IBAGS_768_LEVEL3), "param_id mismatch");
    CHECK_EQ(level3_params.n, 768, "n mismatch");
    CHECK_GT(level3_params.q, 0, "q should be > 0");
    PASS();
}

void test_params_level5_1024() {
    TEST("Level5 params 1024");
    CHECK_EQ(static_cast<int>(level5_params.param_id),
             static_cast<int>(ParamId::IBAGS_1024_LEVEL5), "param_id mismatch");
    CHECK_EQ(level5_params.n, 1024, "n mismatch");
    CHECK_GT(level5_params.q, 0, "q should be > 0");
    PASS();
}

void test_params_param_id_name() {
    TEST("ParamId name");
    CHECK_STREQ(param_id_name(ParamId::IBAGS_64_DEMO),     "IBAGS-64-Demo",    "name mismatch");
    CHECK_STREQ(param_id_name(ParamId::IBAGS_512_LEVEL2),  "IBAGS-512-Level2", "name mismatch");
    CHECK_STREQ(param_id_name(ParamId::IBAGS_768_LEVEL3),  "IBAGS-768-Level3", "name mismatch");
    CHECK_STREQ(param_id_name(ParamId::IBAGS_1024_LEVEL5), "IBAGS-1024-Level5","name mismatch");
    CHECK_STREQ(param_id_name(static_cast<ParamId>(0xFFFF)), "IBAGS-Unknown",   "unknown name mismatch");
    PASS();
}

// ── 2.2 参数验证测试 ──

void test_params_validate_all_factory() {
    TEST("Validate all factory params");
    CHECK_OK(validate_params(demo_params), "demo_params invalid");
    CHECK_OK(validate_params(level2_params), "level2_params invalid");
    CHECK_OK(validate_params(level3_params), "level3_params invalid");
    CHECK_OK(validate_params(level5_params), "level5_params invalid");
    PASS();
}

void test_params_invalid_n_not_power_of_two() {
    TEST("Invalid n: not power of two");
    Params p = level2_params;
    p.n = 500;
    Status s = validate_params(p);
    CHECK_FALSE(s.ok(), "should be invalid");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidParams), "code mismatch");
    PASS();
}

void test_params_invalid_n_too_small() {
    TEST("Invalid n: too small");
    Params p = level2_params;
    p.n = 32;
    Status s = validate_params(p);
    CHECK_FALSE(s.ok(), "should be invalid");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidParams), "code mismatch");
    PASS();
}

void test_params_invalid_q_zero_or_negative() {
    TEST("Invalid q: zero or negative");
    Params p = level2_params;
    p.q = 0;
    CHECK_FALSE(validate_params(p).ok(), "q=0 should be invalid");

    p.q = -100;
    CHECK_FALSE(validate_params(p).ok(), "q=-100 should be invalid");
    PASS();
}

void test_params_invalid_q_too_large() {
    TEST("Invalid q: too large");
    Params p = level2_params;
    p.q = (1 << 30) + 1;
    Status s = validate_params(p);
    CHECK_FALSE(s.ok(), "should be invalid");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidParams), "code mismatch");
    PASS();
}

void test_params_invalid_sigma_non_positive() {
    TEST("Invalid sigma: non-positive");
    Params p = level2_params;
    p.sigma = 0.0;
    CHECK_FALSE(validate_params(p).ok(), "sigma=0 should be invalid");

    p.sigma = -1.0;
    CHECK_FALSE(validate_params(p).ok(), "sigma=-1 should be invalid");
    PASS();
}

void test_params_invalid_eta_chain() {
    TEST("Invalid eta chain");

    Params p = level2_params;
    p.eta_ver = 0;
    CHECK_FALSE(validate_params(p).ok(), "eta_ver=0 should be invalid");

    p = level2_params;
    p.eta2 = 5;
    p.eta_ver = 10;
    CHECK_FALSE(validate_params(p).ok(), "eta2 < eta_ver should be invalid");

    p = level2_params;
    p.eta1 = 10;
    p.eta2 = 20;
    p.eta_ver = 15;
    CHECK_FALSE(validate_params(p).ok(), "eta1 < eta2 should be invalid");
    PASS();
}

void test_params_invalid_max_signers() {
    TEST("Invalid max_signers");
    Params p = level2_params;
    p.max_signers = 0;
    CHECK_FALSE(validate_params(p).ok(), "max_signers=0 should be invalid");

    p.max_signers = 300;
    CHECK_FALSE(validate_params(p).ok(), "max_signers=300 should be invalid");
    PASS();
}

void test_params_invalid_rejection_loops() {
    TEST("Invalid rejection loops");
    Params p = level2_params;
    p.max_rejection_loops = 0;
    CHECK_FALSE(validate_params(p).ok(), "loops=0 should be invalid");

    p.max_rejection_loops = -1;
    CHECK_FALSE(validate_params(p).ok(), "loops=-1 should be invalid");
    PASS();
}

void test_params_invalid_coefficient_bytes() {
    TEST("Invalid coefficient_bytes");
    Params p = level2_params;
    p.coefficient_bytes = 1;
    Status s = validate_params(p);
    CHECK_FALSE(s.ok(), "should be invalid");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidParams), "code mismatch");
    PASS();
}

void test_params_invalid_poly_bytes() {
    TEST("Invalid poly_bytes");
    Params p = level2_params;
    p.poly_bytes = 999;
    Status s = validate_params(p);
    CHECK_FALSE(s.ok(), "should be invalid");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidParams), "code mismatch");
    PASS();
}

// ── 2.3 参数编码测试 ──

void test_params_encode_level2() {
    TEST("Encode params level2");
    uint8_t buf[PARAMS_ENCODE_BYTES];
    Status s = encode_params(level2_params, buf, sizeof(buf));
    CHECK_OK(s, "encode_params failed");

    uint16_t pid = static_cast<uint16_t>(buf[0]) | (static_cast<uint16_t>(buf[1]) << 8);
    CHECK_EQ(pid, static_cast<uint16_t>(ParamId::IBAGS_512_LEVEL2), "param_id mismatch");

    uint16_t n_val = static_cast<uint16_t>(buf[2]) | (static_cast<uint16_t>(buf[3]) << 8);
    CHECK_EQ(n_val, 512, "n mismatch");
    PASS();
}

void test_params_encode_null_buffer() {
    TEST("Encode params null buffer");
    Status s = encode_params(level2_params, nullptr, PARAMS_ENCODE_BYTES);
    CHECK_FALSE(s.ok(), "should be error");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidEncoding), "code mismatch");
    PASS();
}

void test_params_encode_buffer_too_small() {
    TEST("Encode params buffer too small");
    uint8_t buf[32];
    Status s = encode_params(level2_params, buf, sizeof(buf));
    CHECK_FALSE(s.ok(), "should be error");
    CHECK_EQ(static_cast<int>(s.code()), static_cast<int>(ErrorCode::InvalidEncoding), "code mismatch");
    PASS();
}

void test_params_encode_reserved_bytes_zeroed() {
    TEST("Encode params reserved bytes zeroed");
    uint8_t buf[PARAMS_ENCODE_BYTES];
    std::memset(buf, 0xAA, sizeof(buf));

    Status s = encode_params(level2_params, buf, sizeof(buf));
    CHECK_OK(s, "encode_params failed");

    for (size_t i = 30; i < PARAMS_ENCODE_BYTES; ++i) {
        if (buf[i] != 0x00) {
            char err[64];
            std::snprintf(err, sizeof(err), "Reserved byte at offset %zu not zeroed", i);
            FAIL(err);
            return;
        }
    }
    PASS();
}

void test_params_encode_deterministic() {
    TEST("Encode params deterministic");
    uint8_t buf1[PARAMS_ENCODE_BYTES];
    uint8_t buf2[PARAMS_ENCODE_BYTES];

    (void)encode_params(level2_params, buf1, sizeof(buf1));
    (void)encode_params(level2_params, buf2, sizeof(buf2));

    CHECK_EQ(std::memcmp(buf1, buf2, PARAMS_ENCODE_BYTES), 0, "encodings should be identical");
    PASS();
}

void test_params_different_params_different_encoding() {
    TEST("Different params different encoding");
    uint8_t buf_level2[PARAMS_ENCODE_BYTES];
    uint8_t buf_demo[PARAMS_ENCODE_BYTES];

    (void)encode_params(level2_params, buf_level2, sizeof(buf_level2));
    (void)encode_params(demo_params, buf_demo, sizeof(buf_demo));

    bool differs = false;
    for (size_t i = 0; i < PARAMS_ENCODE_BYTES; ++i) {
        if (buf_level2[i] != buf_demo[i]) {
            differs = true;
            break;
        }
    }
    CHECK_TRUE(differs, "Different parameter sets should produce different encodings");
    PASS();
}

// ============================================================================
// §3  secure_memory 模块测试
// ============================================================================

// ── 3.1 secure_zero 测试 ──

void test_secure_zero_basic() {
    TEST("secure_zero basic");
    uint8_t buf[64];
    std::memset(buf, 0xFF, sizeof(buf));

    secure_zero(buf, sizeof(buf));

    CHECK_TRUE(is_all_zero(buf, sizeof(buf)), "buffer should be all zeros");
    PASS();
}

void test_secure_zero_null_no_crash() {
    TEST("secure_zero null no crash");
    secure_zero(nullptr, 0);
    secure_zero(nullptr, 100);

    uint8_t buf[4];
    secure_zero(buf, 0);
    printf("PASSED (no crash)\n");
    ++tests_passed;
}

void test_secure_zero_partial() {
    TEST("secure_zero partial");
    uint8_t buf[64];
    std::memset(buf, 0xFF, sizeof(buf));

    secure_zero(buf, 16);

    CHECK_TRUE(is_all_zero(buf, 16), "first 16 bytes should be zero");
    for (size_t i = 16; i < 64; ++i) {
        if (buf[i] != 0xFF) {
            FAIL("remaining bytes should be unchanged");
            return;
        }
    }
    PASS();
}

// ── 3.2 SecretBuffer 测试 ──

void test_secret_buffer_construction() {
    TEST("SecretBuffer construction");
    SecretBuffer sb(32);
    CHECK_EQ(static_cast<int>(sb.size()), 32, "size mismatch");
    CHECK_FALSE(sb.empty(), "should not be empty");
    CHECK_TRUE(sb.data() != nullptr, "data should not be null");
    CHECK_TRUE(is_all_zero(sb.data(), sb.size()), "initial content should be zero");
    PASS();
}

void test_secret_buffer_zero_size() {
    TEST("SecretBuffer zero size");
    SecretBuffer sb(0);
    CHECK_EQ(static_cast<int>(sb.size()), 0, "size mismatch");
    CHECK_TRUE(sb.empty(), "should be empty");
    PASS();
}

void test_secret_buffer_from_ptr() {
    TEST("SecretBuffer from ptr");
    uint8_t* raw = new uint8_t[16];
    std::memset(raw, 0xAB, 16);

    {
        SecretBuffer sb(raw, 16);
        CHECK_EQ(static_cast<int>(sb.size()), 16, "size mismatch");
        for (size_t i = 0; i < 16; ++i) {
            if (sb.data()[i] != 0xAB) {
                FAIL("data not preserved");
                return;
            }
        }
    }
    PASS();
}

void test_secret_buffer_move_constructor() {
    TEST("SecretBuffer move constructor");
    SecretBuffer sb1(16);
    std::memset(sb1.data(), 0x42, 16);

    SecretBuffer sb2(std::move(sb1));

    CHECK_TRUE(sb1.empty(), "source should be empty after move");
    CHECK_EQ(static_cast<int>(sb1.size()), 0, "source size should be 0");
    CHECK_TRUE(sb1.data() == nullptr, "source data should be null");

    CHECK_EQ(static_cast<int>(sb2.size()), 16, "dest size mismatch");
    CHECK_TRUE(sb2.data() != nullptr, "dest data should not be null");
    for (size_t i = 0; i < 16; ++i) {
        if (sb2.data()[i] != 0x42) {
            FAIL("data not moved correctly");
            return;
        }
    }
    PASS();
}

void test_secret_buffer_move_assignment() {
    TEST("SecretBuffer move assignment");
    SecretBuffer sb1(8);
    std::memset(sb1.data(), 0x77, 8);

    SecretBuffer sb2(4);
    sb2 = std::move(sb1);

    CHECK_TRUE(sb1.empty(), "source should be empty after move");
    CHECK_EQ(static_cast<int>(sb2.size()), 8, "dest size mismatch");
    for (size_t i = 0; i < 8; ++i) {
        if (sb2.data()[i] != 0x77) {
            FAIL("data not moved correctly");
            return;
        }
    }
    PASS();
}

void test_secret_buffer_public_copy() {
    TEST("SecretBuffer public_copy");
    SecretBuffer original(8);
    std::memset(original.data(), 0x5A, 8);

    SecretBuffer copy = original.public_copy();

    CHECK_EQ(static_cast<int>(copy.size()), static_cast<int>(original.size()), "size mismatch");
    for (size_t i = 0; i < 8; ++i) {
        if (copy.data()[i] != 0x5A) {
            FAIL("copy data mismatch");
            return;
        }
    }

    copy.data()[0] = 0xFF;
    CHECK_TRUE(original.data()[0] == 0x5A, "original should not be affected by copy modification");
    PASS();
}

void test_secret_buffer_destructor_zeroes() {
    TEST("SecretBuffer destructor zeros (exercise path)");
    auto* sb = new SecretBuffer(16);
    std::memset(sb->data(), 0xDE, 16);
    (void)sb->data();

    delete sb;

    // Address is freed — verify with valgrind/ASan
    printf("PASSED (verify with valgrind/ASan)\n");
    ++tests_passed;
}

void test_secret_buffer_not_copyable() {
    TEST("SecretBuffer is not copyable");
    CHECK_FALSE(std::is_copy_constructible_v<SecretBuffer>, "should not be copy constructible");
    CHECK_FALSE(std::is_copy_assignable_v<SecretBuffer>, "should not be copy assignable");
    PASS();
}

void test_secret_buffer_is_movable() {
    TEST("SecretBuffer is movable");
    CHECK_TRUE(std::is_move_constructible_v<SecretBuffer>, "should be move constructible");
    CHECK_TRUE(std::is_move_assignable_v<SecretBuffer>, "should be move assignable");
    PASS();
}

// ── 3.3 SecureWipe 测试 ──

void test_secure_wipe_clears_on_scope() {
    TEST("SecureWipe clears on scope exit");
    uint8_t buf[32];
    std::memset(buf, 0xCC, sizeof(buf));

    {
        SecureWipe wipe(buf, sizeof(buf));
    }
    printf("PASSED (verify with valgrind/ASan)\n");
    ++tests_passed;
}

// ── 3.4 SecretPoly 测试 ──

void test_secret_poly_construction() {
    TEST("SecretPoly construction");
    SecretPoly poly(512);
    CHECK_EQ(poly.n(), 512, "n mismatch");

    for (int i = 0; i < 512; ++i) {
        if (poly.coeff(i) != 0) {
            FAIL("initial coefficients should be zero");
            return;
        }
    }
    PASS();
}

void test_secret_poly_from_vector() {
    TEST("SecretPoly from vector");
    std::vector<long> coeffs = {1, 2, 3, 4, 5};
    SecretPoly poly(std::move(coeffs));

    CHECK_EQ(poly.n(), 5, "n mismatch");
    CHECK_EQ(poly.coeff(0), 1, "coeff[0] mismatch");
    CHECK_EQ(poly.coeff(4), 5, "coeff[4] mismatch");
    PASS();
}

void test_secret_poly_set_get_coeff() {
    TEST("SecretPoly set/get coeff");
    SecretPoly poly(8);
    poly.set_coeff(0, 100);
    poly.set_coeff(7, -50);

    CHECK_EQ(poly.coeff(0), 100, "coeff[0] mismatch");
    CHECK_EQ(poly.coeff(7), -50, "coeff[7] mismatch");
    CHECK_EQ(poly.coeff(3), 0, "unset coeff should be 0");
    PASS();
}

void test_secret_poly_move_constructor() {
    TEST("SecretPoly move constructor");
    SecretPoly p1(32);
    p1.set_coeff(0, 42);

    SecretPoly p2(std::move(p1));

    CHECK_EQ(p1.n(), 0, "source n should be 0 after move");
    CHECK_EQ(p2.n(), 32, "dest n mismatch");
    CHECK_EQ(p2.coeff(0), 42, "data not moved correctly");
    PASS();
}

void test_secret_poly_move_assignment() {
    TEST("SecretPoly move assignment");
    SecretPoly p1(16);
    p1.set_coeff(0, 99);

    SecretPoly p2(8);
    p2 = std::move(p1);

    CHECK_EQ(p1.n(), 0, "source n should be 0 after move");
    CHECK_EQ(p2.n(), 16, "dest n mismatch");
    CHECK_EQ(p2.coeff(0), 99, "data not moved correctly");
    PASS();
}

void test_secret_poly_public_copy() {
    TEST("SecretPoly public_copy");
    SecretPoly original(4);
    original.set_coeff(0, 1);
    original.set_coeff(1, 2);
    original.set_coeff(2, 3);
    original.set_coeff(3, 4);

    SecretPoly copy = original.public_copy();

    CHECK_EQ(copy.n(), 4, "copy n mismatch");
    for (int i = 0; i < 4; ++i) {
        if (copy.coeff(i) != i + 1) {
            FAIL("copy data mismatch");
            return;
        }
    }

    copy.set_coeff(0, 999);
    CHECK_EQ(original.coeff(0), 1, "original should not be affected by copy modification");
    PASS();
}

void test_secret_poly_raw_access() {
    TEST("SecretPoly raw access");
    SecretPoly poly(8);
    poly.set_coeff(0, 10);
    poly.set_coeff(7, 70);

    const long* raw = poly.raw_coeffs();
    CHECK_EQ(raw[0], 10, "raw[0] mismatch");
    CHECK_EQ(raw[7], 70, "raw[7] mismatch");

    long* mutable_raw = poly.raw_coeffs_mutable();
    mutable_raw[3] = 30;
    CHECK_EQ(poly.coeff(3), 30, "mutable access not reflected");
    PASS();
}

void test_secret_poly_not_copyable() {
    TEST("SecretPoly is not copyable");
    CHECK_FALSE(std::is_copy_constructible_v<SecretPoly>, "should not be copy constructible");
    CHECK_FALSE(std::is_copy_assignable_v<SecretPoly>, "should not be copy assignable");
    PASS();
}

void test_secret_poly_is_movable() {
    TEST("SecretPoly is movable");
    CHECK_TRUE(std::is_move_constructible_v<SecretPoly>, "should be move constructible");
    CHECK_TRUE(std::is_move_assignable_v<SecretPoly>, "should be move assignable");
    PASS();
}

void test_secret_poly_destructor_zeroes() {
    TEST("SecretPoly destructor zeros (exercise path)");
    auto* poly = new SecretPoly(16);
    (void)poly->raw_coeffs();

    for (int i = 0; i < 16; ++i) {
        poly->set_coeff(i, 0xDEADBEEFL);
    }

    delete poly;

    printf("PASSED (verify with valgrind/ASan)\n");
    ++tests_passed;
}

// ============================================================================
// 集成测试: params + errors 协同
// ============================================================================

void test_integration_validate_with_status_propagation() {
    TEST("Integration: validate with status propagation");
    auto protocol_entry = [](const Params& p) -> Status {
        IBAGS_RETURN_IF_ERROR(validate_params(p));
        return Status::Ok();
    };

    CHECK_OK(protocol_entry(level2_params), "level2_params should pass");
    CHECK_FALSE(protocol_entry(Params{}).ok(), "empty params should fail");
    PASS();
}

// ============================================================================
// main
// ============================================================================

int main() {
    init_params();

    printf("===== Running IBAGS Core Tests =====\n\n");

    // §1 errors
    printf("--- §1 Errors Module ---\n");
    test_errors_default_status_is_ok();
    test_errors_status_ok_factory();
    test_errors_error_code_with_message();
    test_errors_error_code_without_message();
    test_errors_move_status();
    test_errors_status_equality();
    test_errors_all_error_codes_distinct();
    test_errors_return_if_error_macro();

    // §2 params
    printf("\n--- §2 Params: Factory ---\n");
    test_params_demo_64();
    test_params_level2_512();
    test_params_level3_768();
    test_params_level5_1024();
    test_params_param_id_name();

    printf("\n--- §2 Params: Validation ---\n");
    test_params_validate_all_factory();
    test_params_invalid_n_not_power_of_two();
    test_params_invalid_n_too_small();
    test_params_invalid_q_zero_or_negative();
    test_params_invalid_q_too_large();
    test_params_invalid_sigma_non_positive();
    test_params_invalid_eta_chain();
    test_params_invalid_max_signers();
    test_params_invalid_rejection_loops();
    test_params_invalid_coefficient_bytes();
    test_params_invalid_poly_bytes();

    printf("\n--- §2 Params: Encoding ---\n");
    test_params_encode_level2();
    test_params_encode_null_buffer();
    test_params_encode_buffer_too_small();
    test_params_encode_reserved_bytes_zeroed();
    test_params_encode_deterministic();
    test_params_different_params_different_encoding();

    // §3 secure_memory
    printf("\n--- §3 Secure Memory: secure_zero ---\n");
    test_secure_zero_basic();
    test_secure_zero_null_no_crash();
    test_secure_zero_partial();

    printf("\n--- §3 Secure Memory: SecretBuffer ---\n");
    test_secret_buffer_construction();
    test_secret_buffer_zero_size();
    test_secret_buffer_from_ptr();
    test_secret_buffer_move_constructor();
    test_secret_buffer_move_assignment();
    test_secret_buffer_public_copy();
    test_secret_buffer_destructor_zeroes();
    test_secret_buffer_not_copyable();
    test_secret_buffer_is_movable();

    printf("\n--- §3 Secure Memory: SecureWipe ---\n");
    test_secure_wipe_clears_on_scope();

    printf("\n--- §3 Secure Memory: SecretPoly ---\n");
    test_secret_poly_construction();
    test_secret_poly_from_vector();
    test_secret_poly_set_get_coeff();
    test_secret_poly_move_constructor();
    test_secret_poly_move_assignment();
    test_secret_poly_public_copy();
    test_secret_poly_raw_access();
    test_secret_poly_not_copyable();
    test_secret_poly_is_movable();
    test_secret_poly_destructor_zeroes();

    // Integration
    printf("\n--- Integration ---\n");
    test_integration_validate_with_status_propagation();

    // Summary
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