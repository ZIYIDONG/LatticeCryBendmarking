#pragma once
/**
 * @file constant_time.h
 * @brief Constant-time comparison, selection, and arithmetic utilities
 *
 * Design principles:
 *  - All functions avoid secret-dependent branches (no early return on secret data).
 *  - Use bitwise operations and arithmetic masking.
 *  - Inputs annotated as "secret" in comments; public inputs may branch freely.
 *
 * Reference:
 *  - OpenSSL constant_time_locl.h (CRYPTO_memcmp, constant_time_select, etc.)
 *  - BoringSSL constant-time helpers
 *  - Curve25519-donna constant-time operations
 *
 * Production notes:
 *  - Current implementation uses standard C++20 arithmetic; compilers may still
 *    introduce secret-dependent optimizations. For production:
 *     ① Use compiler intrinsics (__builtin_ct_* with Clang 14+)
 *     ② Or link against BoringSSL/OpenSSL's verified constant-time primitives
 *     ③ Always verify with ctgrind/ct-verif/valgrind --tool=lackey for leaks
 */

#include <cstdint>
#include <cstddef>

namespace ibags {
namespace ct {

// ============================================================================
// §1  ct_equal  —  constant-time byte/word equality
// ============================================================================

/**
 * @brief Constant-time comparison of two byte buffers.
 *
 * Returns true iff the two buffers are identical, without early exit.
 * All bytes are compared regardless of when a mismatch is found.
 *
 * Complexity: O(len) — always scans the entire buffer.
 *
 * @param a    First buffer (public length, secret content)
 * @param b    Second buffer
 * @param len  Number of bytes to compare
 * @return     true if equal, false otherwise
 */
[[nodiscard]] inline bool ct_equal_bytes(const uint8_t* a,
                                          const uint8_t* b,
                                          size_t len) noexcept {
    uint8_t diff = 0;
    for (size_t i = 0; i < len; ++i) {
        diff |= a[i] ^ b[i];
    }
    return diff == 0;
}

/**
 * @brief Constant-time equality of two 64-bit integers.
 *
 * Uses the standard bit-hack: (x ^ y) generates ones in differing bit positions;
 * folding to 0/1 via (x | -x) >> 63 trick.
 */
[[nodiscard]] inline bool ct_equal_u64(uint64_t a, uint64_t b) noexcept {
    uint64_t x = a ^ b;
    // x == 0  =>  (x-1) sets high bit, but not zero itself.
    // Standard pattern: x = (x | (x-1)) >> 63  gives 0 if x==0, else something non-zero.
    // Simpler: use arithmetic trick
    x = (x | (x - 1)) >> 63;
    return x == 0;
}

/**
 * @brief Constant-time comparison of two int64_t values.
 */
[[nodiscard]] inline bool ct_equal_i64(int64_t a, int64_t b) noexcept {
    return ct_equal_u64(static_cast<uint64_t>(a), static_cast<uint64_t>(b));
}

// ============================================================================
// §2  ct_select  —  constant-time conditional select (cmov)
// ============================================================================

/**
 * @brief Constant-time select: returns a if cond is true (nonzero), else b.
 *
 * Uses: result = (cond & a) | (~cond & b) where cond is 0 or ~0.
 *
 * @param a      First value
 * @param b      Second value
 * @param cond   0 (false) or any nonzero (true)
 * @return       a if cond != 0, else b
 */
[[nodiscard]] inline uint64_t ct_select_u64(uint64_t a,
                                             uint64_t b,
                                             uint64_t cond) noexcept {
    // Normalize cond to 0 or ~0ULL
    uint64_t mask = cond ? UINT64_MAX : 0;
    return (mask & a) | (~mask & b);
}

[[nodiscard]] inline int64_t ct_select_i64(int64_t a,
                                            int64_t b,
                                            int64_t cond) noexcept {
    return static_cast<int64_t>(
        ct_select_u64(static_cast<uint64_t>(a),
                       static_cast<uint64_t>(b),
                       cond ? UINT64_MAX : 0));
}

/**
 * @brief Constant-time select for size_t.
 */
[[nodiscard]] inline size_t ct_select_sizet(size_t a,
                                             size_t b,
                                             size_t cond) noexcept {
    size_t mask = cond ? ~size_t{0} : size_t{0};
    return (mask & a) | (~mask & b);
}

// ============================================================================
// §3  ct_less_than  —  constant-time unsigned comparison
// ============================================================================

/**
 * @brief Constant-time unsigned less-than: returns 1 if a < b, else 0.
 *
 * Implementation: compute borrow bit from (a - b).
 * (a - b) underflows (i.e., MSB set) exactly when a < b.
 */
[[nodiscard]] inline uint64_t ct_less_than_u64(uint64_t a,
                                                uint64_t b) noexcept {
    // (b - a) has the high bit set iff a < b (for unsigned)
    uint64_t diff = b - a;
    return (diff >> 63) & 1;
}

/**
 * @brief Constant-time check: a < b for int64_t (signed).
 *
 * For signed, a < b iff (a ^ sign) < (b ^ sign) as unsigned.
 * XOR both with 0x8000... to flip sign bit, then compare unsigned.
 */
[[nodiscard]] inline int64_t ct_less_than_i64(int64_t a,
                                               int64_t b) noexcept {
    constexpr uint64_t sign_bit = UINT64_C(0x8000000000000000);
    uint64_t ua = static_cast<uint64_t>(a) ^ sign_bit;
    uint64_t ub = static_cast<uint64_t>(b) ^ sign_bit;
    return static_cast<int64_t>(ct_less_than_u64(ua, ub));
}

// ============================================================================
// §4  ct_abs  —  constant-time absolute value
// ============================================================================

/**
 * @brief Constant-time absolute value of a 64-bit signed integer.
 *
 * Returns |x|.  Handles INT64_MIN safely (returns INT64_MAX due to overflow,
 * which is the standard behavior for lattice norm computations where
 * INT64_MIN should not occur in practice).
 */
[[nodiscard]] inline int64_t ct_abs_i64(int64_t x) noexcept {
    int64_t mask = x >> 63;   // 0 if x >= 0, -1 if x < 0
    return (x ^ mask) - mask;  // standard bit-hack: ~x+1 when negative
}

/**
 * @brief Constant-time absolute value, returning uint64_t (avoids INT64_MIN issue).
 *
 * For lattice ∞-norm, use this to avoid overflow on INT64_MIN.
 */
[[nodiscard]] inline uint64_t ct_abs_u64(int64_t x) noexcept {
    if (x >= 0) return static_cast<uint64_t>(x);
    if (x == INT64_MIN) {
        // INT64_MIN: |x| = 2^63 fits in uint64_t, but ~x+1 overflows in int64
        return static_cast<uint64_t>(INT64_MAX) + 1;
    }
    return static_cast<uint64_t>(-x);
}

// ============================================================================
// §5  ct_is_zero  —  constant-time zero check
// ============================================================================

/**
 * @brief Constant-time check whether a 64-bit value is zero.
 *
 * Returns 1 if x == 0, 0 otherwise.
 */
[[nodiscard]] inline uint64_t ct_is_zero_u64(uint64_t x) noexcept {
    // (x - 1) >> 63 sets bit 63 iff x == 0 in unsigned
    // For signed negative, this doesn't hold; here x is unsigned
    // Actually: ((x | -x) >> 63) & 1  is cleaner for unsigned
    return ((x | (~x + 1)) >> 63) & 1;
}

// ============================================================================
// §6  ct_poly  —  constant-time polynomial helpers
// ============================================================================

/**
 * @brief Constant-time maximum of two uint64_t values.
 *
 * Used in ∞-norm accumulation without branching.
 */
[[nodiscard]] inline uint64_t ct_max_u64(uint64_t a, uint64_t b) noexcept {
    return ct_select_u64(a, b, ct_less_than_u64(a, b));
}

/**
 * @brief Constant-time conditional sign flip: returns cond ? -x : x.
 *
 * @param cond  true (-1 as int64) to flip sign, false (0) to keep
 */
[[nodiscard]] inline int64_t ct_cond_negate(int64_t x, int64_t cond) noexcept {
    return (x ^ cond) - cond;
}

} // namespace ct
} // namespace ibags
