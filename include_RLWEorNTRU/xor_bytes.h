#pragma once
/**
 * @file xor_bytes.h
 * @brief Constant-time byte XOR for PID = RID XOR mask derivation
 *
 * Used in identity pseudonym generation:
 *   PID_V = RID_V XOR H1_mask(a_j)
 *   PID_R = RID_R XOR H1_mask(b_j)
 *
 * See: "Secure Lattice-Based Aggregate Signature Scheme for VANETs" §3.3.1
 */

#include "xof.h"
#include "errors.h"
#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>

namespace ibags {

/**
 * @brief Constant-time XOR of two equal-length byte spans.
 *
 * Writes result to out[0..a.size-1].
 * Returns InvalidEncoding if lengths differ or buffer too small.
 *
 * @param a        First operand
 * @param b        Second operand
 * @param out      Output buffer
 * @param out_len  Output buffer length
 * @return         Status::Ok() or error
 */
[[nodiscard]] Status xor_bytes(ByteSpan a, ByteSpan b,
                                uint8_t* out, size_t out_len) noexcept;

/**
 * @brief Convenience: XOR and return as vector.
 *
 * Throws std::invalid_argument on length mismatch (for test/one-shot use).
 */
[[nodiscard]] std::vector<uint8_t> xor_bytes_vec(ByteSpan a, ByteSpan b);

} // namespace ibags
