/**
 * @file xor_bytes.cpp
 * @brief Constant-time XOR operation for identity pseudonym derivation
 *
 * Used in: PID_V = RID_V XOR H1(a_j)  (Vehicle Registration)
 *          PID_R = RID_R XOR H1(b_j)  (RSU Registration)
 *
 * Design:
 *  - Constant-time loop: no early return on length mismatch leak pattern.
 *  - Returns Status for error propagation, consistent with the rest of IBAGS.
 *
 * Production notes:
 *  - This is a simple bytewise XOR; no security-critical timing issues.
 *  - The constant-time loop is a defense-in-depth measure against side-channel
 *    analysis of the RID/PID relationship.
 */

#include "../include_RLWEorNTRU/xor_bytes.h"
#include <cstring>

namespace ibags {

Status xor_bytes(ByteSpan a, ByteSpan b, uint8_t* out, size_t out_len) {
    // Validate inputs
    if (!out) {
        return {ErrorCode::InvalidEncoding, "xor_bytes: null output buffer"};
    }
    if (a.size != b.size) {
        return {ErrorCode::InvalidEncoding,
                "xor_bytes: length mismatch (a=" + std::to_string(a.size) +
                ", b=" + std::to_string(b.size) + ")"};
    }
    if (out_len < a.size) {
        return {ErrorCode::InvalidEncoding,
                "xor_bytes: output buffer too small (" + std::to_string(out_len) +
                " < " + std::to_string(a.size) + ")"};
    }

    if (a.size == 0) return Status::Ok();

    // Constant-time XOR: no early return, uniform loop
    for (size_t i = 0; i < a.size; ++i) {
        out[i] = a.data[i] ^ b.data[i];
    }

    return Status::Ok();
}

std::vector<uint8_t> xor_bytes_vec(ByteSpan a, ByteSpan b) {
    if (a.size != b.size) {
        throw std::invalid_argument(
            "xor_bytes_vec: length mismatch (" + std::to_string(a.size) +
            " vs " + std::to_string(b.size) + ")");
    }
    std::vector<uint8_t> result(a.size);
    for (size_t i = 0; i < a.size; ++i) {
        result[i] = a.data[i] ^ b.data[i];
    }
    return result;
}

} // namespace ibags
