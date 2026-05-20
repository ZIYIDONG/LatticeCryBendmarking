/**
 * @file encode.cpp
 * @brief 规范序列化实现 — 多项式、身份标识、消息的编解码
 */

#include "../include_RLWEorNTRU/encode.h"

#include <cstring>
#include <algorithm>

namespace ibags {

// ============================================================================
// Polynomial Encode/Decode
// ============================================================================

Status poly_encode(const Poly& poly,
                   const Params& pp,
                   uint8_t* out,
                   size_t out_len)
{
    const size_t expected = static_cast<size_t>(pp.n) * COEFF_ENCODE_BYTES;
    if (out_len < expected) {
        return Status(ErrorCode::InvalidEncoding,
                      "poly_encode: output buffer too small");
    }

    if (!poly.is_canonical(pp)) {
        return Status(ErrorCode::InvalidPolynomial,
                      "poly_encode: polynomial not canonical");
    }

    for (int i = 0; i < pp.n; ++i) {
        int64_t coeff = poly[i];
        size_t base = static_cast<size_t>(i) * COEFF_ENCODE_BYTES;

        // 128-bit big-endian 编码 (系数通常 < 2^30，高位填零)
        for (int b = 0; b < 16; ++b) {
            int shift = (15 - b) * 8;
            out[base + b] = static_cast<uint8_t>(
                (static_cast<uint64_t>(coeff) >> shift) & 0xFF);
        }
    }

    return Status::Ok();
}

Status poly_decode(const uint8_t* encoded,
                   size_t len,
                   const Params& pp,
                   Poly* out)
{
    if (!out) {
        return Status(ErrorCode::InvalidEncoding,
                      "poly_decode: null output pointer");
    }

    const size_t expected = static_cast<size_t>(pp.n) * COEFF_ENCODE_BYTES;
    if (len < expected) {
        return Status(ErrorCode::InvalidEncoding,
                      "poly_decode: input too short");
    }

    std::vector<int64_t> coeffs(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        uint64_t val = 0;
        size_t base = static_cast<size_t>(i) * COEFF_ENCODE_BYTES;

        for (int b = 0; b < 16; ++b) {
            val = (val << 8) | encoded[base + b];
        }

        // 检查系数范围
        if (val >= static_cast<uint64_t>(pp.q)) {
            return Status(ErrorCode::InvalidEncoding,
                          "poly_decode: coefficient out of range");
        }

        coeffs[i] = static_cast<int64_t>(val);
    }

    *out = Poly(pp.n, std::move(coeffs));
    return Status::Ok();
}

// ============================================================================
// Identity Encode/Decode
// ============================================================================

Status encode_identity(const Identity& id,
                       uint8_t* out,
                       size_t* out_len)
{
    if (!out || !out_len) {
        return Status(ErrorCode::InvalidEncoding,
                      "encode_identity: null parameter");
    }

    size_t id_bytes = id.id.size();
    if (id_bytes > MAX_IDENTITY_BYTES) {
        return Status(ErrorCode::InvalidEncoding,
                      "encode_identity: identity too long");
    }

    // Format: domain_tag (1B) || len (1B) || id_bytes
    size_t total = 1 + 1 + id_bytes;

    out[0] = static_cast<uint8_t>(id.type);
    out[1] = static_cast<uint8_t>(id_bytes);
    if (id_bytes > 0) {
        std::memcpy(out + 2, id.id.data(), id_bytes);
    }

    *out_len = total;
    return Status::Ok();
}

Status decode_identity(const uint8_t* encoded,
                       size_t len,
                       Identity* out)
{
    if (!encoded || !out) {
        return Status(ErrorCode::InvalidEncoding,
                      "decode_identity: null parameter");
    }

    // Need at least: domain_tag (1B) + length (1B)
    if (len < 2) {
        return Status(ErrorCode::InvalidEncoding,
                      "decode_identity: input too short");
    }

    uint8_t type_byte = encoded[0];
    uint8_t id_len = encoded[1];

    if (len < 2 + static_cast<size_t>(id_len)) {
        return Status(ErrorCode::InvalidEncoding,
                      "decode_identity: truncated input");
    }

    // 验证类型
    switch (type_byte) {
        case static_cast<uint8_t>(IdentityType::VEHICLE_RID):
        case static_cast<uint8_t>(IdentityType::VEHICLE_PID):
        case static_cast<uint8_t>(IdentityType::RSU_RID):
        case static_cast<uint8_t>(IdentityType::RSU_PID):
            break;
        default:
            return Status(ErrorCode::InvalidEncoding,
                          "decode_identity: unknown identity type");
    }

    out->type = static_cast<IdentityType>(type_byte);
    out->id.assign(
        reinterpret_cast<const char*>(encoded + 2),
        id_len);

    return Status::Ok();
}

// ============================================================================
// Message Encode/Decode
// ============================================================================

Status encode_message(const Message& msg,
                      uint8_t* out,
                      size_t* out_len)
{
    if (!out || !out_len) {
        return Status(ErrorCode::InvalidEncoding,
                      "encode_message: null parameter");
    }

    size_t payload_len = msg.payload.size();
    if (payload_len > MAX_MESSAGE_PAYLOAD_BYTES) {
        return Status(ErrorCode::InvalidEncoding,
                      "encode_message: payload too long");
    }

    // Format: domain_tag (1B) || payload_len (2B, LE) || payload
    size_t total = 1 + 2 + payload_len;

    out[0] = static_cast<uint8_t>(msg.type);

    // payload_len: 2-byte little-endian
    out[1] = static_cast<uint8_t>(payload_len & 0xFF);
    out[2] = static_cast<uint8_t>((payload_len >> 8) & 0xFF);

    if (payload_len > 0) {
        std::memcpy(out + 3, msg.payload.data(), payload_len);
    }

    *out_len = total;
    return Status::Ok();
}

Status decode_message(const uint8_t* encoded,
                      size_t len,
                      Message* out)
{
    if (!encoded || !out) {
        return Status(ErrorCode::InvalidEncoding,
                      "decode_message: null parameter");
    }

    // Need at least: domain_tag (1B) + payload_len (2B)
    if (len < 3) {
        return Status(ErrorCode::InvalidEncoding,
                      "decode_message: input too short");
    }

    uint8_t type_byte = encoded[0];
    uint16_t payload_len = static_cast<uint16_t>(encoded[1])
                         | (static_cast<uint16_t>(encoded[2]) << 8);

    if (len < 3 + static_cast<size_t>(payload_len)) {
        return Status(ErrorCode::InvalidEncoding,
                      "decode_message: truncated input");
    }

    // 验证类型
    switch (type_byte) {
        case static_cast<uint8_t>(MessageType::V2V):
        case static_cast<uint8_t>(MessageType::V2I):
        case static_cast<uint8_t>(MessageType::RSU_BROADCAST):
            break;
        default:
            return Status(ErrorCode::InvalidEncoding,
                          "decode_message: unknown message type");
    }

    out->type = static_cast<MessageType>(type_byte);
    out->payload.assign(
        reinterpret_cast<const char*>(encoded + 3),
        payload_len);

    return Status::Ok();
}

} // namespace ibags