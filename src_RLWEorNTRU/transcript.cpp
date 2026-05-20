/**
 * @file transcript.cpp
 * @brief Sign/Aggregate 转录本编码实现
 */

#include "../include_RLWEorNTRU/transcript.h"

#include <stdexcept>

namespace ibags {

// ============================================================================
// H2 Transcript Builder — Signing Challenge
// ============================================================================

Xof build_h2_transcript(const Params& pp,
                        ByteSpan params_encoded,
                        ByteSpan pk_encoded,
                        ByteSpan alpha_encoded,
                        ByteSpan msg_encoded)
{
    Xof xof(domain::H2_CHALLENGE);

    // 1. 参数集 domain separation
    xof.absorb_params_encoding(params_encoded);

    // 2. 公钥多项式
    xof.absorb_with_length_prefix(pk_encoded);

    // 3. alpha_j 多项式
    xof.absorb_with_length_prefix(alpha_encoded);

    // 4. 消息
    xof.absorb_with_length_prefix(msg_encoded);

    // XOF 未 finalize，调用方负责 finalize + squeeze
    return xof;
}

// ============================================================================
// H3 Transcript Builder — Aggregate Coefficients
// ============================================================================

Xof build_h3_transcript(const Params& pp,
                        ByteSpan params_encoded,
                        const std::vector<std::vector<uint8_t>>& alpha_encodeds,
                        const std::vector<std::vector<uint8_t>>& t_encodeds,
                        const std::vector<std::vector<uint8_t>>& msg_encodeds)
{
    const size_t M = alpha_encodeds.size();
    // 输入一致性检查
    if (t_encodeds.size() != M || msg_encodeds.size() != M) {
        throw std::invalid_argument(
            "build_h3_transcript: list sizes mismatch");
    }

    Xof xof(domain::H3_AGG_COEFF);

    // 1. 参数集 domain separation
    xof.absorb_params_encoding(params_encoded);

    // 2. 签名者数量 (用于 transcript 结构化)
    xof.absorb_u32(static_cast<uint32_t>(M));

    // 3. M 组 (alpha_j, t_j, mu_j)
    for (size_t j = 0; j < M; ++j) {
        // 索引前缀 (结构化 transcript，确保顺序敏感)
        xof.absorb_u32(static_cast<uint32_t>(j));

        xof.absorb_with_length_prefix(
            ByteSpan(alpha_encodeds[j].data(), alpha_encodeds[j].size()));

        xof.absorb_with_length_prefix(
            ByteSpan(t_encodeds[j].data(), t_encodeds[j].size()));

        xof.absorb_with_length_prefix(
            ByteSpan(msg_encodeds[j].data(), msg_encodeds[j].size()));
    }

    return xof;
}

// ============================================================================
// 辅助函数
// ============================================================================

void transcript_absorb_poly(Xof& xof, const Poly& poly, const Params& pp)
{
    // 编码多项式
    size_t poly_bytes = static_cast<size_t>(pp.n) * COEFF_ENCODE_BYTES;
    std::vector<uint8_t> encoded(poly_bytes);
    Status s = poly_encode(poly, pp, encoded.data(), encoded.size());
    if (!s.ok()) {
        throw std::runtime_error(
            std::string("transcript_absorb_poly: encode failed: ")
            + s.message());
    }

    xof.absorb_with_length_prefix(
        ByteSpan(encoded.data(), encoded.size()));
}

void transcript_absorb_identity(Xof& xof, const Identity& id)
{
    uint8_t buf[MAX_IDENTITY_BYTES + 2];
    size_t encoded_len = 0;
    Status s = encode_identity(id, buf, &encoded_len);
    if (!s.ok()) {
        throw std::runtime_error(
            std::string("transcript_absorb_identity: encode failed: ")
            + s.message());
    }

    xof.absorb_with_length_prefix(ByteSpan(buf, encoded_len));
}

void transcript_absorb_message(Xof& xof, const Message& msg)
{
    uint8_t buf[MAX_MESSAGE_PAYLOAD_BYTES + 3];
    size_t encoded_len = 0;
    Status s = encode_message(msg, buf, &encoded_len);
    if (!s.ok()) {
        throw std::runtime_error(
            std::string("transcript_absorb_message: encode failed: ")
            + s.message());
    }

    xof.absorb_with_length_prefix(ByteSpan(buf, encoded_len));
}

} // namespace ibags