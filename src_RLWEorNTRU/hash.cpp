/**
 * @file hash.cpp
 * @brief 核心哈希函数实现 — H1, H2, H3, hash_id_mask
 */

#include "../include_RLWEorNTRU/hash.h"

#include <cstring>
#include <set>
#include <stdexcept>

namespace ibags {

// ============================================================================
// 采样工具函数
// ============================================================================

void sample_ring_element(Xof& xof, const Params& pp, Poly* out)
{
    if (!out) return;

    std::vector<int64_t> coeffs(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        // 使用模约减 (q << 2^64, 偏倚可忽略)
        uint64_t r = xof.squeeze_u64();
        coeffs[i] = static_cast<int64_t>(r % static_cast<uint64_t>(pp.q));
    }
    *out = Poly(pp.n, std::move(coeffs));
}

void sample_challenge_polynomial(Xof& xof, const Params& pp,
                                  int kappa, Poly* out)
{
    if (!out) return;

    // 初始化为零多项式
    std::vector<int64_t> coeffs(pp.n, 0);

    if (kappa <= 0) {
        *out = Poly(pp.n, std::move(coeffs));
        return;
    }
    if (kappa > pp.n) {
        kappa = pp.n;
    }

    // 使用 set 确保位置唯一
    std::set<int> positions;
    while (static_cast<int>(positions.size()) < kappa) {
        uint64_t r = xof.squeeze_u64();
        int pos = static_cast<int>(r % static_cast<uint64_t>(pp.n));
        if (positions.find(pos) == positions.end()) {
            positions.insert(pos);

            // 随机选择 ±1 符号
            uint8_t sign_byte = xof.squeeze_u8();
            int64_t sign = (sign_byte & 1) ? 1 : -1;
            coeffs[pos] = sign;
        }
    }

    *out = Poly(pp.n, std::move(coeffs));
}

// ============================================================================
// H1: hash-to-ring — 伪名生成
// ============================================================================

Poly hash_to_ring(const Params& pp,
                  ByteSpan params_encoded,
                  const Identity& rid)
{
    // 编码身份
    uint8_t id_buf[MAX_IDENTITY_BYTES + 2];
    size_t id_encoded_len = 0;
    Status s = encode_identity(rid, id_buf, &id_encoded_len);
    if (!s.ok()) {
        throw std::runtime_error(
            std::string("hash_to_ring: encode_identity failed: ")
            + s.message());
    }

    // 构建 XOF
    Xof xof(domain::H1_TO_RING);

    // absorb: params || identity
    xof.absorb_params_encoding(params_encoded);
    xof.absorb_with_length_prefix(ByteSpan(id_buf, id_encoded_len));

    xof.finalize();

    // 采样环元素
    Poly t = Poly::zero(pp);
    sample_ring_element(xof, pp, &t);
    return t;
}

// ============================================================================
// hash_id_mask: 签名者身份掩码
// ============================================================================

Poly hash_id_mask(const Params& pp,
                  ByteSpan params_encoded,
                  const Identity& pid)
{
    // 编码身份
    uint8_t id_buf[MAX_IDENTITY_BYTES + 2];
    size_t id_encoded_len = 0;
    Status s = encode_identity(pid, id_buf, &id_encoded_len);
    if (!s.ok()) {
        throw std::runtime_error(
            std::string("hash_id_mask: encode_identity failed: ")
            + s.message());
    }

    // 构建 XOF
    Xof xof(domain::H1_ID_MASK);

    // absorb: params || identity
    xof.absorb_params_encoding(params_encoded);
    xof.absorb_with_length_prefix(ByteSpan(id_buf, id_encoded_len));

    xof.finalize();

    // 采样环元素
    Poly m = Poly::zero(pp);
    sample_ring_element(xof, pp, &m);
    return m;
}

// ============================================================================
// H2: 签名挑战 (Fiat-Shamir)
// ============================================================================

Poly hash_challenge(const Params& pp,
                    ByteSpan params_encoded,
                    ByteSpan pk_encoded,
                    ByteSpan alpha_encoded,
                    ByteSpan msg_encoded)
{
    // 使用 transcript builder 构造 XOF
    Xof xof = build_h2_transcript(
        pp, params_encoded, pk_encoded, alpha_encoded, msg_encoded);

    xof.finalize();

    // 采样挑战多项式 (稀疏 ±1)
    Poly c = Poly::zero(pp);
    sample_challenge_polynomial(xof, pp, pp.kappa, &c);
    return c;
}

// ============================================================================
// H3: 聚合系数
// ============================================================================

std::vector<uint64_t> hash_agg_coeffs(
    const Params& pp,
    ByteSpan params_encoded,
    const std::vector<std::vector<uint8_t>>& alpha_encodeds,
    const std::vector<std::vector<uint8_t>>& t_encodeds,
    const std::vector<std::vector<uint8_t>>& msg_encodeds)
{
    const size_t M = alpha_encodeds.size();
    if (t_encodeds.size() != M || msg_encodeds.size() != M) {
        throw std::invalid_argument(
            "hash_agg_coeffs: list sizes mismatch");
    }

    // 使用 transcript builder 构造 XOF
    Xof xof = build_h3_transcript(
        pp, params_encoded, alpha_encodeds, t_encodeds, msg_encodeds);

    xof.finalize();

    // 为每个签名者生成聚合系数
    std::vector<uint64_t> betas(M);
    for (size_t j = 0; j < M; ++j) {
        // 拒绝采样: squeeze → mod q → accept if < q
        for (;;) {
            uint64_t r = xof.squeeze_u64();
            if (r < static_cast<uint64_t>(pp.q)) {
                betas[j] = r;
                break;
            }
        }
    }

    return betas;
}

} // namespace ibags