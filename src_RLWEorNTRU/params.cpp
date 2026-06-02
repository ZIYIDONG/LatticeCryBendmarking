/**
 * @file params.cpp
 * @brief IBAGS 参数集模块实现
 */

#include "../include_RLWEorNTRU/params.h"
#include <cstring>
#include <cmath>
#include <algorithm>

namespace ibags {

// ============================================================================
// 内部常量
// ============================================================================

/// 最小环维度（低于此值视为不安全）
inline constexpr int MIN_RING_DIM = 64;

/// 最大模数（防止溢出且限制参数膨胀）
inline constexpr int MAX_MODULUS = (1 << 30);   // 2^30

/// 最大签名者数（协议实际限制）
inline constexpr int MAX_SIGNERS = 256;

// ============================================================================
// §1  ParamId 辅助函数
// ============================================================================

const char* param_id_name(ParamId id) noexcept {
    switch (id) {
        case ParamId::IBAGS_64_DEMO:     return "IBAGS-64-Demo";
        case ParamId::IBAGS_512_LEVEL2:  return "IBAGS-512-Level2";
        case ParamId::IBAGS_1024_LEVEL3: return "IBAGS-1024-Level3";
        case ParamId::IBAGS_1024_LEVEL5: return "IBAGS-1024-Level5";
        default:                          return "IBAGS-Unknown";
    }
}

// ============================================================================
// §2  Params 工厂函数
// ============================================================================

/// NIST 安全级别映射 (签名方案对标 Dilithium)
///   LATTICE_LEVEL_1  → IBAGS Level 2 (n=512, Dilithium2, 签名无 Level 1)
///   LATTICE_LEVEL_3  → IBAGS Level 3 (n=1024, Dilithium3)
///   LATTICE_LEVEL_5  → IBAGS Level 5 (n=1024, Dilithium5)
///   无定义           → Demo (n=64)
Params default_params() {
#if defined(LATTICE_LEVEL_5)
    return Params::params_level5_1024();
#elif defined(LATTICE_LEVEL_3)
    return Params::params_level3_1024();
#elif defined(LATTICE_LEVEL_1)
    return Params::params_level2_512();
#else  // Demo
    return Params::params_demo_64();
#endif
}

Params Params::params_demo_64() {
    Params p;
    p.param_id            = ParamId::IBAGS_64_DEMO;
    p.n                   = 64;
    p.q                   = 7681;            // 13-bit NTT-friendly prime (7681 ≡ 1 mod 128)
    p.sigma               = 1.0;
    p.eta1                = 4;
    p.eta2                = 2;
    p.eta_ver             = 2;
    p.max_signers         = 8;
    p.max_rejection_loops = 10;
    p.kappa               = 8;
    p.coefficient_bytes   = 2;               // ceil(log2(8191)/8) = 2
    p.poly_bytes          = p.n * p.coefficient_bytes;
    return p;
}

Params Params::params_level2_512() {
    Params p;
    p.param_id            = ParamId::IBAGS_512_LEVEL2;
    p.n                   = 512;
    p.q                   = 8404993;         // 24-bit NTT-friendly prime (≡ 1 mod 1024)
    p.sigma               = 1.0;
    p.eta1                = 40;
    p.eta2                = 20;
    p.eta_ver             = 20;
    p.max_signers         = 32;
    p.max_rejection_loops = 20;
    p.kappa               = 60;              // τ: challenge Hamming weight
    p.coefficient_bytes   = 3;               // ceil(log2(8404993)/8) = 3
    p.poly_bytes          = p.n * p.coefficient_bytes;
    return p;
}

Params Params::params_level3_1024() {
    Params p;
    p.param_id            = ParamId::IBAGS_1024_LEVEL3;
    p.n                   = 1024;
    p.q                   = 16900097;        // 25-bit NTT-friendly prime (16900097 ≡ 1 mod 2048, verified prime)
    p.sigma               = 1.0;
    p.eta1                = 55;
    p.eta2                = 28;
    p.eta_ver             = 28;
    p.max_signers         = 48;
    p.max_rejection_loops = 25;
    p.kappa               = 80;              // τ: challenge Hamming weight
    p.coefficient_bytes   = 4;               // ceil(log2(16777259)/8) = 4
    p.poly_bytes          = p.n * p.coefficient_bytes;
    return p;
}

Params Params::params_level5_1024() {
    Params p;
    p.param_id            = ParamId::IBAGS_1024_LEVEL5;
    p.n                   = 1024;
    p.q                   = 4206593;         // 23-bit NTT-friendly prime (≡ 1 mod 2048)
    p.sigma               = 1.0;
    p.eta1                = 60;
    p.eta2                = 30;
    p.eta_ver             = 30;
    p.max_signers         = 64;
    p.max_rejection_loops = 30;
    p.kappa               = 100;             // τ: challenge Hamming weight
    p.coefficient_bytes   = 3;               // ceil(log2(4206593)/8) = 3
    p.poly_bytes          = p.n * p.coefficient_bytes;
    return p;
}

// ============================================================================
// §3  参数验证
// ============================================================================

Status validate_params(const Params& p) {
    // ── 1. n 必须为 2 的幂且 ≥ MIN_RING_DIM ──
    if (p.n < MIN_RING_DIM || (p.n & (p.n - 1)) != 0) {
        return {ErrorCode::InvalidParams,
                "n must be a power of 2 and >= " + std::to_string(MIN_RING_DIM)};
    }

    // ── 2. q 必须为正数且 ≤ MAX_MODULUS ──
    if (p.q <= 0 || p.q > MAX_MODULUS) {
        return {ErrorCode::InvalidParams,
                "q must be in (0, " + std::to_string(MAX_MODULUS) + "]"};
    }

    // ── 3. sigma 必须为正 ──
    if (p.sigma <= 0.0) {
        return {ErrorCode::InvalidParams, "sigma must be positive"};
    }

    // ── 4. eta 链合理性: 1 ≤ eta_ver ≤ eta2 ≤ eta1 ──
    if (p.eta_ver < 1 || p.eta2 < p.eta_ver || p.eta1 < p.eta2) {
        return {ErrorCode::InvalidParams,
                "eta values must satisfy: 1 <= eta_ver <= eta2 <= eta1"};
    }

    // ── 5. max_signers ∈ [1, MAX_SIGNERS] ──
    if (p.max_signers < 1 || p.max_signers > MAX_SIGNERS) {
        return {ErrorCode::InvalidParams,
                "max_signers must be in [1, " + std::to_string(MAX_SIGNERS) + "]"};
    }

    // ── 6. max_rejection_loops 必须为正 ──
    if (p.max_rejection_loops <= 0) {
        return {ErrorCode::InvalidParams,
                "max_rejection_loops must be positive"};
    }

    // ── 7. kappa 必须为正且 ≤ n ──
    if (p.kappa <= 0 || p.kappa > p.n) {
        return {ErrorCode::InvalidParams,
                "kappa must be in (0, n]"};
    }

    // ── 8. coefficient_bytes 必须足够容纳 mod q ──
    int required_bytes = 0;
    long tmp = p.q - 1;               // 最大值 = q - 1
    while (tmp > 0) {
        required_bytes++;
        tmp >>= 8;
    }
    if (p.coefficient_bytes < required_bytes) {
        return {ErrorCode::InvalidParams,
                "coefficient_bytes too small: need " + std::to_string(required_bytes) +
                " to represent mod " + std::to_string(p.q)};
    }

    // ── 9. poly_bytes 一致性检查 ──
    if (p.poly_bytes != p.n * p.coefficient_bytes) {
        return {ErrorCode::InvalidParams,
                "poly_bytes must equal n * coefficient_bytes"};
    }

    return Status::Ok();
}

// ============================================================================
// §4  参数编码 (Domain Separation)
// ============================================================================

Status encode_params(const Params& p,
                     uint8_t* out,
                     size_t   out_len) noexcept {
    if (out == nullptr) {
        return {ErrorCode::InvalidEncoding, "encode_params: out is null"};
    }
    if (out_len < PARAMS_ENCODE_BYTES) {
        return {ErrorCode::InvalidEncoding,
                "encode_params: out_len < PARAMS_ENCODE_BYTES"};
    }

    // 先清零整个缓冲区
    std::memset(out, 0, PARAMS_ENCODE_BYTES);

    // 小端序写入各字段
    size_t pos = 0;

    // Bytes 0–1: param_id (uint16_t LE)
    uint16_t pid = static_cast<uint16_t>(p.param_id);
    out[pos++] = static_cast<uint8_t>(pid & 0xFF);
    out[pos++] = static_cast<uint8_t>((pid >> 8) & 0xFF);

    // Bytes 2–3: n (uint16_t LE)
    uint16_t n_val = static_cast<uint16_t>(p.n);
    out[pos++] = static_cast<uint8_t>(n_val & 0xFF);
    out[pos++] = static_cast<uint8_t>((n_val >> 8) & 0xFF);

    // Bytes 4–7: q (uint32_t LE)
    uint32_t q_val = static_cast<uint32_t>(p.q);
    out[pos++] = static_cast<uint8_t>(q_val & 0xFF);
    out[pos++] = static_cast<uint8_t>((q_val >> 8) & 0xFF);
    out[pos++] = static_cast<uint8_t>((q_val >> 16) & 0xFF);
    out[pos++] = static_cast<uint8_t>((q_val >> 24) & 0xFF);

    // Bytes 8–15: sigma (double, IEEE 754 LE)
    {
        uint64_t sigma_bits;
        static_assert(sizeof(double) == sizeof(uint64_t),
                      "double must be 8 bytes");
        std::memcpy(&sigma_bits, &p.sigma, sizeof(sigma_bits));
        for (int i = 0; i < 8; ++i) {
            out[pos++] = static_cast<uint8_t>(sigma_bits & 0xFF);
            sigma_bits >>= 8;
        }
    }

    // Bytes 16–17: eta1
    uint16_t e1 = static_cast<uint16_t>(p.eta1);
    out[pos++] = static_cast<uint8_t>(e1 & 0xFF);
    out[pos++] = static_cast<uint8_t>((e1 >> 8) & 0xFF);

    // Bytes 18–19: eta2
    uint16_t e2 = static_cast<uint16_t>(p.eta2);
    out[pos++] = static_cast<uint8_t>(e2 & 0xFF);
    out[pos++] = static_cast<uint8_t>((e2 >> 8) & 0xFF);

    // Bytes 20–21: eta_ver
    uint16_t ev = static_cast<uint16_t>(p.eta_ver);
    out[pos++] = static_cast<uint8_t>(ev & 0xFF);
    out[pos++] = static_cast<uint8_t>((ev >> 8) & 0xFF);

    // Bytes 22–23: max_signers
    uint16_t ms = static_cast<uint16_t>(p.max_signers);
    out[pos++] = static_cast<uint8_t>(ms & 0xFF);
    out[pos++] = static_cast<uint8_t>((ms >> 8) & 0xFF);

    // Bytes 24–25: max_rejection_loops
    uint16_t mr = static_cast<uint16_t>(p.max_rejection_loops);
    out[pos++] = static_cast<uint8_t>(mr & 0xFF);
    out[pos++] = static_cast<uint8_t>((mr >> 8) & 0xFF);

    // Bytes 26–27: coefficient_bytes
    uint16_t cb = static_cast<uint16_t>(p.coefficient_bytes);
    out[pos++] = static_cast<uint8_t>(cb & 0xFF);
    out[pos++] = static_cast<uint8_t>((cb >> 8) & 0xFF);

    // Bytes 28–29: kappa
    uint16_t k_val = static_cast<uint16_t>(p.kappa);
    out[pos++] = static_cast<uint8_t>(k_val & 0xFF);
    out[pos++] = static_cast<uint8_t>((k_val >> 8) & 0xFF);

    // Bytes 30–63: 保留 (已由 memset 清零)

    return Status::Ok();
}

} // namespace ibags