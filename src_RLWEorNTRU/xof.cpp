/**
 * @file xof.cpp
 * @brief SHAKE256 XOF 包装器实现 (OpenSSL EVP 后端)
 */

#include "../include_RLWEorNTRU/xof.h"

#include <openssl/evp.h>

#include <cstring>
#include <stdexcept>

namespace ibags {

// ============================================================================
// 构造/析构/移动
// ============================================================================

Xof::Xof(std::string_view domain_label)
    : ctx_(nullptr)
    , finalized_(false)
    , moved_from_(false)
{
    ctx_ = EVP_MD_CTX_new();
    if (!ctx_) {
        throw std::runtime_error("Xof: EVP_MD_CTX_new failed");
    }

    int rc = EVP_DigestInit_ex(
        static_cast<EVP_MD_CTX*>(ctx_),
        EVP_shake256(),
        nullptr);
    if (rc != 1) {
        EVP_MD_CTX_free(static_cast<EVP_MD_CTX*>(ctx_));
        ctx_ = nullptr;
        throw std::runtime_error("Xof: EVP_DigestInit_ex(SHAKE256) failed");
    }

    // 立即 absorb domain label
    absorb(ByteSpan(domain_label));
}

Xof::~Xof()
{
    if (ctx_ && !moved_from_) {
        EVP_MD_CTX_free(static_cast<EVP_MD_CTX*>(ctx_));
    }
}

Xof::Xof(Xof&& other) noexcept
    : ctx_(other.ctx_)
    , finalized_(other.finalized_)
    , moved_from_(false)
{
    other.moved_from_ = true;
}

Xof& Xof::operator=(Xof&& other) noexcept
{
    if (this != &other) {
        if (ctx_ && !moved_from_) {
            EVP_MD_CTX_free(static_cast<EVP_MD_CTX*>(ctx_));
        }
        ctx_        = other.ctx_;
        finalized_  = other.finalized_;
        moved_from_ = false;
        other.moved_from_ = true;
    }
    return *this;
}

// ============================================================================
// Absorb 阶段
// ============================================================================

void Xof::absorb(ByteSpan bytes)
{
    if (finalized_) {
        throw std::runtime_error("Xof::absorb: already finalized");
    }
    if (bytes.size == 0) return;

    int rc = EVP_DigestUpdate(
        static_cast<EVP_MD_CTX*>(ctx_),
        bytes.data,
        bytes.size);
    if (rc != 1) {
        throw std::runtime_error("Xof::absorb: EVP_DigestUpdate failed");
    }
}

void Xof::absorb_with_length_prefix(ByteSpan bytes)
{
    uint32_t len = static_cast<uint32_t>(bytes.size);
    absorb_u32(len);
    absorb(bytes);
}

void Xof::absorb_u8(uint8_t value)
{
    uint8_t buf[1] = { value };
    absorb(ByteSpan(buf, 1));
}

void Xof::absorb_u16(uint16_t value)
{
    uint8_t buf[2];
    buf[0] = static_cast<uint8_t>(value & 0xFF);
    buf[1] = static_cast<uint8_t>((value >> 8) & 0xFF);
    absorb(ByteSpan(buf, 2));
}

void Xof::absorb_u32(uint32_t value)
{
    uint8_t buf[4];
    buf[0] = static_cast<uint8_t>(value & 0xFF);
    buf[1] = static_cast<uint8_t>((value >> 8) & 0xFF);
    buf[2] = static_cast<uint8_t>((value >> 16) & 0xFF);
    buf[3] = static_cast<uint8_t>((value >> 24) & 0xFF);
    absorb(ByteSpan(buf, 4));
}

void Xof::absorb_u64(uint64_t value)
{
    uint8_t buf[8];
    buf[0] = static_cast<uint8_t>(value & 0xFF);
    buf[1] = static_cast<uint8_t>((value >> 8) & 0xFF);
    buf[2] = static_cast<uint8_t>((value >> 16) & 0xFF);
    buf[3] = static_cast<uint8_t>((value >> 24) & 0xFF);
    buf[4] = static_cast<uint8_t>((value >> 32) & 0xFF);
    buf[5] = static_cast<uint8_t>((value >> 40) & 0xFF);
    buf[6] = static_cast<uint8_t>((value >> 48) & 0xFF);
    buf[7] = static_cast<uint8_t>((value >> 56) & 0xFF);
    absorb(ByteSpan(buf, 8));
}

void Xof::absorb_params_encoding(ByteSpan params_encoded)
{
    absorb_with_length_prefix(params_encoded);
}

// ============================================================================
// Finalize
// ============================================================================

void Xof::finalize()
{
    finalized_ = true;
    // SHAKE256 无需显式 finalize (可由 squeeze 隐式触发)
    // 设置标志防止后续 absorb
}

// ============================================================================
// Squeeze 阶段
// ============================================================================

void Xof::squeeze(uint8_t* out, size_t out_len)
{
    if (!finalized_) {
        throw std::runtime_error("Xof::squeeze: not finalized");
    }
    if (out_len == 0) return;

    // 创建临时 context 副本进行 squeeze，保留原始 context 用于后续 squeeze
    // OpenSSL 的 EVP_DigestFinalXOF 会消耗 context，所以使用副本
    EVP_MD_CTX* ctx_copy = EVP_MD_CTX_new();
    if (!ctx_copy) {
        throw std::runtime_error("Xof::squeeze: EVP_MD_CTX_new failed");
    }

    int rc = EVP_MD_CTX_copy_ex(ctx_copy, static_cast<EVP_MD_CTX*>(ctx_));
    if (rc != 1) {
        EVP_MD_CTX_free(ctx_copy);
        throw std::runtime_error("Xof::squeeze: EVP_MD_CTX_copy_ex failed");
    }

    rc = EVP_DigestFinalXOF(ctx_copy, out, out_len);
    if (rc != 1) {
        EVP_MD_CTX_free(ctx_copy);
        throw std::runtime_error("Xof::squeeze: EVP_DigestFinalXOF failed");
    }

    EVP_MD_CTX_free(ctx_copy);
}

std::vector<uint8_t> Xof::squeeze(size_t out_len)
{
    std::vector<uint8_t> result(out_len);
    squeeze(result.data(), out_len);
    return result;
}

uint8_t Xof::squeeze_u8()
{
    uint8_t value = 0;
    squeeze(&value, 1);
    return value;
}

uint64_t Xof::squeeze_u64()
{
    uint8_t buf[8] = {0};
    squeeze(buf, 8);
    uint64_t value = 0;
    for (int i = 0; i < 8; ++i) {
        value |= static_cast<uint64_t>(buf[i]) << (i * 8);
    }
    return value;
}

} // namespace ibags