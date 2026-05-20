#pragma once
/**
 * @file xof.h
 * @brief SHAKE256 eXtendable Output Function (XOF) 包装器
 *
 * API:
 *  - absorb():       absorb byte strings (domain separation, inputs)
 *  - absorb_u16/64(): absorb integer fields
 *  - finalize():    transition from absorb to squeeze phase
 *  - squeeze(size): output variable-length bytes
 *  - squeeze_u8():  output single uint8_t
 *  - squeeze_u64(): output single uint64_t
 *
 * 后端：
 *  - PoC 使用 OpenSSL EVP SHAKE256 (EVP_MD_CTX + EVP_DigestInit/Update/FinalXOF)
 *  - 可替换为 Intel Cryptography Primitives SHAKE 或其他后端
 *
 * 典型用法 (H1 hash-to-ring):
 *  @code
 *    Xof xof(domain::H1_TO_RING);
 *    xof.absorb(inputs);
 *    xof.finalize();
 *    for each coeff:
 *      do { c = xof.squeeze_u64() % q; } while (c >= q);  // rejection
 *  @endcode
 *
 * 参考：
 *  - OpenSSL EVP_MD-SHAKE 文档
 *  - Dilithium shake256.h
 *  - PQClean shake API
 */

#include "domain.h"

#include <cstdint>
#include <cstddef>
#include <string_view>
#include <vector>
#include <string>

// 前向声明 OpenSSL 类型（无需在头文件中暴露 OpenSSL header）
struct evp_md_ctx_st;

namespace ibags {

// ============================================================================
// ByteSpan — 轻量级字节序列视图
// ============================================================================

/// 与 std::span<const uint8_t> 语义一致的字节视图
/// PoC: 避免 C++20 std::span 在某些环境缺失的问题
struct ByteSpan {
    const uint8_t* data;
    size_t         size;

    ByteSpan() : data(nullptr), size(0) {}
    ByteSpan(const uint8_t* d, size_t s) : data(d), size(s) {}

    // 从 string_view 隐式构造（用于 domain label）
    ByteSpan(std::string_view sv)
        : data(reinterpret_cast<const uint8_t*>(sv.data()))
        , size(sv.size()) {}

    // 从 string 隐式构造
    ByteSpan(const std::string& s)
        : data(reinterpret_cast<const uint8_t*>(s.data()))
        , size(s.size()) {}

    // 从 vector<uint8_t> 隐式构造
    ByteSpan(const std::vector<uint8_t>& v)
        : data(v.data()), size(v.size()) {}

    [[nodiscard]] bool empty() const noexcept { return size == 0; }
};

// ============================================================================
// Xof — SHAKE256 可扩展输出函数
// ============================================================================

class Xof {
public:
    /// 构造 SHAKE256 XOF，自动 absorb domain label
    explicit Xof(std::string_view domain_label);

    ~Xof();

    // 禁止拷贝
    Xof(const Xof&) = delete;
    Xof& operator=(const Xof&) = delete;

    // 允许移动
    Xof(Xof&&) noexcept;
    Xof& operator=(Xof&&) noexcept;

    // ── Absorb 阶段 ──
    // 调用 finalize() 前可多次 absorb

    /// absorb 字节序列
    void absorb(ByteSpan bytes);

    /// absorb 带长度前缀的字节序列 (length || bytes)
    /// length 使用 uint32_t little-endian 编码
    void absorb_with_length_prefix(ByteSpan bytes);

    /// absorb uint8_t
    void absorb_u8(uint8_t value);

    /// absorb uint16_t (little-endian)
    void absorb_u16(uint16_t value);

    /// absorb uint32_t (little-endian)
    void absorb_u32(uint32_t value);

    /// absorb uint64_t (little-endian)
    void absorb_u64(uint64_t value);

    /// absorb 参数集规范编码 (用于 domain separation)
    void absorb_params_encoding(ByteSpan params_encoded);

    // ── Finalize ──

    /// 从 absorb 阶段切换到 squeeze 阶段
    /// 在第一次 squeeze() 前必须调用
    void finalize();

    // ── Squeeze 阶段 ──

    /// squeeze 指定字节数
    /// 必须在 finalize() 后调用
    void squeeze(uint8_t* out, size_t out_len);

    /// squeeze 到 vector (便捷方法)
    std::vector<uint8_t> squeeze(size_t out_len);

    /// squeeze 单个 uint8_t
    uint8_t squeeze_u8();

    /// squeeze 单个 uint64_t
    uint64_t squeeze_u64();

private:
    evp_md_ctx_st* ctx_;     // OpenSSL EVP_MD_CTX (PIMPL)
    bool            finalized_;
    bool            moved_from_;  // 被移动状态
};

} // namespace ibags