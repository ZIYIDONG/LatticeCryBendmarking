#pragma once
/**
 * @file encode.h
 * @brief 规范序列化 — 多项式、身份标识、消息的编解码
 *
 * API:
 *  - poly_encode(p, out) / poly_decode(bytes, pp) → Poly
 *    * 128-bit 系数编码 (每个系数 16 字节, big-endian, canonical [0, q))
 *    * 总输出 = n * 16 字节
 *
 *  - encode_identity(pp, out) / decode_identity(bytes) → Identity
 *    * 支持 Vehicle RID/PID 和 RSU RID/PID 四种身份类型
 *    * 自描述格式: domain_tag (1B) || length (1B) || id_bytes
 *
 *  - encode_message(pp, out) / decode_message(bytes) → Message
 *    * 支持 V2V, V2I, RSU_Broadcast 三种消息类型
 *    * payload 直接以原始字节存储
 *
 * 设计原则:
 *  - Descriptor—类型标签（枚举）作为自描述头部
 *  - 所有编码格式用字节序约定（big-endian 用于系数）
 *  - 编解码函数返回 Status，支持错误路径
 *
 * 参考:
 *  - PQClean encode/decode (Dilithium packing)
 *  - liboqs OQS_*_pack_* 风格
 */

#include "params.h"
#include "poly.h"
#include "errors.h"
#include "domain.h"

#include <cstdint>
#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

namespace ibags {

// ============================================================================
// 常量
// ============================================================================

/// 每个系数编码字节数 (128-bit)
inline constexpr int COEFF_ENCODE_BYTES = 16;

/// 最大身份字节长度
inline constexpr size_t MAX_IDENTITY_BYTES = 256;

/// 最大消息载荷字节长度
inline constexpr size_t MAX_MESSAGE_PAYLOAD_BYTES = 65535; // 64 KiB

// ============================================================================
// Descriptor — 身份/消息类型标签
// ============================================================================

/// 身份描述符类型
enum class IdentityType : uint8_t {
    VEHICLE_RID = 0x01,
    VEHICLE_PID = 0x02,
    RSU_RID     = 0x11,
    RSU_PID     = 0x12,
};

/// 消息描述符类型
enum class MessageType : uint8_t {
    V2V           = 0x01,
    V2I           = 0x02,
    RSU_BROADCAST = 0x11,
};

// ============================================================================
// Identity — 带类型的身份
// ============================================================================

struct Identity {
    IdentityType type;
    std::string  id;  // 原始身份字节串
};

// ============================================================================
// Message — 带类型的消息
// ============================================================================

struct Message {
    MessageType type;
    std::string payload; // 原始消息载荷
};

// ============================================================================
// Polynomial Encode/Decode
// ============================================================================

/**
 * @brief 将多项式编码为字节数组
 *
 * 格式: n × 16B big-endian canonical coefficients
 * 输出缓冲区大小: pp.n * COEFF_ENCODE_BYTES
 *
 * @param poly   待编码多项式 (必须 canonical)
 * @param pp     参数集
 * @param out    输出缓冲区
 * @param out_len 输出缓冲区长度
 * @return       Ok() 或错误
 */
Status poly_encode(const Poly& poly,
                   const Params& pp,
                   uint8_t* out,
                   size_t out_len);

/**
 * @brief 从字节数组解码多项式
 *
 * @param encoded 输入字节 (长度须 ≥ pp.n * COEFF_ENCODE_BYTES)
 * @param len     输入长度
 * @param pp      参数集
 * @param out     输出多项式
 * @return        Ok() 或错误
 */
Status poly_decode(const uint8_t* encoded,
                   size_t len,
                   const Params& pp,
                   Poly* out);

// ============================================================================
// Identity Encode/Decode
// ============================================================================

/**
 * @brief 编码身份
 *
 * 格式: domain_tag (1B) || len (1B) || id_bytes
 * 其中 domain_tag 同时存储 IdentityType 和 domain label 信息
 *
 * @param id      身份
 * @param out     输出缓冲区
 * @param out_len 输出缓冲区长度
 * @return        Ok() 或实际写入字节数(通过输出参数)或错误
 */
Status encode_identity(const Identity& id,
                       uint8_t* out,
                       size_t* out_len);

/**
 * @brief 解码身份
 *
 * @param encoded  输入字节
 * @param len      输入长度
 * @param out      输出身份
 * @return         Ok() 或错误
 */
Status decode_identity(const uint8_t* encoded,
                       size_t len,
                       Identity* out);

// ============================================================================
// Message Encode/Decode
// ============================================================================

/**
 * @brief 编码消息
 *
 * 格式: domain_tag (1B) || payload_len (2B, LE) || payload
 *
 * @param msg     消息
 * @param out     输出缓冲区
 * @param out_len 输出缓冲区长度
 * @return        Ok() 或错误
 */
Status encode_message(const Message& msg,
                      uint8_t* out,
                      size_t* out_len);

/**
 * @brief 解码消息
 *
 * @param encoded 输入字节
 * @param len     输入长度
 * @param out     输出消息
 * @return        Ok() 或错误
 */
Status decode_message(const uint8_t* encoded,
                      size_t len,
                      Message* out);

} // namespace ibags