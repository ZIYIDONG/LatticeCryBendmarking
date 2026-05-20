#pragma once
/**
 * @file domain.h
 * @brief IBAGS 协议中所有哈希域分隔标签的集中定义
 *
 * 设计原则：
 *  - 每个 cryptographic hash 用途必须有唯一的 domain separation label
 *  - 标签前缀 "IBAGS-v1/" 确保协议版本分离
 *  - 标签作为 byte string 写入 XOF/Hash 输入头部
 *  - 同一输入不同 domain 必须产生不同的哈希输出
 *
 * 参考：
 *  - Merlin/STROBE transcript domain separation
 *  - Dilithium seed expansion label 风格 (e.g., "seed-dilithium-3-ed25519")
 */

#include <cstdint>
#include <cstddef>
#include <string_view>

namespace ibags {
namespace domain {

/**
 * @brief 返回 domain label 作为 byte span 用于 XOF absorb
 *
 * 所有 label 是以常量字符串形式存储，无需堆分配。
 */

// ── H1: Hash-to-Ring for pseudonym derivation ──
// t_j = H1(PID_Vj) 的 domain
inline constexpr std::string_view H1_TO_RING = "IBAGS-v1/H1/TO-RING";

// ── H1: ID mask derivation ──
// PID = RID XOR H1_mask(a_j) 的 mask 生成 domain
// 注意: 使用独立 domain 区分 hash_to_ring_H1 和 hash_id_mask
inline constexpr std::string_view H1_ID_MASK = "IBAGS-v1/H1/ID-MASK";

// ── H2: Challenge polynomial derivation ──
// C_j = H2(PK_Vj, alpha_j, mu_j)
inline constexpr std::string_view H2_CHALLENGE = "IBAGS-v1/H2/CHALLENGE";

// ── H3: Aggregate coefficients derivation ──
// beta_1,...,beta_M = H3(alpha_1,t_1,mu_1,...,alpha_M,t_M,mu_M)
inline constexpr std::string_view H3_AGG_COEFF = "IBAGS-v1/H3/AGG-COEFF";

// ── Gaussian sampler seed expansion ──
// GaussFunction(seed, R, M) -> alpha_1,...,alpha_M
inline constexpr std::string_view GAUSS_FUNCTION = "IBAGS-v1/GAUSS-FUNCTION";

// ── Tracing record transcript ──
// 用于追踪表条目哈希
inline constexpr std::string_view TRACE_RECORD = "IBAGS-v1/TRACE-RECORD";

// ── CSPRNG seed derivation ──
// 从种子确定性派生伪随机流
inline constexpr std::string_view CSPRNG_SEED = "IBAGS-v1/CSPRNG/SEED";

// ── Sign transcript domain prefix ──
inline constexpr std::string_view SIGN_TRANSCRIPT = "IBAGS-v1/SIGN-TRANSCRIPT";

// ── Aggregate transcript domain prefix ──
inline constexpr std::string_view AGG_TRANSCRIPT = "IBAGS-v1/AGG-TRANSCRIPT";

// ── Identity encoding domain prefixes ──
inline constexpr std::string_view ID_VEHICLE_RID = "IBAGS-v1/ID/VEHICLE-RID";
inline constexpr std::string_view ID_VEHICLE_PID = "IBAGS-v1/ID/VEHICLE-PID";
inline constexpr std::string_view ID_RSU_RID      = "IBAGS-v1/ID/RSU-RID";
inline constexpr std::string_view ID_RSU_PID      = "IBAGS-v1/ID/RSU-PID";

// ── Message encoding domain prefixes ──
inline constexpr std::string_view MSG_DOMAIN_V2V          = "IBAGS-v1/MSG/V2V";
inline constexpr std::string_view MSG_DOMAIN_V2I          = "IBAGS-v1/MSG/V2I";
inline constexpr std::string_view MSG_DOMAIN_RSU_BROADCAST = "IBAGS-v1/MSG/RSU-BROADCAST";

} // namespace domain
} // namespace ibags