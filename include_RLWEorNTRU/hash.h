#pragma once
/**
 * @file hash.h
 * @brief 核心哈希函数 — H1, H2, H3, hash_id_mask
 *
 * 四个核心哈希原语:
 *
 *  H1: hash-to-ring (伪名生成)
 *    t_j = H1(RID_j, pp)
 *    输入: Vehicle RID + 参数集
 *    输出: 环元素 (多项式 t_j ∈ R_q)
 *
 *  hash_id_mask: 签名者身份掩码
 *    id_mask = H_idmask(PID_j, pp)
 *    输入: Vehicle PID + 参数集
 *    输出: 多项式 m_j ∈ R_q (用于签名多项式构造)
 *
 *  H2: 签名挑战 (Fiat-Shamir)
 *    c_j = H2(PK_Vj, alpha_j, mu_j, pp)
 *    输入: 公钥 + alpha 承诺 + 消息 + 参数集
 *    输出: 挑战多项式 c_j ∈ C (二进制/稀疏多项式)
 *
 *  H3: 聚合系数
 *    beta_j = H3({ (alpha_j, t_j, mu_j) }_{j=1..M}, pp)
 *    输入: 所有签名者的承诺+伪名+消息
 *    输出: 聚合系数 beta_j ∈ Z_q (或小系数多项式)
 *
 * 设计原则:
 *  - 所有哈希使用 SHAKE256 XOF (通过 xof.h)
 *  - Domain separation: 每个函数使用独立的 domain label
 *  - 输出: H1/H2/H3 产生环元素/多项式，保证在边界内
 *  - 无外部随机源依赖 (确定性，仅依赖 XOF)
 *
 * 参考:
 *  - iBAGS specification: hash-to-ring, Fiat-Shamir, aggregation
 *  - Dilithium: shake256 expand_a, expand_mask, sample_in_ball
 *  - Falcon: hash_to_point (SHAKE-256 → ring element)
 */

#include "xof.h"
#include "csprng.h"
#include "encode.h"
#include "transcript.h"
#include "domain.h"
#include "params.h"
#include "poly.h"
#include "errors.h"

#include <cstdint>
#include <cstddef>
#include <string_view>
#include <vector>

namespace ibags {

// ============================================================================
// H1: hash-to-ring — 伪名生成
// ============================================================================

/**
 * @brief H1 hash-to-ring: 从 RID 生成伪名环元素
 *
 * 算法概述:
 *   1. 构造 XOF: domain = H1_TO_RING
 *   2. absorb: encode_params(pp) || encode_identity(RID)
 *   3. finalize
 *   4. rejection-sample to ring elements until k valid P_k are produced
 *
 * 环元素采样策略:
 *   - 每个系数使用拒绝采样: squeeze 64-bit → modulo q → accept if < q
 *   - 输出多项式 t ∈ R_q
 *
 * @param pp              参数集
 * @param params_encoded  参数集规范编码
 * @param rid             车辆真实身份 (RID)
 * @return                伪名多项式 t
 */
Poly hash_to_ring(const Params& pp,
                  ByteSpan params_encoded,
                  const Identity& rid);

// ============================================================================
// hash_id_mask: 签名者身份掩码
// ============================================================================

/**
 * @brief 从 PID 生成签名者身份掩码多项式
 *
 * 算法:
 *   1. 构造 XOF: domain = H_ID_MASK
 *   2. absorb: encode_params(pp) || encode_identity(PID)
 *   3. finalize
 *   4. squeeze n 个系数 ∈ R_q
 *
 * @param pp              参数集
 * @param params_encoded  参数集规范编码
 * @param pid             车辆假名 (PID)
 * @return                掩码多项式 m
 */
Poly hash_id_mask(const Params& pp,
                  ByteSpan params_encoded,
                  const Identity& pid);

// ============================================================================
// H2: 签名挑战 (Fiat-Shamir)
// ============================================================================

/**
 * @brief H2 签名挑战: 从 transcript 生成挑战多项式
 *
 * 挑战多项式空间 C = { c ∈ R_q : ‖c‖_∞ ≤ 1, Hamming_weight(c) ≤ κ }
 *
 * 算法:
 *   1. 使用 build_h2_transcript 构造 XOF
 *   2. finalize
 *   3. rejection-sample 挑战多项式:
 *      - 重复 squeeze u64 直到获得 κ 个不同的非零位置
 *      - 每个非零位置随机选择 ±1
 *
 * 挑战多项式规范:
 *   c = Σ_{pos} sign(pos) · X^{pos}
 *
 * @param pp              参数集
 * @param params_encoded  参数集规范编码
 * @param pk_encoded      公钥多项式编码
 * @param alpha_encoded   alpha 承诺编码
 * @param msg_encoded     消息编码
 * @return                挑战多项式 c
 */
Poly hash_challenge(const Params& pp,
                    ByteSpan params_encoded,
                    ByteSpan pk_encoded,
                    ByteSpan alpha_encoded,
                    ByteSpan msg_encoded);

// ============================================================================
// H3: 聚合系数
// ============================================================================

/**
 * @brief H3 聚合系数: 从 transcript 生成范数有界的聚合权重
 *
 * 算法:
 *   1. 使用 build_h3_transcript 构造 XOF
 *   2. finalize
 *   3. 对每个签名者 j 生成聚合系数 beta_j ∈ Z_q
 *
 * 输出: M 个标量, 每个 ∈ [0, q)
 *
 * @param pp              参数集
 * @param params_encoded  参数集规范编码
 * @param alpha_encodeds  alpha_j 编码列表 (长度 M)
 * @param t_encodeds      t_j 编码列表 (长度 M)
 * @param msg_encodeds    mu_j 编码列表 (长度 M)
 * @return                M 个聚合系数
 */
std::vector<uint64_t> hash_agg_coeffs(const Params& pp,
                                       ByteSpan params_encoded,
                                       const std::vector<std::vector<uint8_t>>& alpha_encodeds,
                                       const std::vector<std::vector<uint8_t>>& t_encodeds,
                                       const std::vector<std::vector<uint8_t>>& msg_encodeds);

// ============================================================================
// 采样工具函数
// ============================================================================

/**
 * @brief 从 XOF 产生环元素的所有系数 (rejection sampling per coefficient)
 *
 * 每个系数通过 squeeze_u64() % q 产生，拒绝 ≥ q 的值
 *
 * @param xof     已 finalized 的 XOF
 * @param pp      参数集
 * @param out     输出多项式
 */
void sample_ring_element(Xof& xof, const Params& pp, Poly* out);

/**
 * @brief 从 XOF 产生稀疏/二进制多项式 (Fiat-Shamir 挑战)
 *
 * 产生 κ 个非零位置 (±1)，其余位置为 0
 *
 * @param xof     已 finalized 的 XOF
 * @param pp      参数集
 * @param kappa   非零系数数
 * @param out     输出多项式
 */
void sample_challenge_polynomial(Xof& xof, const Params& pp,
                                  int kappa, Poly* out);

} // namespace ibags