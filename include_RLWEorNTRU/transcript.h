#pragma once
/**
 * @file transcript.h
 * @brief Sign/Aggregate 转录本编码 — H2 和 H3 的输入数据构造
 *
 * 设计原则:
 *  - 遵循 crypto domain separation: 每条消息都在指定 domain 下序列化
 *  - 将协议中传递的复杂参数（多项式、身份等）编码为 XOF 可吸收的字节序列
 *  - 返回构建好的 XOF 对象（已吸收所有输入），调用方在外部 finalize/squeeze
 *
 * H2 输入格式 (Signing 阶段):
 *   [domain: H2_CHALLENGE]
 *   encode_params(pp)
 *   poly_encode(PK_Vj)
 *   poly_encode(alpha_j)
 *   encode_message(mu_j)
 *
 * H3 输入格式 (Aggregation 阶段):
 *   [domain: H3_AGG_COEFF]
 *   encode_params(pp)
 *   for j = 1..M:
 *     poly_encode(alpha_j)
 *     poly_encode(t_j)          // H1 输出: 伪名环元
 *     encode_message(mu_j)
 *
 * 参考:
 *  - Dilithium 的 chall_3 上下文哈希构建
 *  - Falcon 的 transcript hash (G, c0, c1)
 *  - Merlin transcript protocol (STROBE)
 */

#include "xof.h"
#include "encode.h"
#include "domain.h"
#include "params.h"
#include "poly.h"
#include "errors.h"

#include <cstdint>
#include <cstddef>
#include <vector>

namespace ibags {

// ============================================================================
// H2 Transcript Builder — Signing Challenge
// ============================================================================

/**
 * @brief 构造 H2 Signing 挑战的 XOF transcript
 *
 * 吸收顺序:
 *   1. H2_CHALLENGE domain label (在 Xof 构造时已吸收)
 *   2. Params 规范编码 (domain separation)
 *   3. PK_Vj (多项式编码)
 *   4. alpha_j (多项式编码)
 *   5. mu_j (消息编码)
 *
 * @param pp            参数集
 * @param params_encoded 参数集规范编码 (来自 encode_params)
 * @param pk_encoded     公钥多项式编码
 * @param alpha_encoded  alpha 多项式编码
 * @param msg_encoded    消息编码
 * @return               已构造的 XOF (处于 absorb 模式，未 finalized)
 *
 * 调用方随后应为每个系数调用 xof.finalize() + xof.squeeze_u64() 提取挑战
 */
Xof build_h2_transcript(const Params& pp,
                        ByteSpan params_encoded,
                        ByteSpan pk_encoded,
                        ByteSpan alpha_encoded,
                        ByteSpan msg_encoded);

// ============================================================================
// H3 Transcript Builder — Aggregate Coefficients
// ============================================================================

/**
 * @brief 构造 H3 聚合系数计算的 XOF transcript
 *
 * 吸收顺序:
 *   1. H3_AGG_COEFF domain label (在 Xof 构造时已吸收)
 *   2. Params 规范编码 (domain separation)
 *   3. M 组 (alpha_j, t_j, mu_j)
 *
 * @param pp              参数集
 * @param params_encoded   参数集规范编码
 * @param alpha_encodeds   alpha_j 多项式编码列表 (长度 M)
 * @param t_encodeds       t_j 多项式编码列表 (长度 M)
 * @param msg_encodeds     mu_j 消息编码列表 (长度 M)
 * @return                 已构造的 XOF (处于 absorb 模式，未 finalized)
 *
 * 调用方随后应为每个系数调用 xof.finalize() + xof.squeeze_u64() 提取聚合系数
 */
Xof build_h3_transcript(const Params& pp,
                        ByteSpan params_encoded,
                        const std::vector<std::vector<uint8_t>>& alpha_encodeds,
                        const std::vector<std::vector<uint8_t>>& t_encodeds,
                        const std::vector<std::vector<uint8_t>>& msg_encodeds);

// ============================================================================
// 辅助函数
// ============================================================================

/**
 * @brief 将多项式编码为字节并拼接到 XOF
 */
void transcript_absorb_poly(Xof& xof, const Poly& poly, const Params& pp);

/**
 * @brief 将身份编码为字节并拼接到 XOF
 */
void transcript_absorb_identity(Xof& xof, const Identity& id);

/**
 * @brief 将消息编码为字节并拼接到 XOF
 */
void transcript_absorb_message(Xof& xof, const Message& msg);

} // namespace ibags