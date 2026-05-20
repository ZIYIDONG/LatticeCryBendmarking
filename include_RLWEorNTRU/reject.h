#pragma once
/**
 * @file reject.h
 * @brief 拒绝采样 — 均匀分布到固定范围内
 *
 * 协议层用法:
 *   - hash_to_uniform: 将 XOF 输出压缩到 [0, q-1] (用于 H1 环哈希)
 *   - coeff_reject_small: 将 XOF 输出压缩到 [-B, B] (用于 H2 挑战)
 *
 * 核心原语:
 *   1. hash_to_uniform_coeff(xof, q): 采样 Z_q 元素 (通过拒绝采样)
 *   2. coeff_reject_small(xof, B):    采样小系数 ∈ [-B, B]
 *
 * 设计原则:
 *   - 确定性 (仅依赖 XOF)
 *   - 偏倚最小 (rejection sampling)
 *   - 恒定时间 (通过条件赋值实现)
 *
 * 参考:
 *   - Dilithium poly_uniform (reject_eta)
 *   - Falcon hash_to_point (拒绝采样到 Z_q)
 *   - PQClean 的 rej_uniform
 */

#include "xof.h"
#include "params.h"
#include "poly.h"
#include "errors.h"

#include <cstdint>
#include <vector>

namespace ibags {

// ============================================================================
// Z_q 拒绝采样
// ============================================================================

/**
 * @brief 从 XOF 采样均匀 Z_q 元素 (coeff-level)
 *
 * 实现:
 *   1. 构造 2^24 范围内的掩码 (24 bits)
 *   2. 从 XOF 读取 uint64_t，取 24 bits
 *   3. 如果 < q，接受；否则拒绝
 *
 * 常数因子:
 *   q ≤ 2^27 → 24-bit 窗口 (NIST L5 q ~ 2^32 → 必要时扩展)
 *
 * 偏差控制: 拒绝采样保证严格的均匀分布
 *
 * @param xof      处于 squeeze 模式的 XOF
 * @param q        模数
 * @param q_shift  log2(q) 上限 (用于掩码)
 * @return         落在 [0, q-1] 的整数
 */
int64_t hash_to_uniform_coeff_reject(Xof& xof, int64_t q, int q_shift);

/**
 * @brief 从 XOF 采样整个多项式 a ∈ Z_q[X]/(X^n+1)
 *
 * 为 n 个系数各调用 hash_to_uniform_coeff_reject 一次
 *
 * @param xof      处于 squeeze 模式的 XOF
 * @param pp       参数集
 * @return         a(x) 多项式 (in Z_q)
 */
Poly hash_to_uniform(Xof& xof, const Params& pp);

// ============================================================================
// 挑战采样: 拒绝到 [-B, B]
// ============================================================================

/**
 * @brief 从 XOF 采样小挑战系数在 [-B, B]
 *
 * 算法:
 *   1. 计算桶数量: ceil(log2(2*B+1))
 *   2. 从 XOF 读取 uint64_t，掩码
 *   3. 如果 c > 2*B，拒绝 (偏倚控制)
 *   4. 返回 c - B
 *
 * 用途: H2 挑战多项式系数生成
 *
 * @param xof         处于 squeeze 模式的 XOF
 * @param B           挑战界 (θ)
 * @param mask_bits   掩码位数
 * @return            系数 ∈ [-B, B]
 */
int64_t coeff_reject_small(Xof& xof, int B, int mask_bits);

} // namespace ibags