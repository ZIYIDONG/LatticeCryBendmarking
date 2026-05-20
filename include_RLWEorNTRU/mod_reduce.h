#pragma once
/**
 * @file mod_reduce.h
 * @brief Barrett/Montgomery 模约减与 centered representation
 *
 * 参考:
 *  - Dilithium reduce.c (Barrett reduction)
 *  - Kyber/ML-KEM reduce.c (Montgomery reduction)
 */

#include <cstdint>

namespace ibags {

// ============================================================================
// Barrett Reduction
// ============================================================================

/// Barrett 模约减：将 int64 规约到 [0, q)
/// 要求 q < 2^30
int64_t barrett_reduce(int64_t a, uint32_t q);

/// Constant-time Barrett 模约减
int64_t barrett_reduce_ct(int64_t a, uint32_t q);

// ============================================================================
// Centered Representation
// ============================================================================

/// 将 [0, q) 映射到 (−q/2, q/2]
int64_t to_centered(int64_t a, uint32_t q);

/// 将 (−q/2, q/2] 映射回 [0, q)
int64_t from_centered(int64_t a, uint32_t q);

} // namespace ibags