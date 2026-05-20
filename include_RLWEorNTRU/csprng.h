#pragma once
/**
 * @file csprng.h
 * @brief CSPRNG 系统 (System + Seed-derived)
 *
 * 设计原则:
 *  - SystemCSPRNG: 使用 OS 熵源 (/dev/urandom, getrandom) 的真随机
 *  - SeedCSPRNG: 从 256-bit seed 确定性派生的伪随机流
 *  - 显式 CSPRNG 接口替代隐式 rand()，确保审计可追踪
 *  - 仅非秘密随机（non-secret randomness）通过此模块提供
 *    （密钥材料使用 secret buffer 模块 + 专用 RNG）
 *
 * 参考：
 *  - NIST SP 800-90A Hash_DRBG / HMAC_DRBG
 *  - PQClean randombytes.c (libsodium/OpenSSL backend)
 *  - Dilithium randombytes 风格
 */

#include "xof.h"

#include <cstdint>
#include <cstddef>
#include <array>
#include <string_view>
#include <vector>

namespace ibags {

// ============================================================================
// 常量
// ============================================================================

/// 种子长度 = 256 位
inline constexpr size_t CSPRNG_SEED_BYTES = 32;

/// 单次最多生成字节数（防止生成器状态泄露）
inline constexpr size_t CSPRNG_MAX_OUTPUT = (1 << 20); // 1 MiB

// ============================================================================
// SystemCSPRNG — OS 熵源真随机
// ============================================================================

/**
 * @class SystemCSPRNG
 * @brief 从 OS 内核熵池获取真随机字节。
 *
 * 后端选择顺序：
 *   1. getrandom() (Linux 3.17+)
 *   2. /dev/urandom
 *   3. OpenSSL RAND_bytes (回退，如果链接)
 *
 * 用法:
 *  @code
 *    SystemCSPRNG rng;
 *    uint8_t seed[32];
 *    rng.randombytes(seed, sizeof(seed));
 *    SeedCSPRNG seeded(seed, sizeof(seed));  // 可重复的测试
 *    // 或
 *    auto seed_vec = rng.randombytes(32);
 *  @endcode
 */
class SystemCSPRNG {
public:
    /// 获取真随机字节
    /// @throws std::runtime_error 如果系统熵源不可用
    void randombytes(uint8_t* out, size_t len);

    /// 便捷: 随机 u64
    uint64_t rand_u64();

    /// 便捷: 范围 [0, bound) 的均匀随机
    uint64_t rand_uint64_bound(uint64_t bound);
};

// ============================================================================
// SeedCSPRNG — 种子确定性伪随机
// ============================================================================

/**
 * @class SeedCSPRNG
 * @brief 从 256-bit seed 确定性派生伪随机流。
 *
 * 实现: SHAKE256 计数器模式 (counter || seed) → output
 * 每个请求从 seed + counter 重新 squeeze。
 *
 * 用法:
 *  @code
 *    SeedCSPRNG csprng(seed, seed_len, domain::GAUSS_FUNCTION);
 *    uint8_t buf[64];
 *    csprng.randombytes(buf, sizeof(buf));
 *  @endcode
 */
class SeedCSPRNG {
public:
    /// 从字节种子构造
    /// @param seed       种子数据 (至少 32 字节)
    /// @param seed_len   种子长度
    /// @param domain_label domain separation label
    SeedCSPRNG(const uint8_t* seed,
               size_t         seed_len,
               std::string_view domain_label);

    /// 从 array 构造
    template<size_t N>
    explicit SeedCSPRNG(const std::array<uint8_t, N>& seed_arr,
                        std::string_view domain_label)
        : SeedCSPRNG(seed_arr.data(), N, domain_label) {}

    /// 从 vector 构造
    explicit SeedCSPRNG(const std::vector<uint8_t>& seed_vec,
                        std::string_view domain_label)
        : SeedCSPRNG(seed_vec.data(), seed_vec.size(), domain_label) {}

    ~SeedCSPRNG() = default;

    // 禁止拷贝
    SeedCSPRNG(const SeedCSPRNG&) = delete;
    SeedCSPRNG& operator=(const SeedCSPRNG&) = delete;

    // 允许移动
    SeedCSPRNG(SeedCSPRNG&&) noexcept = default;
    SeedCSPRNG& operator=(SeedCSPRNG&&) noexcept = default;

    /// 生成伪随机字节
    void randombytes(uint8_t* out, size_t len);

    /// 便捷: 随机 u64
    uint64_t rand_u64();

    /// 便捷: 范围 [0, bound) 的均匀随机 (rejection sampling)
    uint64_t rand_uint64_bound(uint64_t bound);

    /// 当前计数器值 (调试)
    uint64_t counter() const noexcept { return counter_; }

private:
    std::vector<uint8_t> seed_;
    uint64_t             counter_;
    std::string          domain_label_;  // deep-copy 避免 UAF
};

} // namespace ibags