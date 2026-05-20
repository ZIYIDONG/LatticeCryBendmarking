#pragma once
/**
 * @file params.h
 * @brief IBAGS 参数集模块
 *
 * 环: R_p = Z_p[Y] / (Y^m + 1), 默认 m = n = 512
 *
 * 设计原则:
 *  - 参数集为编译期可选的 constexpr 结构体（与 unified_params_plain-LWE.h 风格一致）
 *  - 支持多种参数集，通过 param_id 区分
 *  - 所有输入在构造/加载时必须经过 validate_params() 验证
 *
 * 参考:
 *  - NFLlib 的 ring 参数组织方式 (n, modulus, 多个 modulus 链)
 *  - PQClean 的 parameter header 风格 (单一结构体，命名规范)
 *  - liboqs 的 algorithm/parameter ID 管理风格 (OQS_ALG_* 宏命名)
 *
 * 生产注意事项 (Production Notes):
 *  - PoC: 当前参数集硬编码在源码中
 *  - 生产: 应从受信配置源加载（签名固件 / HSM / 远程 attestation 报告）
 *  - q 的选择: 对 NTT 友好需满足 q ≡ 1 (mod 2n)，当前 q 值为 PoC 暂用值
 */

#include "errors.h"
#include <cstddef>
#include <cstdint>
#include <string>
#include <array>

namespace ibags {

// ============================================================================
// §1  ParamId — 参数集标识符 (liboqs 风格)
// ============================================================================
enum class ParamId : uint16_t {
    /// Demo 级: n=64, q=8191, 与 Plain-LWE Demo 一致，仅测试用
    IBAGS_64_DEMO      = 0x0040,

    /// NIST Level 2: n=512, q=8404993 (24-bit), 对标 Dilithium2
    IBAGS_512_LEVEL2   = 0x0200,

    /// NIST Level 3: n=768, q=16777259 (25-bit), 对标 Dilithium3
    IBAGS_768_LEVEL3   = 0x0300,

    /// NIST Level 5: n=1024, q=4206593 (23-bit), 对标 Dilithium5 / Falcon-1024
    IBAGS_1024_LEVEL5  = 0x0500,

    // ── 保留区间: 0x0600–0x0FFF 用于未来参数集 ──
};

/// 返回参数集 ID 的人类可读名称
[[nodiscard]] const char* param_id_name(ParamId id) noexcept;

// ============================================================================
// §2  Params — 参数集结构体
// ============================================================================

/**
 * @struct Params
 * @brief 存储 IBAGS 方案的所有公有参数。
 *
 * 所有字段在构造后只读。必须通过 default_params_*() 工厂函数创建。
 */
struct Params {
    // ── 标识 ──
    ParamId   param_id;            // 参数集 ID

    // ── 环参数 ──
    int       n;                   // 环维度 / 多项式度数 m (2 的幂)
    int       q;                   // 模数 (NTT-friendly prime 优先)

    // ── 噪声/拒绝采样参数 ──
    double    sigma;               // 离散高斯宽度
    int       eta1;                // 私钥/陷门噪声上界
    int       eta2;                // 签名采样噪声上界
    int       eta_ver;             // 验证噪声界

    // ── 聚合参数 ──
    int       max_signers;         // 最大签名者数量
    int       max_rejection_loops; // 单次签名最大拒绝采样循环数
    int       kappa;               // Fiat-Shamir 挑战多项式非零系数数 (τ)

    // ── 序列化 ──
    int       coefficient_bytes;   // 每个系数的序列化字节数

    // ── 派生字段 ──
    int       poly_bytes;          // 整个多项式的序列化字节数 (n * coefficient_bytes)

    /**
     * @brief 工厂: Demo 级 n=64 参数集
     *
     * q = 8191 (13 位), sigma = 1.0
     * 与 Plain-LWE Demo 级参数一致，仅用于功能测试/CI
     */
    [[nodiscard]] static Params params_demo_64();

    /**
     * @brief 工厂: NIST Level 2 — n=512 参数集
     *
     * q = 8404993 (24 位, NTT-friendly: 8404993 ≡ 1 mod 1024)
     * 经典安全 ≈ 128-bit，对标 Dilithium2
     */
    [[nodiscard]] static Params params_level2_512();

    /**
     * @brief 工厂: NIST Level 3 — n=768 参数集
     *
     * q = 16777259 (25 位, NTT-friendly: 16777259 ≡ 1 mod 1536)
     * 经典安全 ≈ 192-bit，对标 Dilithium3
     */
    [[nodiscard]] static Params params_level3_768();

    /**
     * @brief 工厂: NIST Level 5 — n=1024 参数集
     *
     * q = 4206593 (23 位, NTT-friendly: 4206593 ≡ 1 mod 2048)
     * 经典安全 ≈ 256-bit，对标 Dilithium5 / Falcon-1024
     */
    [[nodiscard]] static Params params_level5_1024();

    // 默认构造（仅内部/测试使用，正常路径请用工厂函数）
    Params() = default;
};

// ============================================================================
// §2b  编译期参数选择（与 Plain-LWE unified_params_plain-LWE.h 一致）
// ============================================================================

/**
 * @brief 返回编译期选定的默认参数集。
 *
 * 通过 CMake -DSECURITY_LEVEL=N 控制：
 *   无参数                       → Demo (n=64,  q=8191)
 *   -DSECURITY_LEVEL=1           → L1   (n=512, q=8404993)
 *   -DSECURITY_LEVEL=3           → L3   (n=768, q=16777259)
 *   -DSECURITY_LEVEL=5           → L5   (n=1024,q=4206593)
 *
 * 与 Plain-LWE 的 unified::default_mp12_params() 模式完全一致。
 */
[[nodiscard]] Params default_params();

// ============================================================================
// §3  参数验证
// ============================================================================

/**
 * @brief 验证参数集的完整性和安全性
 *
 * 检查项:
 *  1. n 是否为 2 的幂且 ≥ MIN_RING_DIM
 *  2. q > 0 且 ≤ MAX_MODULUS
 *  3. sigma > 0
 *  4. 1 ≤ eta_ver ≤ eta2 ≤ eta1
 *  5. max_signers ∈ [1, MAX_SIGNERS]
 *  6. max_rejection_loops > 0
 *  7. kappa > 0 且 kappa ≤ n
 *  8. coefficient_bytes 足够容纳模数
 *
 * @param p  待验证的参数集
 * @return   Ok() 或带描述的无效参数错误
 */
[[nodiscard]] Status validate_params(const Params& p);

// ============================================================================
// §4  参数编码 (Domain Separation)
// ============================================================================

/// 参数编码缓冲区固定大小（足够容纳所有字段）
inline constexpr size_t PARAMS_ENCODE_BYTES = 64;

/**
 * @brief 将参数集规范编码为字节数组，用于哈希的 domain separation
 *
 * 编码格式（小端序）:
 *   Bytes  0– 1: param_id   (uint16_t)
 *   Bytes  2– 3: n          (uint16_t)
 *   Bytes  4– 7: q          (uint32_t)
 *   Bytes  8–15: sigma      (double, IEEE 754)
 *   Bytes 16–17: eta1       (uint16_t)
 *   Bytes 18–19: eta2       (uint16_t)
 *   Bytes 20–21: eta_ver    (uint16_t)
 *   Bytes 22–23: max_signers(uint16_t)
 *   Bytes 24–25: max_rejection_loops (uint16_t)
 *   Bytes 26–27: coefficient_bytes   (uint16_t)
 *   Bytes 28–29: kappa               (uint16_t)
 *   Bytes 30–63: 保留 (填充 0x00)
 *
 * @param p        参数集
 * @param out      输出缓冲区（至少 PARAMS_ENCODE_BYTES 字节）
 * @param out_len  输出缓冲区长度
 * @return         Ok() 或 InvalidEncoding
 */
[[nodiscard]] Status encode_params(const Params& p,
                                   uint8_t* out,
                                   size_t   out_len) noexcept;

} // namespace ibags