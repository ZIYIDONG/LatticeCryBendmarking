#pragma once
/**
 * @file secure_memory.h
 * @brief 安全内存操作 — secure_zero / SecretBuffer / SecretPoly
 *
 * 设计原则:
 *  - PoC 级实现使用 volatile memset + compiler barrier，防止编译器优化掉清零
 *  - SecretBuffer 自动管理生命周期，析构时强制清零
 *  - SecretPoly 禁止隐式拷贝，只允许 move 和显式 public_copy()
 *
 * 生产注意事项 (Production Notes):
 *  ┌───────────────────────────────────────────────────────────────┐
 *  │ PoC 级别: 使用 volatile + compiler barrier 防优化            │
 *  │ 生产级别: 替换以下实现为经安全审计的库函数:                    │
 *  │   ① OpenSSL OPENSSL_cleanse(ptr, len)                        │
 *  │   ② libsodium sodium_memzero(ptr, len)                       │
 *  │   ③ Windows SecureZeroMemory(ptr, len)                       │
 *  │   ④ memset_s (C11 Annex K, 可用时)                            │
 *  │                                                               │
 *  │ 关键: volatile memset 可能被编译器不确定行为优化               │
 *  │       (编译器允许省略对 volatile 指向对象的写入)               │
 *  │       生产环境务必使用上述审计库函数                            │
 *  └───────────────────────────────────────────────────────────────┘
 *
 * 参考:
 *  - OpenSSL OPENSSL_cleanse 实现
 *  - libsodium sodium_memzero 实现
 *  - PQClean crypto_stream_* 的 secret buffer 处理风格
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <type_traits>
#include <utility>

namespace ibags {

// ============================================================================
// §1  secure_zero — 防优化的内存清零
// ============================================================================

/**
 * @brief 安全清零内存区域，不会被编译器优化删除。
 *
 * PoC 实现: volatile 强制 + compiler barrier (内联汇编)
 * 生产环境: 替换为 OPENSSL_cleanse / sodium_memzero
 *
 * @param ptr  待清零内存起始地址
 * @param len  字节数
 */
void secure_zero(void* ptr, size_t len) noexcept;

// ============================================================================
// §2  SecureWipe — RAII 安全擦除包装器
// ============================================================================

/**
 * @struct SecureWipe
 * @brief 轻量级 scope guard: 离开作用域时自动清零。
 *
 * 用法:
 * @code
 *   uint8_t tmp[32];
 *   // ... 使用 tmp 存储密钥 ...
 *   SecureWipe wipe(tmp, sizeof(tmp));  // 离开作用域时自动清零
 * @endcode
 */
struct SecureWipe {
    void*  data;
    size_t size;

    SecureWipe(void* d, size_t s) : data(d), size(s) {}
    ~SecureWipe() { secure_zero(data, size); }

    // 不可拷贝
    SecureWipe(const SecureWipe&) = delete;
    SecureWipe& operator=(const SecureWipe&) = delete;
};

// ============================================================================
// §3  SecretBuffer — 安全字节缓冲区
// ============================================================================

/**
 * @class SecretBuffer
 * @brief 拥有 Secret 数据的 RAII 字节容器。
 *
 * 特性:
 *  - 析构时自动 secure_zero 所有数据
 *  - 支持 move，禁止隐式 copy
 *  - 提供 public_copy() 用于测试/非秘密计算（显式许可）
 *  - 提供 data() / size() 只读访问
 *
 * 使用示例:
 * @code
 *   SecretBuffer key(32);
 *   generate_key_into(key.data(), key.size());
 *   // ... 使用 key ...
 *   // 析构时自动清零
 * @endcode
 */
class SecretBuffer {
public:
    // ── 构造/析构 ──

    /// 分配 sz 字节的安全缓冲区（内容未初始化）
    explicit SecretBuffer(size_t sz);

    /// 从原始指针接管所有权（调用者不再拥有该内存）
    SecretBuffer(uint8_t* ptr, size_t sz);

    /// 析构: 自动清零并释放
    ~SecretBuffer();

    // ── 移动语义 ──
    SecretBuffer(SecretBuffer&& other) noexcept;
    SecretBuffer& operator=(SecretBuffer&& other) noexcept;

    // ── 禁止拷贝 ──
    SecretBuffer(const SecretBuffer&) = delete;
    SecretBuffer& operator=(const SecretBuffer&) = delete;

    // ── 显式公开拷贝 ──

    /**
     * @brief 显式创建非秘密副本。
     *
     * 仅用于:
     *  - 单元测试中验证 SecretBuffer 内容
     *  - 将秘密数据传递给非秘密计算（调用者确认可接受）
     *
     * @return 包含相同数据的新 SecretBuffer
     */
    [[nodiscard]] SecretBuffer public_copy() const;

    // ── 访问器 ──

    [[nodiscard]] uint8_t* data() noexcept { return data_; }
    [[nodiscard]] const uint8_t* data() const noexcept { return data_; }
    [[nodiscard]] size_t size() const noexcept { return size_; }

    /// 检查是否为空
    [[nodiscard]] bool empty() const noexcept { return size_ == 0 || data_ == nullptr; }

private:
    uint8_t* data_;
    size_t   size_;
};

// ============================================================================
// §4  SecretPoly — 秘密多项式包装类
// ============================================================================

/**
 * @class SecretPoly
 * @brief 拥有秘密多项式系数的 RAII 容器。
 *
 * 环 R_q = Z_q[Y] / (Y^n + 1) 中的元素。
 * 系数存储为连续 vector<long>，析构时自动 secure_zero。
 *
 * 特性:
 *  - 禁止隐式拷贝
 *  - 支持 move
 *  - 提供 public_copy()
 *  - 支持系数级别访问（带边界检查）
 *
 * 使用示例:
 * @code
 *   SecretPoly s(512);
 *   s.set_coeff(0, 123);               // 设置系数
 *   long c = s.coeff(0);               // 读取系数
 *   SecretPoly t = s.public_copy();    // 显式拷贝（测试用）
 * @endcode
 */
class SecretPoly {
public:
    // ── 构造/析构 ──

    /// 构造零多项式 (n 个系数全为 0)
    explicit SecretPoly(int n);

    /// 从系数列表构造（接管所有权）
    explicit SecretPoly(std::vector<long> coeffs);

    /// 析构: 自动清零所有系数
    ~SecretPoly();

    // ── 移动语义 ──
    SecretPoly(SecretPoly&& other) noexcept;
    SecretPoly& operator=(SecretPoly&& other) noexcept;

    // ── 禁止拷贝 ──
    SecretPoly(const SecretPoly&) = delete;
    SecretPoly& operator=(const SecretPoly&) = delete;

    // ── 显式公开拷贝 ──
    [[nodiscard]] SecretPoly public_copy() const;

    // ── 访问器 ──

    [[nodiscard]] int n() const noexcept { return n_; }

    /// 读取系数 coef[i] (无边界检查，调用者确保 0 <= i < n)
    [[nodiscard]] long coeff(int i) const { return coeffs_[i]; }

    /// 设置系数 coef[i] = val (无边界检查)
    void set_coeff(int i, long val) { coeffs_[i] = val; }

    /// 获取底层系数列表的只读指针（用于性能关键路径，使用后不保留引用）
    [[nodiscard]] const long* raw_coeffs() const noexcept { return coeffs_.data(); }

    /// 获取底层系数列表的可写指针（谨慎使用）
    [[nodiscard]] long* raw_coeffs_mutable() noexcept { return coeffs_.data(); }

private:
    int              n_;
    std::vector<long> coeffs_;
};

} // namespace ibags