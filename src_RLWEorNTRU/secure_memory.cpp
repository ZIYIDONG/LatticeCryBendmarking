/**
 * @file secure_memory.cpp
 * @brief 安全内存操作实现 — secure_zero / SecretBuffer / SecretPoly
 */

#include "../include_RLWEorNTRU/secure_memory.h"
#include <cstring>
#include <new>
#include <algorithm>

namespace ibags {

// ============================================================================
// §1  secure_zero — 防优化的内存清零
// ============================================================================

// ── PoC 实现 ──
//
// 两层防护:
//   1. 通过 volatile 指针执行 memset，强制编译器生成写入指令
//   2. Compiler barrier (内联 asm) 阻止编译器重排/删除
//
// 局限性 (见 secure_memory.h 生产注意事项):
//   - volatile 写入语义可能被编译器不规则处理
//   - 某些优化器仍可能在同一编译单元内消除它
//   - 不适用于跨函数边界的安全擦除
//
// ── 生产替换 ──
//
// 方案 A (OpenSSL):
//   #include <openssl/crypto.h>
//   OPENSSL_cleanse(ptr, len);
//
// 方案 B (libsodium):
//   #include <sodium.h>
//   sodium_memzero(ptr, len);
//
// 方案 C (手写汇编屏障, MSVC 兼容):
//   volatile unsigned char* p = (volatile unsigned char*)ptr;
//   while (len--) *p++ = 0;
//   // 加上适当的 compiler barrier
//

void secure_zero(void* ptr, size_t len) noexcept {
    if (ptr == nullptr || len == 0) return;

    // volatile 指针确保写入不被优化
    volatile unsigned char* p = static_cast<volatile unsigned char*>(ptr);
    while (len--) {
        *p++ = 0;
    }

    // Compiler barrier: 阻止编译器跨此边界重排或消除
    //
    // GCC/Clang:
    //   __asm__ __volatile__("" : : "r"(ptr) : "memory");
    //
    // MSVC:
    //   _ReadWriteBarrier();
    //
    // 这里使用 GCC/Clang 屏障（Linux 目标）
#if defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__("" : : "r"(ptr) : "memory");
#endif
}

// ============================================================================
// §2  SecretBuffer 实现
// ============================================================================

SecretBuffer::SecretBuffer(size_t sz)
    : data_(nullptr), size_(sz)
{
    if (sz > 0) {
        data_ = new (std::nothrow) uint8_t[sz];
        if (data_ == nullptr) {
            size_ = 0;
            throw std::bad_alloc();
        }
        // 初始化为零（防御性编程：不假设调用者会填充所有字节）
        std::memset(data_, 0, sz);
    }
}

SecretBuffer::SecretBuffer(uint8_t* ptr, size_t sz)
    : data_(ptr), size_(sz)
{
    // 接管所有权: 传入的 ptr 现在由本对象管理
    // 如果 ptr 为空，size_ 必须为 0
}

SecretBuffer::~SecretBuffer() {
    if (data_ != nullptr && size_ > 0) {
        secure_zero(data_, size_);
        delete[] data_;
    }
    data_ = nullptr;
    size_ = 0;
}

// ── 移动构造 ──
SecretBuffer::SecretBuffer(SecretBuffer&& other) noexcept
    : data_(other.data_), size_(other.size_)
{
    // 清空源对象，防止 double-free
    other.data_ = nullptr;
    other.size_ = 0;
}

// ── 移动赋值 ──
SecretBuffer& SecretBuffer::operator=(SecretBuffer&& other) noexcept {
    if (this != &other) {
        // 先清理自身
        if (data_ != nullptr && size_ > 0) {
            secure_zero(data_, size_);
            delete[] data_;
        }
        // 接管
        data_ = other.data_;
        size_ = other.size_;
        other.data_ = nullptr;
        other.size_ = 0;
    }
    return *this;
}

// ── 显式公开拷贝 ──
SecretBuffer SecretBuffer::public_copy() const {
    SecretBuffer copy(size_);
    if (data_ != nullptr && size_ > 0) {
        std::memcpy(copy.data_, data_, size_);
    }
    return copy;
}

// ============================================================================
// §3  SecretPoly 实现
// ============================================================================

SecretPoly::SecretPoly(int n)
    : n_(n), coeffs_(static_cast<size_t>(n), 0L)
{
    // 初始化为零多项式
}

SecretPoly::SecretPoly(std::vector<long> coeffs)
    : n_(static_cast<int>(coeffs.size())), coeffs_(std::move(coeffs))
{
}

SecretPoly::~SecretPoly() {
    if (!coeffs_.empty()) {
        // 对系数内存执行 secure_zero
        // vector<long> 底层是连续内存: reinterpret as bytes
        secure_zero(coeffs_.data(), coeffs_.size() * sizeof(long));
    }
    n_ = 0;
}

// ── 移动构造 ──
SecretPoly::SecretPoly(SecretPoly&& other) noexcept
    : n_(other.n_), coeffs_(std::move(other.coeffs_))
{
    other.n_ = 0;
    // other.coeffs_ 已被 move 置空
}

// ── 移动赋值 ──
SecretPoly& SecretPoly::operator=(SecretPoly&& other) noexcept {
    if (this != &other) {
        // 清理自身
        if (!coeffs_.empty()) {
            secure_zero(coeffs_.data(), coeffs_.size() * sizeof(long));
        }
        // 接管
        n_ = other.n_;
        coeffs_ = std::move(other.coeffs_);
        other.n_ = 0;
    }
    return *this;
}

// ── 显式公开拷贝 ──
SecretPoly SecretPoly::public_copy() const {
    return SecretPoly(coeffs_);  // 拷贝构造 vector<long>
}

} // namespace ibags