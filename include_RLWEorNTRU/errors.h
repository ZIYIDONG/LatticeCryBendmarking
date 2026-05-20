#pragma once
/**
 * @file errors.h
 * @brief 统一错误处理模块 — IBAGS (Identity-Based Aggregate Signature)
 *
 * 设计原则:
 *  - 不使用异常作为主流程错误处理（密码库通常禁用异常）
 *  - Status 类轻量级，支持移动语义，兼容 C++20
 *  - ErrorCode 覆盖所有基础操作可能的错误路径
 *
 * 参考风格:
 *  - Abseil absl::Status / absl::StatusCode
 *  - liboqs OQS_STATUS 返回码
 *
 * 生产注意事项 (Production Notes):
 *  - 当前为 PoC 级实现，使用 std::string 存储消息
 *  - 生产环境建议:
 *    ① 使用固定大小错误消息缓冲区（避免堆分配）
 *    ② 将常见消息编译为静态 const char* 常量
 *    ③ 考虑使用压缩错误码枚举（uint16_t）以减小 ABI 表面积
 */

#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

namespace ibags {

// ============================================================================
// §1  ErrorCode — 覆盖所有基础操作错误路径
// ============================================================================
enum class ErrorCode : int32_t {
    // ── 成功 ──
    Ok = 0,

    // ── 参数/输入验证 (1–9) ──
    InvalidParams        = 1,   // 参数非法（n 非 2 的幂 / q ≤ 0 / eta 不合理）
    InvalidEncoding      = 2,   // 编码/解码错误
    InvalidSignerCount   = 3,   // 签名者数量超限或为零
    DuplicateSigner      = 4,   // 签名者列表中检测到重复身份

    // ── 数学对象验证 (10–19) ──
    InvalidPolynomial    = 10,  // 多项式非法（系数超范围 / 度数错误）
    InvalidPublicKey     = 11,  // 公钥非法（不在环上 / 格式错误）
    InvalidSecretKey     = 12,  // 密钥非法
    InvalidSignature     = 13,  // 签名非法（格式 / 范数 / 验证失败）

    // ── 范数/拒绝采样 (20–24) ──
    NormTooLarge         = 20,  // 范数超出阈值（拒绝采样失败）

    // ── 陷门 / 采样 (30–39) ──
    TrapdoorError        = 30,  // NTRU TrapGen 失败
    SamplingError        = 31,  // 高斯采样失败（拒绝循环超限）

    // ── 内部错误 (90–99) ──
    InternalError        = 90,  // 内部一致性错误（不应发生）
};

// ============================================================================
// §2  Status — 轻量级错误状态
// ============================================================================

/**
 * @class Status
 * @brief 表示操作结果：Ok 或带错误码的失败状态。
 *
 * 用法:
 * @code
 *   Status validate(const Params& p) {
 *     if (p.n <= 0) return {ErrorCode::InvalidParams, "n must be positive"};
 *     return Status::Ok();
 *   }
 *
 *   Status s = some_operation();
 *   if (!s.ok()) { log_error(s.message()); return s; }
 * @endcode
 */
class Status {
public:
    // ── 构造 ──

    /// 默认构造 = Ok
    Status() : code_(ErrorCode::Ok), msg_() {}

    /// 仅错误码（无附加消息）
    explicit Status(ErrorCode c) : code_(c), msg_() {}

    /// 错误码 + 消息
    Status(ErrorCode c, std::string msg) : code_(c), msg_(std::move(msg)) {}

    /// 工厂：快速创建 Ok 状态
    static Status Ok() { return Status(); }

    // ── 查询 ──

    [[nodiscard]] bool ok() const noexcept { return code_ == ErrorCode::Ok; }

    [[nodiscard]] ErrorCode code() const noexcept { return code_; }

    [[nodiscard]] const std::string& message() const noexcept { return msg_; }

    // ── 布尔转换 ──
    explicit operator bool() const noexcept { return ok(); }
    bool operator!() const noexcept { return !ok(); }

    // ── 比较 ──
    bool operator==(const Status& other) const noexcept {
        return code_ == other.code_;
    }
    bool operator!=(const Status& other) const noexcept {
        return !(*this == other);
    }

private:
    ErrorCode    code_;
    std::string  msg_;
};

/**
 * @brief 便捷宏：提前返回错误
 *
 * 使用示例:
 * @code
 *   IBAGS_RETURN_IF_ERROR(validate_params(p));
 * @endcode
 */
#define IBAGS_RETURN_IF_ERROR(expr)          \
    do {                                      \
        ::ibags::Status _s = (expr);          \
        if (!_s.ok()) return _s;              \
    } while (0)

/**
 * @brief 便捷宏：断言必须 Ok，否则返回 InternalError
 */
#define IBAGS_ASSIGN_OR_RETURN(lhs, expr)                    \
    auto _tmp_status_##__LINE__ = (expr);                     \
    if (!_tmp_status_##__LINE__.ok())                         \
        return ::ibags::Status(::ibags::ErrorCode::InternalError, \
                               "_tmp_status_##__LINE__ failed"); \
    lhs = std::move(_tmp_status_##__LINE__##_value)

} // namespace ibags