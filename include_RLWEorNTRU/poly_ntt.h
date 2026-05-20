#pragma once
/**
 * @file poly_ntt.h
 * @brief NTT-domain 多项式表示及 NTT 变换
 *
 * 设计原则:
 *  - NttPoly 是 NTT 域的独立类型，不与普通 Poly 混用
 *  - 仅对 NTT-friendly 的 q 有效 (q ≡ 1 mod 2n)
 *  - 支持正向 NTT、逆向 INTT、逐点乘法
 *
 * 参考:
 *  - NFLlib NTT
 *  - Dilithium ntt.c / invntt.c
 *  - PQClean ML-DSA NTT 实现
 */

#include "params.h"
#include "errors.h"

#include <cstdint>
#include <vector>

namespace ibags {

// ============================================================================
// NttPoly — NTT 域多项式
// ============================================================================

class NttPoly {
public:
    explicit NttPoly(int n);
    explicit NttPoly(const Params& pp);

    int n() const noexcept { return n_; }
    const int64_t* data() const noexcept { return coeffs_.data(); }
    int64_t* data() noexcept { return coeffs_.data(); }
    int64_t operator[](int i) const;
    int64_t& operator[](int i);

    NttPoly(const NttPoly&) = default;
    NttPoly& operator=(const NttPoly&) = default;
    NttPoly(NttPoly&&) noexcept = default;
    NttPoly& operator=(NttPoly&&) noexcept = default;

private:
    int n_;
    std::vector<int64_t> coeffs_;
};

// ============================================================================
// NTT / INTT 接口
// ============================================================================

/// 正向 NTT: normal domain → NTT domain
/// 要求 q 为 NTT-friendly 素数 (q ≡ 1 mod 2n)
Status poly_ntt(const class Poly& a, const Params& pp, NttPoly* out);

/// 逆 NTT: NTT domain → normal domain (canonical)
Status poly_invntt(const NttPoly& a, const Params& pp, class Poly* out);

/// NTT 域逐点乘法
/// 输入必须是 NTT-domain type，输出仍为 NTT-domain
Status poly_pointwise_mul_ntt(const NttPoly& a, const NttPoly& b,
                              const Params& pp, NttPoly* out);

/// 判断 q 是否满足 NTT-friendly 条件
bool is_ntt_friendly(const Params& pp);

/// 检查 NTT 实现的 roundtrip 一致性
/// INTT(NTT(a)) == a (带 tolerance)
Status ntt_roundtrip_test(const class Poly& a, const Params& pp);

} // namespace ibags