/**
 * @file poly.cpp
 * @brief R_p 多项式实现
 */

#include "../include_RLWEorNTRU/poly.h"
#include "../include_RLWEorNTRU/secure_memory.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace ibags {

// ============================================================================
// Poly 类实现
// ============================================================================

Poly::Poly(int n) : n_(n), coeffs_(n, 0) {
    assert(n > 0 && "n must be positive");
}

Poly::Poly(int n, std::vector<int64_t> coeffs)
    : n_(n), coeffs_(std::move(coeffs)) {
    if (static_cast<int>(coeffs_.size()) != n) {
        throw std::invalid_argument(
            "Poly: coeff vector size " + std::to_string(coeffs_.size()) +
            " != n " + std::to_string(n));
    }
}

Poly Poly::zero(const Params& pp) {
    return Poly(pp.n);
}

Poly Poly::from_coeffs(const Params& pp, std::vector<int64_t> coeffs) {
    return Poly(pp.n, std::move(coeffs));
}

Poly Poly::from_canonical(const Params& pp, std::vector<int64_t> coeffs) {
    Poly p(pp.n, std::move(coeffs));
    p.reduce_mod_q(pp);
    return p;
}

int64_t Poly::operator[](int i) const {
    assert(i >= 0 && i < n_ && "index out of range");
    return coeffs_[i];
}

int64_t& Poly::operator[](int i) {
    assert(i >= 0 && i < n_ && "index out of range");
    return coeffs_[i];
}

bool Poly::is_canonical(const Params& pp) const {
    if (n_ != pp.n) return false;
    int64_t q = static_cast<int64_t>(pp.q);
    for (int i = 0; i < n_; ++i) {
        if (coeffs_[i] < 0 || coeffs_[i] >= q) return false;
    }
    return true;
}

bool Poly::is_centered(const Params& pp) const {
    if (n_ != pp.n) return false;
    int64_t half_q = static_cast<int64_t>(pp.q) / 2;
    for (int i = 0; i < n_; ++i) {
        if (coeffs_[i] < -half_q || coeffs_[i] > half_q) return false;
    }
    return true;
}

void Poly::reduce_mod_q(const Params& pp) {
    assert(n_ == pp.n && "dimension mismatch");
    for (int i = 0; i < n_; ++i) {
        coeffs_[i] = barrett_reduce(coeffs_[i], pp.q);
    }
}

Poly Poly::normalize_centered(const Params& pp) const {
    Poly result = *this;
    for (int i = 0; i < n_; ++i) {
        result.coeffs_[i] = to_centered(result.coeffs_[i], pp.q);
    }
    return result;
}

// ============================================================================
// 自由函数实现
// ============================================================================

void poly_reduce(Poly& poly, const Params& pp) {
    poly.reduce_mod_q(pp);
}

Poly poly_ring_reduce_raw(const Params& pp, const std::vector<int64_t>& conv) {
    // conv 长度为 2n-1，使用 Y^n = -1 折叠
    size_t len = conv.size();
    size_t n = static_cast<size_t>(pp.n);
    Poly result(pp.n);

    for (size_t i = 0; i < n; ++i) {
        int64_t sum = 0;
        if (i < len) sum = conv[i];
        // 对于 Y^{n + i} 项: 系数乘 -1 加到位置 i
        if (i + n < len) {
            sum -= conv[i + n]; // Y^n = -1, 所以 Y^{n+i} -> -Y^i
        }
        result[i] = sum;
    }

    // 将结果规约到 [0, q)
    result.reduce_mod_q(pp);
    return result;
}

// ──── 算术运算 ────

Poly poly_add(const Poly& a, const Poly& b, const Params& pp) {
    assert(a.n() == b.n() && a.n() == pp.n);
    Poly result(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        result[i] = barrett_reduce(a[i] + b[i], pp.q);
    }
    return result;
}

Poly poly_sub(const Poly& a, const Poly& b, const Params& pp) {
    assert(a.n() == b.n() && a.n() == pp.n);
    Poly result(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        result[i] = barrett_reduce(a[i] - b[i], pp.q);
    }
    return result;
}

Poly poly_neg(const Poly& a, const Params& pp) {
    assert(a.n() == pp.n);
    Poly result(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        result[i] = barrett_reduce(-a[i], pp.q);
    }
    return result;
}

Poly poly_mul_naive(const Poly& a, const Poly& b, const Params& pp) {
    assert(a.n() == b.n() && a.n() == pp.n);
    int n = pp.n;

    // 卷积结果长度 2n-1
    std::vector<int64_t> conv(2 * n - 1, 0);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            conv[i + j] += a[i] * b[j];
        }
    }

    return poly_ring_reduce_raw(pp, conv);
}

// ──── 范数 ────

uint64_t poly_norm_inf(const Poly& a, const Params& pp) {
    assert(a.n() == pp.n);
    int64_t half_q = static_cast<int64_t>(pp.q) / 2;
    uint64_t max_abs = 0;

    for (int i = 0; i < pp.n; ++i) {
        // 先转换到 centered 表示
        int64_t c = barrett_reduce(a[i], pp.q);
        if (c > half_q) c -= static_cast<int64_t>(pp.q);

        // 安全计算绝对值，防止 INT64_MIN 溢出
        uint64_t abs_val;
        if (c == INT64_MIN) {
            abs_val = static_cast<uint64_t>(INT64_MAX) + 1;
        } else if (c < 0) {
            abs_val = static_cast<uint64_t>(-c);
        } else {
            abs_val = static_cast<uint64_t>(c);
        }
        if (abs_val > max_abs) max_abs = abs_val;
    }
    return max_abs;
}

bool poly_norm_bound_check(const Poly& a, uint64_t bound, const Params& pp) {
    assert(a.n() == pp.n);
    int64_t half_q = static_cast<int64_t>(pp.q) / 2;
    bool result = true; // 默认为 true，不早停

    for (int i = 0; i < pp.n; ++i) {
        int64_t c = barrett_reduce(a[i], pp.q);
        if (c > half_q) c -= static_cast<int64_t>(pp.q);

        uint64_t abs_val;
        if (c == INT64_MIN) {
            abs_val = static_cast<uint64_t>(INT64_MAX) + 1;
        } else if (c < 0) {
            abs_val = static_cast<uint64_t>(-c);
        } else {
            abs_val = static_cast<uint64_t>(c);
        }

        // constant-time: 不早停，用逻辑 AND
        bool within = (abs_val <= bound);
        result = result && within;
    }
    return result;
}

bool poly_equal_ct(const Poly& a, const Poly& b) {
    if (a.n() != b.n()) return false;

    // constant-time 比较: XOR 累积差异
    int64_t diff = 0;
    for (int i = 0; i < a.n(); ++i) {
        diff |= a[i] ^ b[i];
    }
    return diff == 0;
}

// ============================================================================
// poly_inv — 扩展欧几里得算法 (NTRU 风格)
// ============================================================================

namespace {

/// 整数模逆: a^{-1} mod m, 若不可逆返回 -1
[[nodiscard]] int64_t mod_inv_int64(int64_t a, int64_t m) {
    if (m <= 0) return -1;
    a = ((a % m) + m) % m;
    int64_t t = 0, newt = 1;
    int64_t r = m, newr = a;
    while (newr != 0) {
        int64_t quotient = r / newr;
        int64_t tmp = newt;
        newt = t - quotient * newt;
        t = tmp;
        tmp = newr;
        newr = r - quotient * newr;
        r = tmp;
    }
    if (r > 1) return -1;
    if (t < 0) t += m;
    return t % m;
}

/// 多项式度数（最高非零系数下标），零多项式返回 -1
[[nodiscard]] int poly_deg(const std::vector<int64_t>& a) {
    for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) {
        if (a[i] != 0) return i;
    }
    return -1;
}

/// 规约多项式系数到 [0, q)
void poly_mod_q(std::vector<int64_t>& a, int64_t q) {
    for (auto& c : a) {
        c = ((c % q) + q) % q;
    }
}

/// 删除多项式尾部零系数
void poly_trim(std::vector<int64_t>& a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}

/// 多项式乘以标量 mod q
void poly_scale_mod_q(std::vector<int64_t>& a, int64_t s, int64_t q) {
    s = ((s % q) + q) % q;
    for (auto& c : a) {
        c = (c * s) % q;
        if (c < 0) c += q;
    }
}

/// 多项式加法: a += b (mod q), a 自动扩容
void poly_add_raw(std::vector<int64_t>& a, const std::vector<int64_t>& b, int64_t q) {
    if (b.size() > a.size()) a.resize(b.size(), 0);
    for (size_t i = 0; i < b.size(); ++i) {
        a[i] = (a[i] + b[i]) % q;
        if (a[i] < 0) a[i] += q;
    }
}

/// 多项式减法: a -= b (mod q), a 自动扩容
void poly_sub_raw(std::vector<int64_t>& a, const std::vector<int64_t>& b, int64_t q) {
    if (b.size() > a.size()) a.resize(b.size(), 0);
    for (size_t i = 0; i < b.size(); ++i) {
        a[i] = (a[i] - b[i]) % q;
        if (a[i] < 0) a[i] += q;
    }
}

/// 多项式乘法（朴素卷积），结果长度 len(a)+len(b)-1
[[nodiscard]] std::vector<int64_t> poly_mul_raw(
        const std::vector<int64_t>& a,
        const std::vector<int64_t>& b,
        int64_t q) {
    if (a.empty() || b.empty()) return {};
    std::vector<int64_t> res(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            if (b[j] == 0) continue;
            int64_t prod = (a[i] * b[j]) % q;
            if (prod < 0) prod += q;
            res[i + j] = (res[i + j] + prod) % q;
            if (res[i + j] < 0) res[i + j] += q;
        }
    }
    poly_trim(res);
    return res;
}

/// 环约减: 将长度 ≤ 2n-1 的多项式对 x^n + 1 取模
[[nodiscard]] std::vector<int64_t> poly_ring_reduce_vec(
        const std::vector<int64_t>& conv, int n, int64_t q) {
    std::vector<int64_t> res(n, 0);
    size_t len = conv.size();
    for (int i = 0; i < n; ++i) {
        int64_t sum = 0;
        if (static_cast<size_t>(i) < len) sum = conv[i];
        // x^{n+i} = -x^i, 所以系数乘 -1
        int j = i + n;
        if (static_cast<size_t>(j) < len) {
            sum = (sum - conv[j]) % q;
        }
        res[i] = (sum % q + q) % q;
    }
    return res;
}

} // anonymous namespace

Status poly_inv(const Poly& a, const Params& pp, Poly* inv) {
    if (!inv) {
        return {ErrorCode::InvalidParams, "poly_inv: inv is null"};
    }
    if (a.n() != pp.n) {
        return {ErrorCode::InvalidParams, "poly_inv: dimension mismatch"};
    }

    int n = pp.n;
    int64_t q = static_cast<int64_t>(pp.q);

    // ── 扩展欧几里得算法在 Z_q[x] 中计算 a(x)^{-1} mod (q, x^n+1) ──
    // 参考: IEEE P1363.1 / NTRU inversion / Falcon keygen poly_small_invert
    //
    // 算法:
    //   r0 ← x^n + 1 (模多项式)
    //   r1 ← a (mod q)
    //   t0 ← 0, t1 ← 1
    //   while r1 ≠ 0:
    //     (quotient, remainder) = poly_div(r0, r1)
    //     t = t0 - quotient * t1
    //     r0 ← r1, r1 ← remainder
    //     t0 ← t1, t1 ← t
    //   若 r0 为非零常数 λ，则 inv = λ^{-1} * t0 mod (q, x^n+1)
    //   否则，a 不可逆

    // 规约输入系数到 [0, q)
    std::vector<int64_t> r0(n + 1, 0);
    r0[0] = 1; r0[n] = 1;          // r0 = x^n + 1
    std::vector<int64_t> r1 = a.coeffs();
    poly_mod_q(r1, q);
    poly_trim(r1);

    // 检查零多项式
    if (r1.empty()) {
        return {ErrorCode::InvalidPolynomial,
                "poly_inv: zero polynomial is not invertible"};
    }

    // t0 = 0, t1 = 1
    std::vector<int64_t> t0;
    std::vector<int64_t> t1 = {1};  // t1 = 1

    const int max_iter = 10 * n; // 安全上限

    for (int iter = 0; iter < max_iter; ++iter) {
        if (r1.empty()) break; // r1 = 0, 算法终止

        // ── 多项式除法: r0 / r1 ──
        // quotient, remainder 使得 r0 = quotient * r1 + remainder, deg(remainder) < deg(r1)
        std::vector<int64_t> remainder = r0;
        std::vector<int64_t> quotient;
        int deg_r1 = poly_deg(r1);
        int64_t lead_r1 = r1[static_cast<size_t>(deg_r1)];
        int64_t inv_lead_r1 = mod_inv_int64(lead_r1, q);
        if (inv_lead_r1 < 0) {
            return {ErrorCode::InvalidPolynomial,
                    "poly_inv: leading coefficient of r1 not invertible mod q"};
        }

        while (true) {
            int deg_rem = poly_deg(remainder);
            if (deg_rem < deg_r1) break;
            int deg_diff = deg_rem - deg_r1;
            int64_t factor = (remainder[static_cast<size_t>(deg_rem)] * inv_lead_r1) % q;
            if (factor < 0) factor += q;

            // 记录商项
            if (static_cast<int>(quotient.size()) <= deg_diff) {
                quotient.resize(static_cast<size_t>(deg_diff) + 1, 0);
            }
            quotient[static_cast<size_t>(deg_diff)] = factor;

            // remainder -= factor * x^{deg_diff} * r1
            for (int k = 0; k <= deg_r1; ++k) {
                size_t idx = static_cast<size_t>(k + deg_diff);
                int64_t sub = (factor * r1[static_cast<size_t>(k)]) % q;
                if (sub < 0) sub += q;
                remainder[idx] = (remainder[idx] - sub) % q;
                if (remainder[idx] < 0) remainder[idx] += q;
            }
            poly_trim(remainder);
        }

        // ── 更新 t: t = t0 - quotient * t1 ──
        std::vector<int64_t> qt1 = poly_mul_raw(quotient, t1, q);
        std::vector<int64_t> t_new = t0;
        poly_sub_raw(t_new, qt1, q);
        poly_trim(t_new);
        // 环约减 t_new 到度数 < n
        if (static_cast<int>(t_new.size()) > n) {
            t_new = poly_ring_reduce_vec(t_new, n, q);
        }
        t_new.resize(static_cast<size_t>(n), 0);
        poly_mod_q(t_new, q);

        // ── 移位: r0 ← r1, r1 ← remainder, t0 ← t1, t1 ← t_new ──
        r0 = std::move(r1);
        r1 = std::move(remainder);
        t0 = std::move(t1);
        t1 = std::move(t_new);
    }

    // 检查 r0 是否为非零常数
    int deg_r0 = poly_deg(r0);
    if (deg_r0 != 0) {
        return {ErrorCode::InvalidPolynomial,
                "poly_inv: polynomial is not invertible in Z_" +
                std::to_string(q) + "[x]/(x^" + std::to_string(n) + "+1)"};
    }

    int64_t lambda = r0[0];
    int64_t inv_lambda = mod_inv_int64(lambda, q);
    if (inv_lambda < 0) {
        return {ErrorCode::InvalidPolynomial,
                "poly_inv: gcd(" + std::to_string(lambda) + ", " +
                std::to_string(q) + ") != 1, not invertible"};
    }

    // inv = lambda^{-1} * t0 mod (q, x^n+1)
    poly_scale_mod_q(t0, inv_lambda, q);
    t0.resize(static_cast<size_t>(n), 0);

    *inv = Poly(pp.n, std::move(t0));
    inv->reduce_mod_q(pp);
    return Status::Ok();
}

// ============================================================================
// poly_is_invertible_mod_q — 仅检测是否可逆（不计算逆）
// ============================================================================

bool poly_is_invertible_mod_q(const Poly& a, const Params& pp) {
    Poly dummy(pp.n);
    Status st = poly_inv(a, pp, &dummy);
    if (!st.ok()) {
        // 安全擦除临时变量
        secure_zero(dummy.data(), static_cast<size_t>(pp.n) * sizeof(int64_t));
    }
    return st.ok();
}

} // namespace ibags
