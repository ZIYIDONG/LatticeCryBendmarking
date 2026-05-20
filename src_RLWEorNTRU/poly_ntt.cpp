/**
 * @file poly_ntt.cpp
 * @brief NTT / INTT / pointwise multiplication 实现
 *
 * 使用标准 Cooley–Tukey (NTT) 和 Gentleman–Sande (INTT) 算法。
 * 支持模 Montgomery/Barrett 约减。
 *
 * NTT-friendly 素数的要求: q ≡ 1 (mod 2n)
 *   即: 存在 primitive 2n-th root of unity ω ∈ Z_q
 *
 * 环: Z_q[X]/(X^n + 1)，使用反循环卷积 (negacyclic convolution)。
 * 通过预扭/后扭因子将反循环卷积转换为标准循环卷积。
 *
 * 参考:
 *  - Dilithium ntt.c / invntt.c
 *  - Kyber ntt.c
 *  - "Number Theoretic Transform and Its Applications in Lattice-based Cryptography"
 */

#include "../include_RLWEorNTRU/poly_ntt.h"
#include "../include_RLWEorNTRU/poly.h"
#include "../include_RLWEorNTRU/mod_reduce.h"

#include <cassert>
#include <cstring>
#include <algorithm>
#include <cmath>

namespace ibags {

// ============================================================================
// NttPoly 实现
// ============================================================================

NttPoly::NttPoly(int n) : n_(n), coeffs_(n, 0) {
    assert(n > 0);
}

NttPoly::NttPoly(const Params& pp) : n_(pp.n), coeffs_(pp.n, 0) {}

int64_t NttPoly::operator[](int i) const {
    assert(i >= 0 && i < n_);
    return coeffs_[i];
}

int64_t& NttPoly::operator[](int i) {
    assert(i >= 0 && i < n_);
    return coeffs_[i];
}

// ============================================================================
// is_ntt_friendly
// ============================================================================

bool is_ntt_friendly(const Params& pp) {
    int64_t n2 = static_cast<int64_t>(2) * static_cast<int64_t>(pp.n);
    return (pp.q % static_cast<uint32_t>(n2)) == 1;
}

// ============================================================================
// 内部: 计算 primitive 2n-th root of unity in Z_q
// ============================================================================

namespace {

/// 快速幂: a^e mod q
int64_t mod_pow(int64_t a, int64_t e, uint32_t q) {
    int64_t result = 1;
    int64_t base = a % static_cast<int64_t>(q);
    while (e > 0) {
        if (e & 1) {
            result = barrett_reduce(result * base, q);
        }
        base = barrett_reduce(base * base, q);
        e >>= 1;
    }
    return result;
}

/// 寻找 primitive 2n-th root of unity
/// 返回 ψ 满足 ψ^n = -1, ψ^(2n) = 1
Status find_root_of_unity(int n, uint32_t q, int64_t* root) {
    if (!root) {
        return {ErrorCode::InvalidParams, "find_root_of_unity: null output"};
    }

    int64_t n2 = static_cast<int64_t>(2) * n;
    if (q % static_cast<uint32_t>(n2) != 1) {
        return {ErrorCode::InvalidParams,
                "q is not NTT-friendly: q must be ≡ 1 (mod 2n)"};
    }

    int64_t phi = static_cast<int64_t>(q) - 1;
    std::vector<int64_t> factors;

    int64_t tmp = phi;
    for (int64_t p = 2; p * p <= tmp; ++p) {
        if (tmp % p == 0) {
            factors.push_back(p);
            while (tmp % p == 0) tmp /= p;
        }
    }
    if (tmp > 1) factors.push_back(tmp);

    int64_t g = 2;
    while (true) {
        bool is_primitive = true;
        for (int64_t f : factors) {
            if (mod_pow(g, phi / f, q) == 1) {
                is_primitive = false;
                break;
            }
        }
        if (is_primitive) break;
        ++g;
        if (g >= static_cast<int64_t>(q)) {
            return {ErrorCode::InternalError,
                    "Cannot find primitive root mod " + std::to_string(q)};
        }
    }

    // ψ = g^{(q-1)/(2n)} is primitive 2n-th root
    int64_t exponent = phi / n2;
    *root = mod_pow(g, exponent, q);
    return Status::Ok();
}

/// 计算 twiddle factors:
/// zetas[i] = ψ^i (i = 0..2n-1)
/// zetas_inv[i] = ψ^{-i} (i = 0..2n-1)
/// 其中 ψ 是 primitive 2n-th root of unity
Status compute_zetas(int n, const Params& pp,
                     std::vector<int64_t>& zetas,
                     std::vector<int64_t>& zetas_inv) {
    int64_t psi;
    auto status = find_root_of_unity(n, pp.q, &psi);
    if (!status.ok()) return status;

    int m = 2 * n;
    zetas.resize(m);
    zetas_inv.resize(m);

    // ψ^i
    int64_t w = 1;
    for (int i = 0; i < m; ++i) {
        zetas[i] = w;
        w = barrett_reduce(w * psi, pp.q);
    }

    // ψ^{-i} = ψ^{m-i}
    zetas_inv[0] = 1;
    for (int i = 1; i < m; ++i) {
        zetas_inv[i] = zetas[m - i];
    }

    return Status::Ok();
}

/// Bit-reversal permutation (in-place)
void bit_reverse(int64_t* a, int n, int log_n) {
    for (int i = 0; i < n; ++i) {
        int rev = 0;
        for (int j = 0; j < log_n; ++j) {
            if (i & (1 << j)) {
                rev |= (1 << (log_n - 1 - j));
            }
        }
        if (i < rev) {
            std::swap(a[i], a[rev]);
        }
    }
}

} // anonymous namespace

// ============================================================================
// poly_ntt — Cooley–Tukey 正向 NTT (反循环卷积)
// ============================================================================

Status poly_ntt(const Poly& a, const Params& pp, NttPoly* out) {
    if (!out) {
        return {ErrorCode::InvalidParams, "poly_ntt: out is null"};
    }
    if (a.n() != pp.n) {
        return {ErrorCode::InvalidParams, "poly_ntt: dimension mismatch"};
    }
    if (!is_ntt_friendly(pp)) {
        return {ErrorCode::InvalidParams,
                "poly_ntt: q is not NTT-friendly (q must ≡ 1 mod 2n)"};
    }

    int n = pp.n;
    int log_n = 0;
    while ((1 << log_n) < n) ++log_n;
    int64_t q = static_cast<int64_t>(pp.q);

    // 获取 twiddle factors
    std::vector<int64_t> zetas, zetas_inv;
    auto status = compute_zetas(n, pp, zetas, zetas_inv);
    if (!status.ok()) return status;

    // 复制输入并模约减
    std::vector<int64_t> f(n);
    for (int i = 0; i < n; ++i) {
        f[i] = barrett_reduce(a[i], q);
    }

    // 预扭 (pre-twist): f[i] *= ψ^i mod q
    // 将反循环卷积 (mod X^n+1) 转换为标准循环卷积 (mod X^n-1)
    int64_t psi = zetas[1];  // ψ^1
    int64_t pow_psi = 1;
    for (int i = 0; i < n; ++i) {
        f[i] = barrett_reduce(f[i] * pow_psi, q);
        pow_psi = barrett_reduce(pow_psi * psi, q);
    }

    // Bit-reversal permutation
    bit_reverse(f.data(), n, log_n);

    // Cooley–Tukey NTT (DIT), 使用 ω = ψ^2 作为 n-th primitive root
    for (int len = 2; len <= n; len <<= 1) {
        int half = len / 2;
        for (int start = 0; start < n; start += len) {
            for (int j = 0; j < half; ++j) {
                // ω^{j * n/len} = ψ^{2 * j * n/len}
                int zeta_idx = 2 * j * (n / len);
                int64_t zeta = zetas[zeta_idx];

                int64_t t = barrett_reduce(zeta * f[start + half + j], q);
                f[start + half + j] = barrett_reduce(f[start + j] - t, q);
                f[start + j]        = barrett_reduce(f[start + j] + t, q);
            }
        }
    }

    // 输出到 NttPoly
    *out = NttPoly(n);
    for (int i = 0; i < n; ++i) {
        (*out)[i] = f[i];
    }
    return Status::Ok();
}

// ============================================================================
// poly_invntt — Gentleman–Sande 逆 NTT (反循环卷积)
// ============================================================================

Status poly_invntt(const NttPoly& a, const Params& pp, Poly* out) {
    if (!out) {
        return {ErrorCode::InvalidParams, "poly_invntt: out is null"};
    }
    if (a.n() != pp.n) {
        return {ErrorCode::InvalidParams, "poly_invntt: dimension mismatch"};
    }
    if (!is_ntt_friendly(pp)) {
        return {ErrorCode::InvalidParams,
                "poly_invntt: q is not NTT-friendly"};
    }

    int n = pp.n;
    int log_n = 0;
    while ((1 << log_n) < n) ++log_n;
    int64_t q = static_cast<int64_t>(pp.q);

    // 获取 twiddle factors
    std::vector<int64_t> zetas, zetas_inv;
    auto status = compute_zetas(n, pp, zetas, zetas_inv);
    if (!status.ok()) return status;

    // 复制输入
    std::vector<int64_t> f(n);
    for (int i = 0; i < n; ++i) {
        f[i] = barrett_reduce(a[i], q);
    }

    // Gentleman–Sande INTT, 使用 ω^{-1} = ψ^{-2}
    for (int len = n; len >= 2; len >>= 1) {
        int half = len / 2;
        for (int start = 0; start < n; start += len) {
            for (int j = 0; j < half; ++j) {
                int64_t u = f[start + j];
                int64_t v = f[start + half + j];

                f[start + j]        = barrett_reduce(u + v, q);
                f[start + half + j] = barrett_reduce(u - v, q);

                // ψ^{-2 * j * n/len}
                int zeta_inv_idx = 2 * j * (n / len);
                f[start + half + j] = barrett_reduce(
                    zetas_inv[zeta_inv_idx] * f[start + half + j], q);
            }
        }
    }

    // Bit-reversal permutation
    bit_reverse(f.data(), n, log_n);

    // 缩放因子 n^{-1} mod q
    int64_t n_inv = mod_pow(n, q - 2, q);

    // 后扭 (post-twist): f[i] *= ψ^{-i} mod q
    int64_t psi = zetas[1];  // ψ^1
    int64_t psi_inv = mod_pow(psi, q - 2, q);
    int64_t pow_psi_inv = 1;

    // 输出
    *out = Poly(n);
    for (int i = 0; i < n; ++i) {
        (*out)[i] = barrett_reduce(n_inv * f[i] % q * pow_psi_inv, q);
        pow_psi_inv = barrett_reduce(pow_psi_inv * psi_inv, q);
    }

    return Status::Ok();
}

// ============================================================================
// poly_pointwise_mul_ntt
// ============================================================================

Status poly_pointwise_mul_ntt(const NttPoly& a, const NttPoly& b,
                              const Params& pp, NttPoly* out) {
    if (!out) {
        return {ErrorCode::InvalidParams, "poly_pointwise_mul_ntt: out is null"};
    }
    if (a.n() != b.n() || a.n() != pp.n) {
        return {ErrorCode::InvalidParams,
                "poly_pointwise_mul_ntt: dimension mismatch"};
    }

    int n = pp.n;
    *out = NttPoly(n);
    for (int i = 0; i < n; ++i) {
        (*out)[i] = barrett_reduce(a[i] * b[i], pp.q);
    }
    return Status::Ok();
}

// ============================================================================
// ntt_roundtrip_test
// ============================================================================

Status ntt_roundtrip_test(const Poly& a, const Params& pp) {
    if (!is_ntt_friendly(pp)) {
        return {ErrorCode::InvalidParams,
                "ntt_roundtrip_test: q is not NTT-friendly"};
    }

    NttPoly a_ntt(pp);
    auto status = poly_ntt(a, pp, &a_ntt);
    if (!status.ok()) return status;

    Poly a_recovered(pp.n);
    status = poly_invntt(a_ntt, pp, &a_recovered);
    if (!status.ok()) return status;

    Poly a_canon = a;
    a_canon.reduce_mod_q(pp);

    if (!poly_equal_ct(a_canon, a_recovered)) {
        return {ErrorCode::InternalError,
                "ntt_roundtrip_test: INTT(NTT(a)) != a"};
    }
    return Status::Ok();
}

} // namespace ibags