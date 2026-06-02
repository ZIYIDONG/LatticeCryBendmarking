/**
 * @file poly_ntt.cpp
 * @brief NTT / INTT / pointwise multiplication 实现
 *
 * 使用标准 Cooley–Tukey (NTT) 和 Gentleman–Sande (INTT) 算法。
 * 环: Z_q[X]/(X^n + 1)，使用反循环卷积 (negacyclic convolution)。
 * 通过预扭/后扭因子将反循环卷积转换为标准循环卷积。
 *
 * NttTable 设计:
 *  所有预计算数据 (twiddle factors, bitrev, 约减常数) 在 NttTable::create()
 *  中一次性计算并缓存。NTT/INTT 接口直接使用预计算表，零运行时初始化开销。
 *
 * NTT-friendly 素数的要求: q ≡ 1 (mod 2n)
 *  即: 存在 primitive 2n-th root of unity ψ ∈ Z_q, 满足 ψ^n = -1, ψ^{2n} = 1
 *
 * 参考:
 *  - Dilithium ntt.c / invntt.c
 *  - Kyber ntt.c
 *  - "Number Theoretic Transform and Its Applications in Lattice-based Cryptography"
 *    (Longa & Naehrig, 2016)
 */

#include "../include_RLWEorNTRU/poly_ntt.h"
#include "../include_RLWEorNTRU/poly.h"

#include <cassert>
#include <cstring>
#include <algorithm>
#include <stdexcept>

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
// 内部: 寻找 primitive 2n-th root of unity (用于 NttTable::create)
// ============================================================================

namespace {

/// 快速幂: a^e mod q (使用 Barrett 约减)
int64_t mod_pow(int64_t a, int64_t e, const BarrettConst& bc) {
    int64_t result = 1;
    int64_t base = a;
    while (e > 0) {
        if (e & 1) {
            result = barrett_reduce(result * base, bc);
        }
        base = barrett_reduce(base * base, bc);
        e >>= 1;
    }
    return result;
}

/// 寻找 primitive 2n-th root of unity ψ
/// 返回 ψ 满足 ψ^n = -1, ψ^(2n) = 1
int64_t find_root_of_unity(int n, const BarrettConst& bc) {
    const int n2 = 2 * n;
    const int64_t phi = static_cast<int64_t>(bc.q) - 1;

    // 分解 phi = q - 1 的素因子
    std::vector<int64_t> factors;
    int64_t tmp = phi;
    for (int64_t p = 2; p * p <= tmp; ++p) {
        if (tmp % p == 0) {
            factors.push_back(p);
            while (tmp % p == 0) tmp /= p;
        }
    }
    if (tmp > 1) factors.push_back(tmp);

    // 寻找生成元 g ∈ Z_q^*
    int64_t g = 2;
    while (true) {
        bool is_primitive = true;
        for (int64_t f : factors) {
            if (mod_pow(g, phi / f, bc) == 1) {
                is_primitive = false;
                break;
            }
        }
        if (is_primitive) break;
        ++g;
        if (g >= static_cast<int64_t>(bc.q)) {
            throw std::invalid_argument(
                "Cannot find primitive root mod " + std::to_string(bc.q));
        }
    }

    // ψ = g^{(q-1)/(2n)} is primitive 2n-th root
    int64_t exponent = phi / n2;
    return mod_pow(g, exponent, bc);
}

} // anonymous namespace

// ============================================================================
// NttTable::create — 预计算工厂
// ============================================================================

NttTable NttTable::create(const Params& pp) {
    if (!is_ntt_friendly(pp)) {
        throw std::invalid_argument(
            "NttTable::create: q=" + std::to_string(pp.q) +
            " is not NTT-friendly for n=" + std::to_string(pp.n) +
            " (need q ≡ 1 mod " + std::to_string(2 * pp.n) + ")");
    }

    NttTable tbl;
    tbl.n = pp.n;
    tbl.q = pp.q;

    // log_n
    tbl.log_n = 0;
    while ((1 << tbl.log_n) < pp.n) ++tbl.log_n;

    // 约减常数
    tbl.barrett = make_barrett(pp.q);
    tbl.montgomery = make_montgomery(pp.q);

    // ψ = primitive 2n-th root of unity
    const int64_t psi = find_root_of_unity(pp.n, tbl.barrett);

    // zetas[i] = ψ^i  for i = 0..2n-1
    const int m = 2 * pp.n;
    tbl.zetas.resize(m);
    tbl.zetas_inv.resize(m);

    int64_t w = 1;
    for (int i = 0; i < m; ++i) {
        tbl.zetas[i] = w;
        w = barrett_reduce(w * psi, tbl.barrett);
    }

    // zetas_inv[i] = ψ^{-i} = ψ^{m-i}
    tbl.zetas_inv[0] = 1;
    for (int i = 1; i < m; ++i) {
        tbl.zetas_inv[i] = tbl.zetas[m - i];
    }

    // Bit-reversal permutation table
    tbl.bitrev.resize(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        int rev = 0;
        for (int j = 0; j < tbl.log_n; ++j) {
            if (i & (1 << j)) {
                rev |= (1 << (tbl.log_n - 1 - j));
            }
        }
        tbl.bitrev[i] = rev;
    }

    // n^{-1} mod q (用于 INTT 最终缩放)
    tbl.n_inv = mod_pow(static_cast<int64_t>(pp.n), pp.q - 2, tbl.barrett);

    return tbl;
}

// ============================================================================
// poly_ntt — Cooley–Tukey 正向 NTT (反循环卷积)
// ============================================================================

Status poly_ntt(const Poly& a, const NttTable& tbl, NttPoly* out) {
    if (!out) {
        return {ErrorCode::InvalidParams, "poly_ntt: out is null"};
    }
    if (a.n() != tbl.n) {
        return {ErrorCode::InvalidParams, "poly_ntt: dimension mismatch"};
    }

    const int n = tbl.n;
    const BarrettConst& bc = tbl.barrett;

    // 复制输入并模约减到 canonical
    std::vector<int64_t> f(n);
    for (int i = 0; i < n; ++i) {
        f[i] = barrett_reduce(a[i], bc);
    }

    // ── 预扭 (pre-twist): f[i] *= ψ^i ──
    // 将反循环卷积 (mod X^n+1) 转换为标准循环卷积 (mod X^n-1)
    int64_t psi = tbl.zetas[1];
    int64_t pow_psi = 1;
    for (int i = 0; i < n; ++i) {
        f[i] = barrett_reduce(f[i] * pow_psi, bc);
        pow_psi = barrett_reduce(pow_psi * psi, bc);
    }

    // ── Bit-reversal permutation (使用预计算表) ──
    for (int i = 0; i < n; ++i) {
        int rev = tbl.bitrev[i];
        if (i < rev) {
            std::swap(f[i], f[rev]);
        }
    }

    // ── Cooley–Tukey NTT (DIT), ω = ψ^2 ──
    for (int len = 2; len <= n; len <<= 1) {
        int half = len / 2;
        for (int start = 0; start < n; start += len) {
            for (int j = 0; j < half; ++j) {
                // zeta = ω^{j * n/len} = ψ^{2*j*n/len}
                int zeta_idx = 2 * j * (n / len);
                int64_t zeta = tbl.zetas[zeta_idx];

                int64_t t = barrett_reduce(zeta * f[start + half + j], bc);
                f[start + half + j] = barrett_reduce(f[start + j] - t, bc);
                f[start + j]        = barrett_reduce(f[start + j] + t, bc);
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

Status poly_invntt(const NttPoly& a, const NttTable& tbl, Poly* out) {
    if (!out) {
        return {ErrorCode::InvalidParams, "poly_invntt: out is null"};
    }
    if (a.n() != tbl.n) {
        return {ErrorCode::InvalidParams, "poly_invntt: dimension mismatch"};
    }

    const int n = tbl.n;
    const BarrettConst& bc = tbl.barrett;

    // 复制输入
    std::vector<int64_t> f(n);
    for (int i = 0; i < n; ++i) {
        f[i] = barrett_reduce(a[i], bc);
    }

    // ── Gentleman–Sande INTT, ω^{-1} = ψ^{-2} ──
    for (int len = n; len >= 2; len >>= 1) {
        int half = len / 2;
        for (int start = 0; start < n; start += len) {
            for (int j = 0; j < half; ++j) {
                int64_t u = f[start + j];
                int64_t v = f[start + half + j];

                f[start + j]        = barrett_reduce(u + v, bc);
                f[start + half + j] = barrett_reduce(u - v, bc);

                // ψ^{-2*j*n/len}
                int zeta_inv_idx = 2 * j * (n / len);
                f[start + half + j] = barrett_reduce(
                    tbl.zetas_inv[zeta_inv_idx] * f[start + half + j], bc);
            }
        }
    }

    // ── Bit-reversal permutation (使用预计算表) ──
    for (int i = 0; i < n; ++i) {
        int rev = tbl.bitrev[i];
        if (i < rev) {
            std::swap(f[i], f[rev]);
        }
    }

    // ── 缩放因子 n^{-1} mod q ──
    const int64_t n_inv = tbl.n_inv;

    // ── 后扭 (post-twist): f[i] *= ψ^{-i} ──
    // ψ^{-1} = zetas[2n-1] = zetas_inv[1]
    int64_t psi_inv = tbl.zetas_inv[1];
    int64_t pow_psi_inv = 1;

    *out = Poly(n);
    for (int i = 0; i < n; ++i) {
        int64_t scaled = barrett_reduce(n_inv * f[i], bc);
        (*out)[i] = barrett_reduce(scaled * pow_psi_inv, bc);
        pow_psi_inv = barrett_reduce(pow_psi_inv * psi_inv, bc);
    }

    return Status::Ok();
}

// ============================================================================
// poly_pointwise_mul_ntt
// ============================================================================

Status poly_pointwise_mul_ntt(const NttPoly& a, const NttPoly& b,
                              const NttTable& tbl, NttPoly* out) {
    if (!out) {
        return {ErrorCode::InvalidParams, "poly_pointwise_mul_ntt: out is null"};
    }
    if (a.n() != b.n() || a.n() != tbl.n) {
        return {ErrorCode::InvalidParams,
                "poly_pointwise_mul_ntt: dimension mismatch"};
    }

    const int n = tbl.n;
    const BarrettConst& bc = tbl.barrett;
    *out = NttPoly(n);
    for (int i = 0; i < n; ++i) {
        (*out)[i] = barrett_reduce(a[i] * b[i], bc);
    }
    return Status::Ok();
}

// ============================================================================
// ntt_roundtrip_test
// ============================================================================

Status ntt_roundtrip_test(const Poly& a, const NttTable& tbl) {
    NttPoly a_ntt(tbl.n);
    auto status = poly_ntt(a, tbl, &a_ntt);
    if (!status.ok()) return status;

    Poly a_recovered(tbl.n);
    status = poly_invntt(a_ntt, tbl, &a_recovered);
    if (!status.ok()) return status;

    // 将原始输入规约到 canonical [0, q) 用于比较
    Poly a_canon(tbl.n);
    int64_t q = static_cast<int64_t>(tbl.q);
    for (int i = 0; i < tbl.n; ++i) {
        int64_t c = a[i] % q;
        if (c < 0) c += q;
        a_canon[i] = c;
    }

    if (!poly_equal_ct(a_canon, a_recovered)) {
        return {ErrorCode::InternalError,
                "ntt_roundtrip_test: INTT(NTT(a)) != a"};
    }
    return Status::Ok();
}

// ============================================================================
// poly_mul_ntt — 完整 NTT 乘法管线
// ============================================================================

Status poly_mul_ntt(const Poly& a, const Poly& b,
                    const NttTable& tbl, Poly* out) {
    if (!out) {
        return {ErrorCode::InvalidParams, "poly_mul_ntt: out is null"};
    }
    if (a.n() != b.n() || a.n() != tbl.n) {
        return {ErrorCode::InvalidParams, "poly_mul_ntt: dimension mismatch"};
    }

    NttPoly ntt_a(tbl.n), ntt_b(tbl.n), ntt_c(tbl.n);

    auto s = poly_ntt(a, tbl, &ntt_a);
    if (!s.ok()) return s;
    s = poly_ntt(b, tbl, &ntt_b);
    if (!s.ok()) return s;
    s = poly_pointwise_mul_ntt(ntt_a, ntt_b, tbl, &ntt_c);
    if (!s.ok()) return s;
    return poly_invntt(ntt_c, tbl, out);
}

} // namespace ibags
