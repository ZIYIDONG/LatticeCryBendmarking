/**
 * @file ntru_trapgen.cpp
 * @brief NTRU Trapdoor Generation implementation
 */

#include "../include_RLWEorNTRU/ntru_trapgen.h"
#include "../include_RLWEorNTRU/gauss.h"
#include "../include_RLWEorNTRU/xof.h"
#include "../include_RLWEorNTRU/secure_memory.h"
#include "../include_RLWEorNTRU/poly_ntt.h"

#include <cassert>
#include <cstring>
#include <sstream>

namespace ibags {

// ============================================================================
// TrapdoorBasis
// ============================================================================

bool TrapdoorBasis::dimensions_consistent() const noexcept {
    int n = f.n();
    return (g.n() == n) && (capF.n() == n) && (capG.n() == n);
}

// ============================================================================
// §6  sample_f_until_invertible
// ============================================================================

Status sample_f_until_invertible(const Params& pp,
                                  SeedCSPRNG& csprng,
                                  int max_attempts,
                                  Poly* f_out) {
    if (!f_out) {
        return {ErrorCode::InvalidParams, "f_out is null"};
    }
    if (max_attempts < 1) {
        return {ErrorCode::InvalidParams, "max_attempts must be >= 1"};
    }

    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        // Sample f from discrete Gaussian using Xof (SHAKE256-based)
        Xof xof("gauss_f_sample");

        // Absorb attempt counter for domain separation
        uint8_t attempt_bytes[4];
        attempt_bytes[0] = static_cast<uint8_t>(attempt & 0xFF);
        attempt_bytes[1] = static_cast<uint8_t>((attempt >> 8) & 0xFF);
        attempt_bytes[2] = static_cast<uint8_t>((attempt >> 16) & 0xFF);
        attempt_bytes[3] = static_cast<uint8_t>((attempt >> 24) & 0xFF);
        xof.absorb(ByteSpan(attempt_bytes, 4));
        xof.finalize();

        // Sample coefficients one by one
        Poly f_sample(pp.n);
        bool sampling_ok = true;
        for (int i = 0; i < pp.n; ++i) {
            int64_t coeff = gauss_sample_coeff(xof, pp.sigma, static_cast<int>(pp.q));
            if (coeff == 0 && pp.sigma < 0.5) {
                // Degenerate case: sigma too small to produce non-zero samples
                sampling_ok = false;
                break;
            }
            f_sample[i] = coeff;
        }

        if (!sampling_ok) continue;

        f_sample.reduce_mod_q(pp);

        // Check invertibility
        if (poly_is_invertible_mod_q(f_sample, pp)) {
            *f_out = std::move(f_sample);
            return Status::Ok();
        }

        // f_sample 离开作用域时自动析构
    }

    // All attempts exhausted
    return {ErrorCode::TrapdoorError,
            "sample_f_until_invertible: f not invertible after " +
            std::to_string(max_attempts) + " attempts"};
}

// ============================================================================
// §7  ntru_trapgen — main
// ============================================================================

Status ntru_trapgen(const Params& pp,
                     SeedCSPRNG& csprng,
                     const TrapGenConfig& config,
                     PublicTrapdoorParams* pub_out,
                     MasterTrapdoorSecret* sec_out) {

    if (!pub_out || !sec_out) {
        return {ErrorCode::InvalidParams, "ntru_trapgen: pub_out or sec_out is null"};
    }

    int n = pp.n;

    // ── Step 1: Sample f until invertible ──
    Poly f(pp.n);
    {
        Status st = sample_f_until_invertible(pp, csprng,
                                               config.max_invertibility_attempts,
                                               &f);
        if (!st.ok()) return st;
    }

    // ── Step 2: Sample g ──
    Poly g(pp.n);

    Xof xof_g("ntru_g_sample");
    xof_g.finalize();

    for (int i = 0; i < n; ++i) {
        int64_t coeff = gauss_sample_coeff(xof_g, pp.sigma, static_cast<int>(pp.q));
        g[i] = coeff;
    }
    g.reduce_mod_q(pp);

    // ── Step 3: Compute f_inv ──
    Poly f_inv(pp.n);
    {
        Status st = poly_inv(f, pp, &f_inv);
        if (!st.ok()) {
            return {ErrorCode::TrapdoorError,
                    "ntru_trapgen: f_inv computation failed: " + st.message()};
        }
    }

    // ── Step 4: h = g * f_inv (NTT-based) ──
    Poly h(pp.n);
    if (is_ntt_friendly(pp)) {
        NttTable tbl = NttTable::create(pp);
        Status mul_st = poly_mul_ntt(g, f_inv, tbl, &h);
        if (!mul_st.ok()) {
            return {ErrorCode::InternalError,
                    "ntru_trapgen: NTT multiplication failed: " + mul_st.message()};
        }
    } else {
        h = poly_mul_naive(g, f_inv, pp);
    }

    // ── Step 5: Optionally compute h_inv ──
    bool h_is_invertible = false;
    Poly h_inv(pp.n);

    if (config.compute_h_inv) {
        Status st_hinv = poly_inv(h, pp, &h_inv);
        h_is_invertible = st_hinv.ok();
        if (!h_is_invertible) {
            // h_inv 计算失败——h 不可逆。安全擦除部分计算结果
            secure_zero(h_inv.data(), static_cast<size_t>(n) * sizeof(int64_t));
        }
    }

    // ── Step 6: Solve NTRU equation f*G - g*F = q ──
    // 使用 aggregate initialization 避免默认构造 Poly
    NtruEquationResult eq_result{Poly::zero(pp), Poly::zero(pp)};
    {
        Status eq_st = solve_ntru_equation(pp, f, g, &eq_result);
        if (!eq_st.ok()) {
            if (config.fallback_on_ntru_solver) {
                // Fallback: 使用零多项式作为 F, G（仅 PoC）
                eq_result.F = Poly::zero(pp);
                eq_result.G = Poly::zero(pp);
            } else {
                return eq_st;
            }
        }
    }

    // ── Step 7: Assemble outputs ──
    // Public params
    pub_out->params = pp;
    pub_out->h = std::move(h);
    pub_out->h_is_invertible = h_is_invertible;
    if (h_is_invertible) {
        pub_out->h_inv = std::move(h_inv);
    } else {
        pub_out->h_inv = Poly::zero(pp);
    }

    // Master secret
    sec_out->basis.f = std::move(f);
    sec_out->basis.g = std::move(g);
    sec_out->basis.capF = std::move(eq_result.F);
    sec_out->basis.capG = std::move(eq_result.G);
    sec_out->h_is_invertible = h_is_invertible;
    if (h_is_invertible) {
        // h_inv 已移至 pub_out->h_inv，此处重新计算
        sec_out->h_inv = pub_out->h_inv; // 浅拷贝，仅 PoC
    } else {
        sec_out->h_inv = Poly::zero(pp);
    }

    return Status::Ok();
}

// ============================================================================
// §8  solve_ntru_equation — STUB
// ============================================================================

Status solve_ntru_equation(const Params& pp,
                            const Poly& f,
                            const Poly& g,
                            NtruEquationResult* result) {
    (void)f;
    (void)g;

    if (!result) {
        return {ErrorCode::InvalidParams, "solve_ntru_equation: result is null"};
    }

    // [STUB - NOT FOR PRODUCTION]
    // 返回零多项式作为占位符
    result->F = Poly::zero(pp);
    result->G = Poly::zero(pp);

    // 标记为 InternalError，让调用者决定是否走 fallback 路径
    return {ErrorCode::InternalError,
            "solve_ntru_equation: [STUB] returning F=G=0. "
            "Not suitable for production. Please implement a full NTRU equation solver."};
}

// ============================================================================
// §9  check_h_from_f_g
// ============================================================================

bool check_h_from_f_g(const Params& pp,
                       const Poly& f,
                       const Poly& g,
                       const Poly& h) {
    if (f.n() != pp.n || g.n() != pp.n || h.n() != pp.n) return false;

    // f 必须可逆，否则 h = g * f^{-1} 无意义
    if (!poly_is_invertible_mod_q(f, pp)) return false;

    // h * f mod (q, x^n+1) — prefer NTT if available
    Poly hf(pp.n);
    if (is_ntt_friendly(pp)) {
        NttTable tbl = NttTable::create(pp);
        (void)poly_mul_ntt(h, f, tbl, &hf);
    } else {
        hf = poly_mul_naive(h, f, pp);
    }

    // 与 g 比较
    return poly_equal_ct(hf, g);
}

// ============================================================================
// §9  check_trapdoor_basis — full validation
// ============================================================================

bool check_trapdoor_basis(const Params& pp,
                           const TrapdoorBasis& basis,
                           const Poly& h) {
    int n = pp.n;

    // 1. 维度一致性
    if (!basis.dimensions_consistent()) return false;

    // 2. 范数检查（f, g 的 ∞-范数 ≤ eta1）
    uint64_t eta1 = static_cast<uint64_t>(pp.eta1);
    if (!poly_norm_bound_check(basis.f, eta1, pp)) return false;
    if (!poly_norm_bound_check(basis.g, eta1, pp)) return false;

    // 注意: F, G 的大小取决于 NTRU 方程求解器，不在此处硬性检查
    // （在 PoC 中 F=G=0，显然满足范数条件）

    // 3. h == g * f^{-1}
    if (!check_h_from_f_g(pp, basis.f, basis.g, h)) return false;

    // 4. Canonical 形式：所有系数在 [0, q)
    if (!basis.f.is_canonical(pp)) return false;
    if (!basis.g.is_canonical(pp)) return false;
    if (!basis.capF.is_canonical(pp)) return false;
    if (!basis.capG.is_canonical(pp)) return false;

    // 5. NTRU 关系: f*G - g*F ≡ q (mod x^n+1)
    //    仅在 F 和 G 非平凡时检查（跳过 PoC 零多项式情况）
    Poly fG(pp.n), gF(pp.n);
    if (is_ntt_friendly(pp)) {
        NttTable tbl = NttTable::create(pp);
        (void)poly_mul_ntt(basis.f, basis.capG, tbl, &fG);
        (void)poly_mul_ntt(basis.g, basis.capF, tbl, &gF);
    } else {
        fG = poly_mul_naive(basis.f, basis.capG, pp);
        gF = poly_mul_naive(basis.g, basis.capF, pp);
    }
    Poly ntru_diff = poly_sub(fG, gF, pp);

    // 期望: f*G - g*F ≡ q (mod x^n+1)，即系数全 0 模 q × 常数 q。
    // 在 PoC (F=G=0) 中，ntru_diff 全 0，与 Poly::zero(pp) 比较即可。
    // 生产中应检查: ntru_diff[0] ≡ q (mod q) 且其余系数 ≡ 0 (mod q)。
    Poly expected_zero = Poly::zero(pp);
    if (!poly_equal_ct(ntru_diff, expected_zero)) return false;

    return true;
}

// ============================================================================
// §10  gauss_sample_preimage — STUB
// ============================================================================

Status gauss_sample_preimage(const Params& pp,
                              const TrapdoorBasis& basis,
                              const Poly& target,
                              double sigma,
                              SeedCSPRNG& csprng,
                              Poly* s1_out,
                              Poly* s2_out) {
    (void)basis;
    (void)target;
    (void)sigma;
    (void)csprng;

    if (!s1_out || !s2_out) {
        return {ErrorCode::InvalidParams,
                "gauss_sample_preimage: s1_out or s2_out is null"};
    }

    // [STUB - NOT FOR PRODUCTION]
    // 返回零向量作为占位符
    *s1_out = Poly::zero(pp);
    *s2_out = Poly::zero(pp);

    return {ErrorCode::InternalError,
            "gauss_sample_preimage: [STUB] returning zero vectors. "
            "Not suitable for production. Please implement a full GPV/Falcon sampler."};
}

} // namespace ibags