#pragma once
/**
 * @file ntru_secret.h
 * @brief NTRU trapdoor secret polynomial sampling (f, g)
 *
 * Samples short polynomials f, g ∈ R_p for NTRU TrapGen.
 * These are the master secret component of the lattice trapdoor.
 *
 * Algorithm:
 *   f, g ← D_{Z^n, σ} (coefficient-wise discrete Gaussian)
 *   Coefficients stored in centered representation ∈ (−q/2, q/2]
 *
 * Design:
 *  - Parameterized distribution (sigma from Params).
 *  - Deterministic from SeedCSPRNG for reproducibility.
 *  - Returns SecretPoly (auto-zero on destruction).
 *  - This is only the sampling step of TrapGen; full TrapGen additionally
 *    checks invertibility, computes h = g·f⁻¹, and solves the NTRU equation.
 *
 * Reference:
 *  - Falcon keygen — f, g sampling (Falcon spec §3.1)
 *  - NTRU key generation — f, g distribution
 *
 * Production notes:
 *  - PoC: Uses coefficient-wise discrete Gaussian from ibags::gauss module.
 *    The resulting f, g may not be invertible with the same probability
 *    as Falcon's targeted sampling (Falcon uses tailored bit-width f, g).
 *  - Production: Replace with Falcon-derived sampler that ensures invertibility
 *    and matches the exact NTRU distribution (binary/ternary near center).
 */

#include "params.h"
#include "secure_memory.h"
#include "csprng.h"
#include "gauss.h"
#include "poly.h"
#include "errors.h"

namespace ibags {

/**
 * @brief Sample a short polynomial f for NTRU trapdoor.
 *
 * Coefficients are drawn from a discrete Gaussian with width sigma
 * (from Params), stored in SecretPoly.
 *
 * @param pp      Parameter set (sigma, n).
 * @param csprng  Seeded CSPRNG (deterministic if same seed).
 * @return        SecretPoly containing f coefficients.
 */
[[nodiscard]] SecretPoly sample_ntru_secret_f(const Params& pp,
                                               SeedCSPRNG& csprng);

/**
 * @brief Sample a short polynomial g for NTRU trapdoor.
 *
 * Same distribution as f; separate function for clarity.
 *
 * @param pp      Parameter set.
 * @param csprng  Seeded CSPRNG.
 * @return        SecretPoly containing g coefficients.
 */
[[nodiscard]] SecretPoly sample_ntru_secret_g(const Params& pp,
                                               SeedCSPRNG& csprng);

/**
 * @brief Sample both f and g.
 *
 * Convenience for TrapGen which needs both.
 *
 * @param pp       Parameter set.
 * @param csprng   Seeded CSPRNG.
 * @param f_out    Output: secret polynomial f.
 * @param g_out    Output: secret polynomial g.
 */
void sample_ntru_secret_pair(const Params& pp,
                              SeedCSPRNG& csprng,
                              SecretPoly* f_out,
                              SecretPoly* g_out);

} // namespace ibags
