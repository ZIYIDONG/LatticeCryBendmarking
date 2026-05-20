/**
 * @file ntru_secret.cpp
 * @brief NTRU trapdoor secret polynomial sampling (f, g)
 *
 * Implements coefficient-wise discrete Gaussian sampling of short
 * polynomials f, g ∈ R_p for NTRU TrapGen (Falcon-style key generation).
 *
 * This is ONLY the sampling step. Full TrapGen additionally:
 *   - checks invertibility of f (and optionally g)
 *   - computes h = g · f⁻¹
 *   - solves the NTRU equation for F, G
 *   - assembles and validates the TrapdoorBasis B
 *
 * Reference:
 *   - Falcon keygen — f, g sampling (Falcon spec §3.1)
 *   - NTRU key generation — f, g distribution
 *
 * Production notes:
 *   - PoC: Uses coefficient-wise discrete Gaussian from ibags::gauss module.
 *     The resulting f, g may NOT be invertible with the same probability
 *     as Falcon's targeted sampling (Falcon tailors bit-width of f, g).
 *   - Production: Replace with Falcon-derived sampler that ensures
 *     invertibility and matches exact NTRU distribution (binary/ternary
 *     near center).
 */

#include "../include_RLWEorNTRU/ntru_secret.h"
#include <cassert>
#include <vector>
#include <array>

namespace ibags {

// ============================================================================
// Internal: derive a 32-byte seed from SeedCSPRNG
// ============================================================================

namespace {

/**
 * @brief Extract 32 deterministic bytes from a SeedCSPRNG.
 *
 * SeedCSPRNG has no seed_bytes() accessor, so we request 32 pseudorandom
 * bytes which are deterministic given the same construction parameters.
 * These bytes are then fed to build_gauss_xof to domain-separate f vs g.
 */
std::array<uint8_t, 32> extract_seed(SeedCSPRNG& csprng) {
    std::array<uint8_t, 32> buf{};
    csprng.randombytes(buf.data(), buf.size());
    return buf;
}

} // anonymous namespace

// ============================================================================
// Public API
// ============================================================================

SecretPoly sample_ntru_secret_f(const Params& pp, SeedCSPRNG& csprng) {
    assert(validate_params(pp).ok() && "Params must be valid before sampling");

    // Domain-separated XOF: "IBAGS-v1/NTRU-F" ensures f and g samples are
    // independent even when the CSPRNG state is shared.
    auto seed = extract_seed(csprng);
    ByteSpan seed_span(seed.data(), seed.size());
    Xof xof = build_gauss_xof("IBAGS-v1/NTRU-F", seed_span);
    SmallPoly small = gauss_sample_poly(xof, pp);

    // Convert SmallPoly (int16) → SecretPoly (int64, auto-zero-on-destroy)
    SecretPoly secret(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        secret.set_coeff(i, static_cast<long>(small[i]));
    }

    // Secure-wipe the intermediate SmallPoly (defense-in-depth)
    // safe: SmallPoly owns the buffer; we cast away const because the
    // underlying bytes are writable.
    secure_zero(const_cast<void*>(
                    static_cast<const void*>(small.coeffs().data())),
                static_cast<size_t>(pp.n) * sizeof(int16_t));

    return secret;
}

SecretPoly sample_ntru_secret_g(const Params& pp, SeedCSPRNG& csprng) {
    assert(validate_params(pp).ok() && "Params must be valid before sampling");

    // Independent domain tag: "IBAGS-v1/NTRU-G"
    auto seed = extract_seed(csprng);
    ByteSpan seed_span(seed.data(), seed.size());
    Xof xof = build_gauss_xof("IBAGS-v1/NTRU-G", seed_span);
    SmallPoly small = gauss_sample_poly(xof, pp);

    SecretPoly secret(pp.n);
    for (int i = 0; i < pp.n; ++i) {
        secret.set_coeff(i, static_cast<long>(small[i]));
    }

    secure_zero(const_cast<void*>(
                    static_cast<const void*>(small.coeffs().data())),
                static_cast<size_t>(pp.n) * sizeof(int16_t));

    return secret;
}

void sample_ntru_secret_pair(const Params& pp,
                              SeedCSPRNG& csprng,
                              SecretPoly* f_out,
                              SecretPoly* g_out) {
    assert(f_out && g_out && "Output pointers must be non-null");
    assert(validate_params(pp).ok() && "Params must be valid before sampling");

    *f_out = sample_ntru_secret_f(pp, csprng);
    *g_out = sample_ntru_secret_g(pp, csprng);
}

} // namespace ibags