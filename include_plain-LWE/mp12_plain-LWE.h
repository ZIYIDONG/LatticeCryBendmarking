#pragma once
/**
 * MP12 Trapdoor Generation (Micciancio-Peikert 2012)
 * "Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller"
 * EUROCRYPT 2012 / IACR ePrint 2011/501
 *
 * Core construction:
 *   Public matrix A = [Ā | G - Ā·R]  ∈ Z_q^{n×m}
 *   Trapdoor        R                  ∈ Z^{nk×nk}  (short random)
 *   Gadget matrix   G = I_n ⊗ g^T,    g = (1,b,b²,...,b^{k-1})
 *
 * This file provides:
 *   - Matrix/vector primitives (mod-q arithmetic)
 *   - Gadget G and its basis S_g
 *   - Discrete Gaussian sampler over Z (Knuth-Yao / rejection)
 *   - GenTrap  : generate (A, trapdoor T)
 *   - SamplePre: given T, sample short x s.t. A·x = u (mod q)
 */

#include <cstdint>
#include <cstddef>
#include <vector>
#include <random>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <functional>



namespace mp12 {

/* ─────────────────────────── parameters ─────────────────────────── */
struct Params {
    int n;        // lattice dimension
    long q;       // modulus
    int b;        // gadget base  (typically 2)
    int k;        // k = ceil(log_b q)
    int m_bar;    // columns of Ā  (≥ n)
    int m;        // total columns of A = m_bar + n*k
    double sigma; // Gaussian width for trapdoor entries (≈ 1)
    double s;     // Gaussian width for preimage sampling

    static Params make(int n, long q, int b = 2) {
        Params p;
        p.n     = n;
        p.q     = q;
        p.b     = b;
        p.k     = (int)std::ceil(std::log((double)q) / std::log((double)b));
        p.m_bar = n;                    // Ā is n×n (square, random)
        p.m     = p.m_bar + n * p.k;   // full width of A
#if defined(MP12_SIGMA)
        p.sigma = MP12_SIGMA;
#else
        p.sigma = 1.0;                  // trapdoor entry distribution width
#endif
        // Preimage Gaussian width: s > σ_1(T) · √(n·k + n) · ω(√log n)
        // A conservative choice: s = b·k·n  (provably secure for moderate n)
        p.s     = (double)b * p.k * std::sqrt((double)n) * 5.0;
        return p;
    }
};

/* ─────────────────────────── matrix type ────────────────────────── */
using Vec = std::vector<long>;
using Mat = std::vector<Vec>;   // Mat[row][col]

inline Mat make_mat(int rows, int cols, long fill = 0) {
    return Mat(rows, Vec(cols, fill));
}

inline long mod(long x, long q) {
    return ((x % q) + q) % q;
}

// C = A*B mod q
inline Mat mat_mul_mod(const Mat& A, const Mat& B, long q) {
    int r = A.size(), n = B.size(), c = B[0].size();
    Mat C = make_mat(r, c);
    for (int i = 0; i < r; i++)
        for (int k = 0; k < n; k++) {
            if (A[i][k] == 0) continue;
            for (int j = 0; j < c; j++)
                C[i][j] = mod(C[i][j] + A[i][k] * B[k][j], q);
        }
    return C;
}

// b = A*x mod q
inline Vec mat_vec_mod(const Mat& A, const Vec& x, long q) {
    int r = A.size(), c = A[0].size();
    Vec b(r, 0);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            b[i] = mod(b[i] + A[i][j] * x[j], q);
    return b;
}

// A - B mod q
inline Mat mat_sub_mod(const Mat& A, const Mat& B, long q) {
    int r = A.size(), c = A[0].size();
    Mat C = make_mat(r, c);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            C[i][j] = mod(A[i][j] - B[i][j], q);
    return C;
}

// horizontal concat [A | B]
inline Mat mat_hcat(const Mat& A, const Mat& B) {
    int r = A.size(), ca = A[0].size(), cb = B[0].size();
    Mat C = make_mat(r, ca + cb);
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < ca; j++) C[i][j]      = A[i][j];
        for (int j = 0; j < cb; j++) C[i][ca + j] = B[i][j];
    }
    return C;
}

// vector norm squared (over integers)
inline double vec_norm2(const Vec& v) {
    double s = 0;
    for (auto x : v) s += (double)x * x;
    return s;
}

/* ─────────────────────── Gadget matrix G ────────────────────────── */
/**
 * G = I_n ⊗ g^T  where g = (1, b, b², ..., b^{k-1})
 * Size: n × (n·k)
 * G[i][i*k + j] = b^j
 */
inline Mat gadget_matrix(const Params& p) {
    Mat G = make_mat(p.n, p.n * p.k);
    long bpow = 1;
    for (int j = 0; j < p.k; j++) {
        for (int i = 0; i < p.n; i++)
            G[i][i * p.k + j] = bpow % p.q;
        bpow = (bpow * p.b) % p.q;
    }
    return G;
}

/**
 * Basis S_g for Λ^⊥(g^T) in Z^k
 * This is the "gadget basis" – a (k×k) integer matrix S_g s.t. g^T · S_g = 0 mod q
 * For binary gadget (b=2):
 *   S_g = [ -2  0  0 ...  0 ]
 *         [  1 -2  0 ...  0 ]
 *         [  0  1 -2 ...  0 ]
 *         [  ...             ]
 *         [  0  0  0 ... -2 ]
 *   plus last column compensating for q = 2^k (or nearest).
 * For general base b, S_g has the same structure with (-b) on diagonal.
 *
 * We store the block-diagonal matrix S (n·k × n·k).
 */
inline Mat gadget_basis(const Params& p) {
    /**
     * Correct construction for ARBITRARY q (MP12, Lemma 5).
     *
     * Scalar gadget: g = (1, b, b², ..., b^{k-1})  ∈ Z^k
     * We need k×k integer matrix S_k s.t.  g^T · S_k ≡ 0  (mod q).
     *
     * Columns 0 … k-2  (shift columns):
     *   S_k[j][j]   = -b      ← diagonal
     *   S_k[j+1][j] =  1      ← sub-diagonal
     *   Verification: g^T · col_j = b^j·(-b) + b^{j+1}·1 = 0  ✓
     *
     * Column k-1  (wrap column = base-b digits of q):
     *   S_k[i][k-1] = q_i,   where q = Σ q_i · b^i  is the base-b expansion
     *   Verification: g^T · col_{k-1} = Σ b^i · q_i = q ≡ 0 (mod q)  ✓
     *
     * This works for ANY q, including non-powers-of-b (e.g. prime q).
     * When q is an exact power of b all digits are 0, giving the usual 0 last col.
     */

    // Step 1: base-b digits of q
    std::vector<long> q_digits(p.k, 0);
    long tmp = p.q;
    for (int i = 0; i < p.k; i++) {
        q_digits[i] = tmp % p.b;
        tmp /= p.b;
    }
    // tmp == 0 here because k = ceil(log_b q)

    // Step 2: build scalar k×k block
    Mat Sk = make_mat(p.k, p.k);
    for (int j = 0; j < p.k - 1; j++) {
        Sk[j][j]   = -p.b;   // shift: -b on diagonal
        Sk[j+1][j] =  1;     // shift:  1 on sub-diagonal
    }
    for (int i = 0; i < p.k; i++)
        Sk[i][p.k - 1] = q_digits[i];  // wrap: base-b digits of q

    // Step 3: block-diagonal S = diag(Sk, …, Sk) of size (nk × nk)
    int nk = p.n * p.k;
    Mat S = make_mat(nk, nk);
    for (int blk = 0; blk < p.n; blk++)
        for (int i = 0; i < p.k; i++)
            for (int j = 0; j < p.k; j++)
                S[blk * p.k + i][blk * p.k + j] = Sk[i][j];
    return S;
}

/* ──────────────────── Discrete Gaussian sampler ─────────────────── */
/**
 * Sample x ~ D_{Z, sigma, center} using rejection sampling.
 * Samples from a discrete Gaussian with standard deviation sigma centered at c.
 */
class DGSampler {
public:
    explicit DGSampler(double sigma, double center = 0.0, uint64_t seed = 0)
        : sigma_(sigma), center_(center) {
        if (seed == 0) {
            std::random_device rd;
            rng_.seed(rd());
        } else {
            rng_.seed(seed);
        }
        tail_bound_ = (int)std::ceil(sigma * 12.0); // 12σ tail cut
    }

    long sample() {
        // Rejection sampling: propose from Z in [c-tail, c+tail],
        // accept with probability proportional to exp(-x²/(2σ²))
        std::uniform_int_distribution<long> uniform(
            (long)std::floor(center_) - tail_bound_,
            (long)std::ceil(center_)  + tail_bound_);
        std::uniform_real_distribution<double> unif01(0.0, 1.0);

        while (true) {
            long x = uniform(rng_);
            double dist = (double)x - center_;
            double prob = std::exp(-dist * dist / (2.0 * sigma_ * sigma_));
            if (unif01(rng_) < prob)
                return x;
        }
    }

    Vec sample_vec(int n) {
        Vec v(n);
        for (int i = 0; i < n; i++) v[i] = sample();
        return v;
    }

    Mat sample_mat(int rows, int cols) {
        Mat M = make_mat(rows, cols);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                M[i][j] = sample();
        return M;
    }

private:
    double sigma_, center_;
    int tail_bound_;
    std::mt19937_64 rng_;
};

/* ───────────────────── Uniform mod-q sampler ────────────────────── */
class UniformSampler {
public:
    explicit UniformSampler(long q, uint64_t seed = 0) : q_(q) {
        if (seed == 0) {
            std::random_device rd;
            rng_.seed(rd());
        } else {
            rng_.seed(seed);
        }
    }
    long sample() {
        std::uniform_int_distribution<long> d(0, q_ - 1);
        return d(rng_);
    }
    Mat sample_mat(int rows, int cols) {
        Mat M = make_mat(rows, cols);
        std::uniform_int_distribution<long> d(0, q_ - 1);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                M[i][j] = d(rng_);
        return M;
    }
private:
    long q_;
    std::mt19937_64 rng_;
};

/* ────────────────────────── Trapdoor type ───────────────────────── */
struct Trapdoor {
    Mat R;    // n×(nk) short integer matrix (the actual trapdoor)
    Mat A;    // n×m public matrix
};

/* ─────────────────── GenTrap: Algorithm 1 of MP12 ──────────────── */
/**
 * GenTrap(1^n, 1^m, q):
 *   1. Sample Ā ← U(Z_q^{n × m_bar})
 *   2. Sample R ← D_{Z,σ}^{m_bar × nk}  (short random matrix)
 *   3. G = gadget matrix (n × nk)
 *   4. A_R = Ā·R mod q
 *   5. A = [Ā | G - A_R]   ∈ Z_q^{n × m}
 *
 * Trapdoor: R, because A · [R; I_{nk}] = Ā·R + G - Ā·R = G (mod q)
 * So T = [R; I] is a basis for Λ^⊥(A) relative to G.
 */
inline Trapdoor gen_trap(const Params& p, uint64_t seed = 0) {
    UniformSampler usampler(p.q, seed);
    DGSampler      dsampler(p.sigma, 0.0, seed ? seed + 1 : 0);

    // Step 1: random Ā ∈ Z_q^{n × m_bar}
    Mat A_bar = usampler.sample_mat(p.n, p.m_bar);

    // Step 2: short R ∈ Z^{m_bar × nk}
    int nk = p.n * p.k;
    Mat R = dsampler.sample_mat(p.m_bar, nk);

    // Step 3: G ∈ Z_q^{n × nk}
    Mat G = gadget_matrix(p);

    // Step 4: Ā·R mod q
    Mat AR = mat_mul_mod(A_bar, R, p.q);

    // Step 5: G - Ā·R mod q
    Mat G_minus_AR = mat_sub_mod(G, AR, p.q);

    // A = [Ā | G - Ā·R]
    Mat A = mat_hcat(A_bar, G_minus_AR);

    return Trapdoor{std::move(R), std::move(A)};
}

/* ──────────── G-lattice preimage sampling (SampleG) ────────────── */
/**
 * Given target vector u ∈ Z_q^n, sample z ∈ Z^{nk} s.t. G·z = u (mod q).
 *
 * Uses the gadget basis S_g: solve each scalar equation
 *   g^T · z_i = u_i (mod q)
 * via "balanced" base-b decomposition.
 *
 * Returns z ∈ Z^{nk} (exact preimage, NOT Gaussian; for a Gaussian version
 * one would add a perturbation from the coset; here we give the deterministic
 * short preimage used in the proof-of-concept setting).
 */
inline Vec sample_g(const Params& p, const Vec& u) {
    // For each block i, we want g^T · z_i = u_i (mod q)
    // where g = (1,b,b²,...,b^{k-1}) and z_i ∈ Z^k.
    // Short solution: base-b representation of u_i
    //   z_i[0] = u_i mod b (possibly negative centered)
    //   carry   = (u_i - z_i[0]) / b
    //   z_i[j] = carry mod b, carry >>= b, etc.
    Vec z(p.n * p.k, 0);
    for (int i = 0; i < p.n; i++) {
        long val = u[i];
        for (int j = 0; j < p.k; j++) {
            long digit = val % p.b;
            // Balanced representation: center in [-b/2, b/2)
            if (digit > p.b / 2) digit -= p.b;
            if (digit < -(p.b / 2)) digit += p.b;
            z[i * p.k + j] = digit;
            val = (val - digit) / p.b;
        }
        // val should be 0 now (mod q/b^k ≈ 1); if not, absorb into last digit
        if (val != 0)
            z[i * p.k + p.k - 1] += val;
    }
    return z;
}

/* ──────────────── SamplePre: Algorithm 2 of MP12 ───────────────── */
/**
 * SamplePre(A, T, u, s):
 *   Sample x ← D_{Λ_u^⊥(A), s}  s.t. A·x = u (mod q)
 *
 * Strategy (simplified offline-online decomposition):
 *   1. The trapdoor T = (R, I_{nk}) satisfies A·T = G (mod q).
 *   2. To find x s.t. A·x = u:
 *      a. Find z s.t. G·z = u (mod q)  via SampleG
 *      b. We need y s.t. A·y = 0 and output x = T·z + y  -- but
 *         instead we use the direct Gaussian perturbation:
 *
 *   Full perturbation sampling (MP12 §4):
 *      p  ← SamplePerturb(T, s)   (Gaussian perturb. hiding T)
 *      v  = A·p mod q
 *      z  ← SampleG(G, u - v)
 *      x  = p + T·z
 *
 *   Here we implement the simplified version without full perturbation
 *   (proof-of-concept; add proper perturbation for secure use).
 *
 * Output: x ∈ Z^m s.t. A·x = u (mod q) and ||x|| is small.
 */
inline Vec sample_pre(const Params& p, const Trapdoor& td,
                      const Vec& u, uint64_t seed = 0)
{
    DGSampler dsampler(p.s, 0.0, seed);

    // --- Perturbation step ---
    // p_vec ∈ Z^m, sampled from discrete Gaussian
    // (In the full version this would use the trapdoor Gram-Schmidt;
    //  here we use a simple spherical Gaussian as a placeholder.)
    Vec p_vec = dsampler.sample_vec(p.m);

    // v = A · p_vec mod q
    Vec v = mat_vec_mod(td.A, p_vec, p.q);

    // --- G-lattice sampling ---
    // Compute u - v mod q
    Vec u_minus_v(p.n);
    for (int i = 0; i < p.n; i++)
        u_minus_v[i] = mod(u[i] - v[i], p.q);

    // z ∈ Z^{nk} s.t. G·z = u - v (mod q)
    Vec z = sample_g(p, u_minus_v);

    // --- Combine ---
    // x = p_vec + [R; I_{nk}]·z
    // T = [R; I_{nk}] is (m_bar + nk) × nk, so T·z ∈ Z^m
    Vec Tz(p.m, 0);
    int nk = p.n * p.k;
    // Upper part: R·z  (m_bar rows)
    for (int i = 0; i < p.m_bar; i++)
        for (int j = 0; j < nk; j++)
            Tz[i] += td.R[i][j] * z[j];
    // Lower part: I·z  (nk rows)
    for (int j = 0; j < nk; j++)
        Tz[p.m_bar + j] = z[j];

    Vec x(p.m);
    for (int i = 0; i < p.m; i++)
        x[i] = p_vec[i] + Tz[i];

    return x;
}

/* ──────────────────────── Verification ─────────────────────────── */
inline bool verify(const Params& p, const Mat& A, const Vec& x, const Vec& u) {
    Vec Ax = mat_vec_mod(A, x, p.q);
    for (int i = 0; i < p.n; i++)
        if (Ax[i] != mod(u[i], p.q)) return false;
    return true;
}

} // namespace mp12