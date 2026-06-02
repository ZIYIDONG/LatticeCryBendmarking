#pragma once
/**
 * unienc.h — UniEnc(pp, id, μ_D) 掩码方案
 *
 * 出自格基 ABE/HIBE/同态加密文献,典型形式:
 *
 *   UniEnc(pp, id, μ_D) → (𝒰, 𝒞)
 *
 *   1) R ← {0,1}^{m×m}
 *      A = pp ∈ Z_q^{N×m},  N = (d+1)·n + 1
 *      𝒞 = A·R + μ_D · G    ∈ Z_q^{N×m}
 *
 *   2) 对 R 的每个比特独立加密:
 *      V^(i,j) = Enc(A, R[i][j])  ∈ Z_q^N
 *      𝒰 = (V^(1,1), …, V^(m,m))  ∈ (Z_q^N)^{m²}
 *
 * 基本操作清单(本文件实现):
 *   ① sample_binary_mat   — R ← {0,1}^{m×m}
 *   ② mat_mul_mod         — A·R
 *   ③ scalar_mat_mul_mod  — μ_D · G
 *   ④ mat_add_mod         — AR + μG
 *   ⑤ lwe_encrypt_bit     — 单比特 LWE 加密
 *   ⑥ unienc              — 把 ①–⑤ 串起来
 */

#include "matops_plain-LWE.h"
#include "mp12_plain-LWE.h"
#include <random>
#include <chrono>
#include <vector>
#include <cmath>

namespace unienc {

using matops::Vec;
using matops::Mat;
using matops::mod_pos;
using matops::make_mat;
using matops::mat_mul;
using matops::mat_add;

/* ══════════════════════════════════════════════════
   §1  参数集
   ══════════════════════════════════════════════════ */
struct Params {
    int  n;        // 格维度
    int  d;        // 层级深度(决定 N)
    int  m;        // R 的边长 (= 公钥列数)
    long q;        // 模数
    int  N;        // 密文行数 = (d+1)·n + 1
    double sigma;  // LWE 噪声标准差

    static Params make(int n, int d, int m, long q, double sigma = 3.2) {
        Params p;
        p.n = n; p.d = d; p.m = m; p.q = q;
        p.N = (d + 1) * n + 1;
        if (sigma > 0.0) {
            p.sigma = sigma;
        } else {
#if defined(UNIENC_SIGMA)
            p.sigma = UNIENC_SIGMA;
#else
            p.sigma = 3.2;
#endif
        }
        return p;
    }
};

/* ══════════════════════════════════════════════════
   §2  操作 ① 二进制矩阵采样
   ══════════════════════════════════════════════════ */
/**
 * R ← {0,1}^{rows × cols},每个元素独立均匀
 * 复杂度: O(rows · cols)
 */
inline Mat sample_binary_mat(int rows, int cols, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> bit(0, 1);
    Mat R = make_mat(rows, cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            R[i][j] = bit(rng);
    return R;
}

/* ══════════════════════════════════════════════════
   §3  操作 ③ 标量乘矩阵
   ══════════════════════════════════════════════════ */
/**
 * C = (μ · M) mod q
 * 复杂度: O(rows · cols)
 *
 * 优化: μ == 0 时直接返回零矩阵,μ == 1 时直接返回 M 的副本
 */
inline Mat scalar_mat_mul(long mu, const Mat& M, long q) {
    int r = (int)M.size();
    int c = r ? (int)M[0].size() : 0;
    Mat C = make_mat(r, c);
    long mu_norm = mod_pos(mu, q);
    if (mu_norm == 0) return C;
    if (mu_norm == 1) {
        for (int i = 0; i < r; i++) C[i] = M[i];
        return C;
    }
    unsigned long long mu_barrett = matops::barrett_mu(q);
    for (int i = 0; i < r; i++) {
        const Vec& mi = M[i];
        Vec& ci = C[i];
        for (int j = 0; j < c; j++)
            ci[j] = matops::barrett_reduce_lwe(mu_norm * mi[j], q, mu_barrett);
    }
    return C;
}

/* ══════════════════════════════════════════════════
   §4  Gadget 矩阵 G  (N × m 形状)
   ══════════════════════════════════════════════════ */
/**
 * 这里 G 的形状必须和 𝒞 一样: N × m
 *
 * 标准的 MP12 gadget 是 n × nk,但 UniEnc 用的是稍微不同的形状。
 * 通常做法是: G = I_N ⊗ g^T 的截断版本,
 * 即 G[i][j] = b^j 当 j < k 且 i 对应正确块时
 *
 * 这里我们提供一个简单的"宽 gadget":每行用 g 向量循环填充,
 * 保证有 gadget 性质且形状是 N × m
 */
inline Mat make_gadget(int N, int m, int b, long q) {
    Mat G = make_mat(N, m);
    long bp = 1;
    for (int j = 0; j < m; j++) {
        int row = j % N;
        G[row][j] = bp % q;
        if ((j + 1) % N == 0) bp = (bp * b) % q;
    }
    return G;
}

/* ══════════════════════════════════════════════════
   §5  操作 ⑤ LWE 单比特加密
   ══════════════════════════════════════════════════ */
/**
 * Enc(pp = A, μ ∈ {0,1}) → V ∈ Z_q^N
 *
 * 标准 LWE 单比特加密:
 *   选 s ← {0,1}^m  (短向量)
 *   选 e ← D_{Z, σ}^N  (高斯噪声)
 *   V = A·s + e + μ·⌊q/2⌋·e_1
 *
 * 这里 e_1 是第一个标准基向量,μ 编码到第一个分量
 *
 * 在 UniEnc 上下文中,V^(i,j) 的具体形式可能因方案而异;
 * 这是最常见的"基础版本"
 */
inline Vec lwe_encrypt_bit(const Mat& A, long mu, long q,
                            double sigma, std::mt19937_64& rng) {
    int N = (int)A.size();
    int m = N ? (int)A[0].size() : 0;

    // ① 短随机 s ∈ {0,1}^m
    std::uniform_int_distribution<int> bit(0, 1);
    Vec s(m);
    for (int j = 0; j < m; j++) s[j] = bit(rng);

    // ② V = A·s  (Barrett 约减)
    Vec V(N, 0);
    unsigned long long mu_bar = matops::barrett_mu(q);
    for (int i = 0; i < N; i++) {
        long acc = 0;
        for (int j = 0; j < m; j++) acc += A[i][j] * s[j];
        V[i] = matops::barrett_reduce_lwe(acc, q, mu_bar);
    }

    // ③ 加高斯噪声 e ∈ Z^N  (条件加减)
    std::normal_distribution<double> gauss(0.0, sigma);
    for (int i = 0; i < N; i++) {
        long ei = (long)std::llround(gauss(rng));
        long sum = V[i] + ei;
        V[i] = matops::mod_add(matops::mod_sub(sum, q), q);
    }

    // ④ 把 μ 编码到第一个分量: 加 μ·⌊q/2⌋
    if ((mu & 1) == 1) {
        long sum = V[0] + q / 2;
        V[0] = (sum >= q) ? sum - q : sum;
    }
    return V;
}

/* ══════════════════════════════════════════════════
   §6  UniEnc 主体  (操作 ⑥ 串起 ①–⑤)
   ══════════════════════════════════════════════════ */
/**
 * UniEnc 的输出结构
 *
 * 𝒞: N × m 矩阵           — 主密文
 * 𝒰: 长度 m² 的列向量数组   — 每个 V^(i,j) 是 N × 1
 *     按行主序: U[i*m + j] = V^(i,j)
 */
struct UniEncOutput {
    Mat C;                  // 𝒞 ∈ Z_q^{N×m}
    std::vector<Vec> U;     // 𝒰: m² 个长度为 N 的列向量
    Mat R_internal;         // 内部用 R(调试/验证用,实际方案中保密)
};

/**
 * 主调度函数
 *
 * 输入:
 *   A    — 公共参数 pp,N × m
 *   G    — gadget 矩阵, N × m
 *   mu_D — 消息(标量 ∈ {0, 1, …, q-1})
 *   p    — 参数集
 *   seed — 随机种子
 *
 * 输出: (𝒰, 𝒞)
 */
inline UniEncOutput uni_enc(const Mat& A, const Mat& G, long mu_D,
                            const Params& p, uint64_t seed = 0) {
    std::mt19937_64 rng(seed ? seed : std::random_device{}());

    // ─── 步骤 1: 𝒞 = A·R + μ_D · G ───
    // ① 采样 R ∈ {0,1}^{m×m}
    Mat R = sample_binary_mat(p.m, p.m, rng);
    // ② 矩阵乘 A·R  (N×m · m×m → N×m)
    Mat AR = mat_mul(A, R, p.q);
    // ③ 标量乘 μ_D·G
    Mat muG = scalar_mat_mul(mu_D, G, p.q);
    // ④ 加法
    Mat C = mat_add(AR, muG, p.q);

    // ─── 步骤 2: 加密 R 的每个比特 ───
    // ⑤+⑥ 逐元素加密 + 打包
    std::vector<Vec> U;
    U.reserve((size_t)p.m * p.m);
    for (int i = 0; i < p.m; i++) {
        for (int j = 0; j < p.m; j++) {
            U.push_back(lwe_encrypt_bit(A, R[i][j], p.q, p.sigma, rng));
        }
    }

    return UniEncOutput{ std::move(C), std::move(U), std::move(R) };
}

/* ══════════════════════════════════════════════════
   §7  LWE Secret / Noise Sampling (unified API)
   ══════════════════════════════════════════════════ */
/**
 * sample_lwe_secret — 生成 LWE 私钥向量 s ∈ D_{Z,σ}^m
 *
 * 仿照 OpenFHE 的统一密钥生成接口:
 *   s ← D_{Z,σ}^m  (离散高斯)
 *
 * OpenFHE 对应: DiscreteUniformGeneratorImpl -> ternary
 * 这里使用: mp12::DGSampler  (离散高斯, 12σ tail cut)
 *
 * 参数:
 *   m     — 向量维度
 *   sigma — 高斯宽度 (默认 3.2, 对应 LWE 标准参数)
 *   seed  — 随机种子 (0 = 随机设备)
 */
inline Vec sample_lwe_secret(int m, double sigma = 3.2, uint64_t seed = 0) {
    mp12::DGSampler dg(sigma, 0.0, seed);
    Vec s(m);
    for (int i = 0; i < m; ++i) s[i] = dg.sample();
    return s;
}

/**
 * sample_lwe_noise — 生成 LWE 噪声向量 e ∈ D_{Z,σ}^N
 *
 * 与 sample_lwe_secret 相同的底层采样器,
 * 语义分离是为了代码可读性和未来扩展 (如改用 CBD)。
 */
inline Vec sample_lwe_noise(int N, double sigma = 3.2, uint64_t seed = 0) {
    mp12::DGSampler dg(sigma, 0.0, seed);
    Vec e(N);
    for (int i = 0; i < N; ++i) e[i] = dg.sample();
    return e;
}

/**
 * sample_lwe_secret_ternary — 生成三元 LWE 私钥 s ∈ {−1,0,1}^m
 *
 * 大多数主流库 (OpenFHE, SEAL, Lattigo 等) 默认使用三元分布:
 *   优点: (1) 操作极快 (只需加减) (2) 安全性与高斯等价
 *
 * 使用 std::uniform_int_distribution<int>(0, 2) 然后映射:
 *   0 → -1,  1 → 0,  2 → 1
 */
inline Vec sample_lwe_secret_ternary(int m, uint64_t seed = 0) {
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::uniform_int_distribution<int> dist(0, 2);
    Vec s(m);
    for (int i = 0; i < m; ++i) {
        int v = dist(rng);
        s[i] = (v == 0) ? -1 : (v == 2 ? 1 : 0);
    }
    return s;
}

} // namespace unienc
