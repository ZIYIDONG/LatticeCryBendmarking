#pragma once
/**
 * decrypt.h — 多身份全同态加密阈值解密 (GSW 类型密文)
 *
 * 算法来源:
 *   ξ.Decrpty = (ξ.PartDec, ξ.FinDec)
 *
 * ──────────────────────────────────────────────────────────
 *  ξ.PartDec(ĉ, (id_1,…,id_N), k, sk_{id_k})
 * ──────────────────────────────────────────────────────────
 *   输入:
 *     ĉ   — 扩展密文 Ĉ = [Ĉ_1; …; Ĉ_N]^T,  维度 (N·R) × (N·M)
 *     k   — 当前解密参与者的身份索引 (0-based)
 *     t_k — 身份 id_k 的私钥,  长度 R = (d+1)·n + 1
 *
 *   步骤:
 *     ① 构造 ŵ = [0, …, 0, ⌈q/2⌉],长度 R
 *     ② 计算 G⁻¹(ŵ^T),长度 M = R·k (b 进制位分解)
 *     ③ 提取 Ĉ_k (第 k 个块行, R × (N·M))
 *     ④ 计算 v = t_k · Ĉ_k  (1 × (N·M))
 *     ⑤ 对每个列块 j, 累加 v_j · G⁻¹(ŵ^T) 得标量 γ_k
 *        等价于: γ_k = t_k · Ĉ_k · Ĝ⁻¹(ŵ_ext^T)
 *     ⑥ 采样掩蔽噪声 e^{sm} ← [-2^{dλlogλ}·B_χ, 2^{dλlogλ}·B_χ]
 *     ⑦ ED_k = γ_k + e^{sm}  (mod q)
 *
 * ──────────────────────────────────────────────────────────
 *  ξ.FinDec(ED_1, …, ED_N)
 * ──────────────────────────────────────────────────────────
 *   输入: N 个部分解密份额 ED_i
 *   ① p   = Σ ED_i  (mod q)
 *   ② μ_D = |⌈p / (q/2)⌉|  ∈ {0, 1}
 *
 * 维度约定 (与已有头文件一致):
 *   R = (d+1)·n + 1       (unienc.h 中的 N)
 *   k = ⌈log_b q⌉         (eval.h 的 eval_compute_k)
 *   M = R · k              (满足 m = R·k, eval.h 注释)
 *   G = build_gadget(R, q, b)     ∈ Z_q^{R × M}
 *   G⁻¹ = gadget_inverse          Z_q^{R × c} → Z^{M × c}
 */

#include "matops_plain-LWE.h"
#include "eval_plain-LWE.h"
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>
#include <cstdint>

namespace cryptolib {

/* 复用 matops / eval 的类型 */
using matops::Mat;
using matops::Vec;
using matops::make_mat;

/* ══════════════════════════════════════════════════════════
   §1  多身份解密参数
   ══════════════════════════════════════════════════════════ */
struct MIDParams {
    int  n;          // 格维度
    int  d;          // 电路深度
    long q;          // 模数
    int  b;          // gadget 基 (通常 2)
    int  k;          // ⌈log_b q⌉
    int  N_id;       // 身份数量 N
    int  R;          // 每个块的行数 = (d+1)·n + 1
    int  M;          // 每个块的列数 = R · k
    int  lambda;     // 安全参数 λ
    int  B_chi;      // 原始 LWE 噪声界 B_χ

    static MIDParams make(int n, int d, long q, int N_id,
                          int b = 2, int lambda = 8, int B_chi = 1)
    {
        MIDParams p;
        p.n      = n;
        p.d      = d;
        p.q      = q;
        p.b      = b;
        p.k      = eval_compute_k(q, b);
        p.N_id   = N_id;
        p.R      = (d + 1) * n + 1;
        p.M      = p.R * p.k;      // m = R·k
        p.lambda = lambda;
        p.B_chi  = B_chi;
        return p;
    }
};

/* ══════════════════════════════════════════════════════════
   §2  辅助函数
   ══════════════════════════════════════════════════════════ */

/**
 * 掩蔽噪声界: B_sm = 2^{d·λ·log₂λ} · B_χ
 * 需要足够大以统计隐藏 γ_k 中的小误差
 */
inline long smudging_bound(const MIDParams& p) {
    double exponent = (double)p.d * p.lambda
                      * std::log2(std::max(2, p.lambda));
    /* 防止指数爆炸——演示用时 clamp 到可控范围 */
    if (exponent > 40.0) exponent = 40.0;
    long bound = (long)std::pow(2.0, exponent) * p.B_chi;
    if (bound < 1) bound = 1;
    /* 保守限制: 避免掩蔽噪声接近判决阈值 q/4，从而破坏解密。
       将上限设置为 q/8（更保守），确保掩蔽噪声远小于阈值。 */
    if (bound > p.q / 8) bound = p.q / 8;
    return bound;
}

/**
 * 对称区间均匀采样: x ← [-bound, bound]
 */
inline long sample_symmetric(long bound, std::mt19937_64& rng) {
    if (bound <= 0) return 0;
    std::uniform_int_distribution<long> dist(-bound, bound);
    return dist(rng);
}

/**
 * 构造 ŵ = [0, …, 0, ⌈q/2⌉], 长度 R
 *
 * 这是 GSW 解密的核心选择向量:
 *   t · C · G⁻¹(ŵ^T) ≈ μ · t · ŵ^T = μ · t[R-1] · ⌈q/2⌉
 */
inline Vec build_w_hat(int R, long q) {
    Vec w(R, 0);
    w[R - 1] = (q + 1) / 2;     // ⌈q/2⌉
    return w;
}

/**
 * 提取扩展密文 Ĉ 的第 k 个块行: Ĉ_k ∈ Z_q^{R × (N·M)}
 * Ĉ 总维度: (N·R) × (N·M)
 */
inline Mat extract_block_row(const Mat& C_hat, int k, int R) {
    size_t NM = C_hat[0].size();
    Mat Ck(R, Vec(NM, 0));
    for (int r = 0; r < R; ++r)
        Ck[r] = C_hat[(size_t)k * R + r];
    return Ck;
}

/**
 * 提取第 (a,b) 个 R×M 子块
 */
inline Mat extract_block(const Mat& C_hat, int a, int b_idx,
                         int R, int M) {
    Mat blk(R, Vec(M, 0));
    for (int r = 0; r < R; ++r)
        for (int c = 0; c < M; ++c)
            blk[r][c] = C_hat[(size_t)a * R + r][(size_t)b_idx * M + c];
    return blk;
}

/**
 * 向量-矩阵乘: v (1×R) · A (R×C) → result (1×C),  mod q
 * 用 matops 风格实现
 */
inline Vec vec_mat_mul(const Vec& v, const Mat& A, long q) {
    int R = (int)v.size();
    int C = (int)A[0].size();
    Vec result(C, 0);
    for (int r = 0; r < R; ++r) {
        long vr = v[r];
        if (vr == 0) continue;
        const Vec& row = A[r];
        for (int c = 0; c < C; ++c)
            result[c] = matops::mod_pos(result[c] + vr * row[c], q);
    }
    return result;
}

/**
 * 向量内积 mod q
 */
inline long vec_dot_mod(const Vec& a, const Vec& b, long q) {
    long acc = 0;
    for (size_t i = 0; i < a.size(); ++i)
        acc = matops::mod_pos(acc + a[i] * b[i], q);
    return acc;
}

/**
 * 中心化归约: 将 x ∈ [0, q) 映射到 (-q/2, q/2]
 */
inline long center_mod(long x, long q) {
    long r = matops::mod_pos(x, q);
    if (r > q / 2) r -= q;
    return r;
}

/* ══════════════════════════════════════════════════════════
   §3  ξ.PartDec — 部分解密
   ══════════════════════════════════════════════════════════ */
/**
 * part_dec(C_hat, k, t_k, params, seed) → ED_k
 *
 * 输入:
 *   C_hat  — 扩展密文 (N·R × N·M)
 *   k      — 当前身份索引 (0-based)
 *   t_k    — 私钥向量, 长度 R
 *   params — 参数
 *   seed   — 随机种子 (掩蔽噪声用)
 *
 * 公式:
 *   ŵ = [0, …, 0, ⌈q/2⌉]            (长度 R)
 *   u = G⁻¹(ŵ^T)                     (长度 M, 列向量)
 *   γ_k = Σ_{j=0}^{N-1} (t_k · Ĉ_{k,j}) · u   (mod q)
 *       = t_k · Ĉ_k · Ĝ⁻¹(ŵ_ext^T)
 *   ED_k = γ_k + e^{sm}              (mod q)
 *
 * 其中 Ĝ⁻¹(ŵ_ext^T) 是将 G⁻¹(ŵ^T) 在每个列块重复 N 次,
 * 等价于 Ĝ = I_N ⊗ G 的逆在 ŵ_ext = (ŵ, ŵ, …, ŵ)^T 上的作用.
 */
inline long part_dec(const Mat& C_hat,
                     int k,
                     const Vec& t_k,
                     const MIDParams& params,
                     uint64_t seed = 0)
{
    if ((int)t_k.size() != params.R)
        throw std::invalid_argument("part_dec: |t_k| != R");
    if ((int)C_hat.size() != params.N_id * params.R)
        throw std::invalid_argument("part_dec: C_hat row count mismatch");
    if ((int)C_hat[0].size() != params.N_id * params.M)
        throw std::invalid_argument("part_dec: C_hat col count mismatch");

    const long q = params.q;

    /* ① 构造 ŵ = [0,…,0, ⌈q/2⌉], 长度 R */
    Vec w_hat = build_w_hat(params.R, q);

    /* ② G⁻¹(ŵ^T): 将 ŵ 视为 R×1 矩阵, 调用 gadget_inverse */
    Mat w_col(params.R, Vec(1));
    for (int i = 0; i < params.R; ++i)
        w_col[i][0] = w_hat[i];
    Mat u_col = gadget_inverse(w_col, q, params.b);
    /* u_col: M × 1, 提取为长度 M 的 Vec */
    Vec u(params.M);
    for (int i = 0; i < params.M; ++i)
        u[i] = u_col[i][0];

    /* ③ 提取第 k 个块行: Ĉ_k ∈ Z_q^{R × (N·M)} */
    /* ④ 计算 v = t_k · Ĉ_k  (1 × (N·M)) */
    /* ⑤ 逐列块计算 γ_k = Σ_j v_j · u */
    /*
     *  优化: 不需要构造完整的 Ĉ_k, 逐块操作节省内存
     *  对每个列块 j:
     *    Ĉ_{k,j} ∈ Z_q^{R × M}
     *    v_j = t_k · Ĉ_{k,j}  (1 × M)
     *    γ_k += v_j · u        (标量)
     */
    long gamma_k = 0;
    for (int j = 0; j < params.N_id; ++j) {
        /* 提取 (k, j) 块: R × M */
        Mat Ckj = extract_block(C_hat, k, j, params.R, params.M);

        /* v_j = t_k · Ckj  (1 × M) */
        Vec vj = vec_mat_mul(t_k, Ckj, q);

        /* γ_k += v_j · u  (mod q) */
        long dot = vec_dot_mod(vj, u, q);
        gamma_k = matops::mod_pos(gamma_k + dot, q);
    }

    /* ⑥ 采样掩蔽噪声 e^{sm} */
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    long B_sm = smudging_bound(params);
    long e_sm = sample_symmetric(B_sm, rng);

    /* ⑦ ED_k = γ_k + e^{sm}  (mod q) */
    long ED_k = matops::mod_pos(gamma_k + e_sm, q);

    return ED_k;
}

/* ══════════════════════════════════════════════════════════
   §4  ξ.FinDec — 最终解密
   ══════════════════════════════════════════════════════════ */
/**
 * fin_dec(ED, q) → μ_D ∈ {0, 1}
 *
 * 输入:
 *   ED — N 个部分解密份额
 *   q  — 模数
 *
 * 公式:
 *   p   = Σ_{i=1}^{N} ED_i                 (mod q)
 *   μ_D = |⌈ p / (q/2) ⌉|  ∈ {0, 1}
 *
 * 等价判定:
 *   将 p 中心化到 (-q/2, q/2]
 *   若 |p| ≤ q/4  →  μ = 0  (p 更接近 0)
 *   若 |p| >  q/4  →  μ = 1  (p 更接近 ±q/2)
 */
inline int fin_dec(const std::vector<long>& ED, long q) {
    /* ① 累加所有部分解密份额 */
    long p_sum = 0;
    for (long ed : ED)
        p_sum = matops::mod_pos(p_sum + ed, q);

    /* ② 中心化并判定 */
    long p_centered = center_mod(p_sum, q);
    long abs_p = (p_centered >= 0) ? p_centered : -p_centered;

    /* 阈值 = q/4 */
    long threshold = q / 4;

    return (abs_p > threshold) ? 1 : 0;
}

/* ══════════════════════════════════════════════════════════
   §5  辅助: 完整解密接口 (串联 PartDec + FinDec)
   ══════════════════════════════════════════════════════════ */
/**
 * decrypt_full(C_hat, keys, params) → μ_D
 *
 * 便捷接口: 对所有 N 个身份依次执行 PartDec, 然后调用 FinDec
 *
 * keys[i] = t_i, 长度 R 的私钥向量
 * 要求: Σ t_i[R-1] ≡ 1  (mod q, 近似) 以保证正确解密
 */
inline int decrypt_full(const Mat& C_hat,
                        const std::vector<Vec>& keys,
                        const MIDParams& params,
                        uint64_t base_seed = 0)
{
    if ((int)keys.size() != params.N_id)
        throw std::invalid_argument("decrypt_full: |keys| != N_id");

    std::vector<long> ED(params.N_id);

    for (int i = 0; i < params.N_id; ++i) {
        uint64_t seed_i = base_seed ? (base_seed + 100 * i) : 0;
        ED[i] = part_dec(C_hat, i, keys[i], params, seed_i);
    }

    return fin_dec(ED, params.q);
}

/* ══════════════════════════════════════════════════════════
   §6  调试辅助: 打印解密中间值
   ══════════════════════════════════════════════════════════ */
struct DecryptTrace {
    std::vector<long> gamma;       // 各身份的 γ_k (无噪声)
    std::vector<long> ED;          // 各身份的 ED_k (含掩蔽噪声)
    long p_sum;                    // Σ ED_i mod q
    long p_centered;               // 中心化后的 p
    int  mu;                       // 解密结果
};

/**
 * 带追踪的完整解密: 返回所有中间值供调试
 */
inline DecryptTrace decrypt_with_trace(const Mat& C_hat,
                                        const std::vector<Vec>& keys,
                                        const MIDParams& params,
                                        uint64_t base_seed = 42)
{
    DecryptTrace trace;
    trace.gamma.resize(params.N_id);
    trace.ED.resize(params.N_id);

    const long q = params.q;

    /* ŵ 和 G⁻¹(ŵ^T) — 所有身份共用 */
    Vec w_hat = build_w_hat(params.R, q);
    Mat w_col(params.R, Vec(1));
    for (int i = 0; i < params.R; ++i)
        w_col[i][0] = w_hat[i];
    Mat u_col = gadget_inverse(w_col, q, params.b);
    Vec u(params.M);
    for (int i = 0; i < params.M; ++i)
        u[i] = u_col[i][0];

    std::mt19937_64 rng(base_seed);
    long B_sm = smudging_bound(params);

    for (int k = 0; k < params.N_id; ++k) {
        /* 计算 γ_k */
        long gamma_k = 0;
        for (int j = 0; j < params.N_id; ++j) {
            Mat Ckj = extract_block(C_hat, k, j, params.R, params.M);
            Vec vj  = vec_mat_mul(keys[k], Ckj, q);
            gamma_k = matops::mod_pos(gamma_k + vec_dot_mod(vj, u, q), q);
        }
        trace.gamma[k] = gamma_k;

        /* 掩蔽噪声 */
        long e_sm = sample_symmetric(B_sm, rng);
        trace.ED[k] = matops::mod_pos(gamma_k + e_sm, q);
    }

    /* FinDec */
    trace.p_sum = 0;
    for (int i = 0; i < params.N_id; ++i)
        trace.p_sum = matops::mod_pos(trace.p_sum + trace.ED[i], q);

    trace.p_centered = center_mod(trace.p_sum, q);
    long abs_p = (trace.p_centered >= 0) ? trace.p_centered : -trace.p_centered;
    trace.mu = (abs_p > q / 4) ? 1 : 0;

    return trace;
}

} // namespace cryptolib

