#pragma once
/**
 * MP12 DelTrapGen — Delegated Trapdoor Generation
 * Micciancio & Peikert, EUROCRYPT 2012, Section 5
 *
 * 三个核心算法：
 *
 * ① DelTrapGen (Tag-Based Delegation, MP12 §5)
 *   ─────────────────────────────────────────
 *   给定基础陷门 R（满足 A·[R;I] = G，tag = I_n），
 *   以及新 tag H ∈ GL_n(Z_q)（可逆方阵），
 *   生成满足 A·T_H = H·G (mod q) 的新陷门 T_H。
 *
 *   构造：T_H = [R·(H⊗I_k);  (H⊗I_k)]
 *   其中 H⊗I_k 是 Kronecker 积，维度 nk×nk。
 *
 *   验证：A·T_H = [Ā | G-ĀR]·[R(H⊗Ik); (H⊗Ik)]
 *             = ĀR(H⊗Ik) + (G-ĀR)(H⊗Ik)
 *             = G·(H⊗Ik) = H·G  ✓
 *   （利用恒等式 G·(H⊗I_k) = H·G，证明见下）
 *
 * ② SampleLeft (ABB10 + MP12)
 *   ─────────────────────────
 *   给定 A 的陷门 T_A 和任意矩阵 B，
 *   对扩展矩阵 [A|B] 采样短向量 x，满足 [A|B]·x = u。
 *
 *   策略：x2 ← 高斯随机，u' = u - B·x2，x1 ← SamplePre(A,T_A,u')
 *   输出 x = [x1; x2]
 *
 * ③ SampleRight (ABB10 + MP12)
 *   ──────────────────────────
 *   给定 A 和短矩阵 B_R（满足 B = A·B_R + G），
 *   对扩展矩阵 [A|B] 采样短向量 x，满足 [A|B]·x = u。
 *
 *   策略：扰动 p，G·x2 = u - A·p，x1 = p - B_R·x2
 *
 * 数学背景（Kronecker 积性质）：
 *   G = I_n ⊗ g^T  （n×nk）
 *   G·(H⊗I_k) = (I_n·H)⊗(g^T·I_k) = H⊗g^T = H·G  ✓
 */

#include "mp12.h"
#include <numeric>  // for gcd

namespace mp12 {

/* ══════════════════════════════════════════════════
   §1  数学基础：模运算辅助工具
   ══════════════════════════════════════════════════ */

/**
 * 扩展欧几里得算法：求 a 在模 q 意义下的逆元
 * 要求 gcd(a, q) = 1，否则抛出异常
 */
inline long mod_inv(long a, long q) {
    a = ((a % q) + q) % q;
    long t = 0, newt = 1;
    long r = q, newr = a;
    while (newr != 0) {
        long quotient = r / newr;
        long tmp;
        tmp = newt; newt = t - quotient * newt; t = tmp;
        tmp = newr; newr = r - quotient * newr; r = tmp;
    }
    if (r > 1) throw std::runtime_error("Element not invertible mod q");
    return (t % q + q) % q;
}

/**
 * n×n 方阵在 Z_q 上的求逆（Gauss-Jordan 消元）
 * 要求 det(H) 与 q 互素，否则不可逆
 */
inline Mat mat_inv_mod(const Mat& H, long q) {
    int n = H.size();
    // 构造增广矩阵 [H | I_n]
    Mat aug(n, Vec(2 * n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            aug[i][j] = mod(H[i][j], q);
        aug[i][n + i] = 1;
    }

    for (int col = 0; col < n; col++) {
        // 找主元（非零行）
        int pivot = -1;
        for (int row = col; row < n; row++) {
            if (aug[row][col] != 0) { pivot = row; break; }
        }
        if (pivot < 0) throw std::runtime_error("Matrix not invertible mod q");
        std::swap(aug[col], aug[pivot]);

        // 主元归一
        long inv_p = mod_inv(aug[col][col], q);
        for (int j = 0; j < 2 * n; j++)
            aug[col][j] = mod(aug[col][j] * inv_p, q);

        // 消去同列其他行
        for (int row = 0; row < n; row++) {
            if (row == col || aug[row][col] == 0) continue;
            long factor = aug[row][col];
            for (int j = 0; j < 2 * n; j++)
                aug[row][j] = mod(aug[row][j] - factor * aug[col][j], q);
        }
    }

    // 提取右半部分作为逆矩阵
    Mat inv = make_mat(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = aug[i][n + j];
    return inv;
}

/**
 * 矩阵向量乘（整数，不取模）——用于陷门合并步骤
 * 结果可能很大，调用者自行处理溢出
 */
inline Vec mat_vec_int(const Mat& A, const Vec& x) {
    int r = A.size(), c = A[0].size();
    Vec b(r, 0);
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            b[i] += A[i][j] * x[j];
    return b;
}

/**
 * 整数矩阵乘（不取模）：C = A·B
 */
inline Mat mat_mul_int(const Mat& A, const Mat& B) {
    int r = A.size(), mid = B.size(), c = B[0].size();
    Mat C = make_mat(r, c, 0);
    for (int i = 0; i < r; i++)
        for (int k = 0; k < mid; k++) {
            if (A[i][k] == 0) continue;
            for (int j = 0; j < c; j++)
                C[i][j] += A[i][k] * B[k][j];
        }
    return C;
}

/**
 * Kronecker 积 H ⊗ I_k
 * H: n×n 整数矩阵，I_k: k×k 单位矩阵
 * 结果: nk×nk 块对角矩阵，块(i,j) = H[i][j]·I_k
 *
 * 即：(H⊗I_k)[i*k+a][j*k+b] = H[i][j] · δ_{a,b}
 */
inline Mat kron_H_Ik(const Mat& H, int k) {
    int n = (int)H.size();
    int nk = n * k;
    Mat result = make_mat(nk, nk, 0);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (H[i][j] != 0)
                for (int a = 0; a < k; a++)
                    result[i * k + a][j * k + a] = H[i][j];
    return result;
}

/**
 * 生成随机可逆矩阵 H ∈ GL_n(Z_q)
 * 方法：对角矩阵，每个对角元素从 Z_q* 中随机选取（与 q 互素）
 * 对一般应用已足够；若需完全随机可逆矩阵，可换用 LU 随机生成
 */
inline Mat random_invertible_mat(int n, long q, uint64_t seed = 0) {
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    std::uniform_int_distribution<long> dist(1, q - 1);

    Mat H = make_mat(n, n, 0);
    for (int i = 0; i < n; i++) {
        long d;
        // 找与 q 互素的随机元素
        do { d = dist(rng); } while (std::gcd(d, q) != 1);
        H[i][i] = d;
    }
    return H;
}

/* ══════════════════════════════════════════════════
   §2  Tagged Trapdoor 结构体
   ══════════════════════════════════════════════════ */

/**
 * 带 tag 的陷门结构
 *
 * 陷门关系：A · T_H = H · G  (mod q)
 * 其中 T_H = [R_H; H⊗I_k]，R_H = R_base · (H⊗I_k)
 *
 * 当 tag = I_n 时退化为基础 GenTrap 的陷门关系 A·[R;I] = G
 */
struct TaggedTrapdoor {
    Mat R_H;   // 上半部分陷门：m_bar × nk，= R_base · (H⊗I_k)
    Mat H;     // 当前 tag：n × n，∈ GL_n(Z_q)
    Mat A;     // 公共矩阵（与基础陷门共用，不随 tag 改变）
};

/* ══════════════════════════════════════════════════
   §3  DelTrapGen：Tag 委托算法（MP12 §5）
   ══════════════════════════════════════════════════ */

/**
 * DelTrapGen(td_base, H_new, q)
 *
 * 输入：
 *   td_base  — 基础陷门（GenTrap 输出，tag = I_n）
 *   H_new    — 新 tag，n×n 可逆矩阵 ∈ GL_n(Z_q)
 *   q        — 模数
 *
 * 输出：
 *   TaggedTrapdoor，满足 A · T_{H_new} = H_new · G  (mod q)
 *
 * 构造：
 *   H̃ = H_new ⊗ I_k          （nk×nk Kronecker 积）
 *   R_{H_new} = R_base · H̃    （m_bar×nk，整数乘法，不取模）
 *   T_{H_new} = [R_{H_new}; H̃]
 *
 * 数学验证：
 *   A · T_{H_new}
 *   = [Ā | G-ĀR] · [R·H̃; H̃]
 *   = Ā·R·H̃ + (G-ĀR)·H̃
 *   = G·H̃
 *   = G·(H_new⊗I_k)
 *   = H_new·G   （利用 G·(H⊗I_k) = H·G）  ✓
 */
inline TaggedTrapdoor del_trap_gen(const Params& p,
                                   const Trapdoor& td_base,
                                   const Mat& H_new) {
    // Step 1: 计算 Kronecker 积 H̃ = H_new ⊗ I_k  (nk×nk)
    Mat H_tilde = kron_H_Ik(H_new, p.k);

    // Step 2: R_{H_new} = R_base · H̃  (m_bar×nk，整数矩阵乘)
    Mat R_H = mat_mul_int(td_base.R, H_tilde);

    // A 保持不变（与 tag 无关）
    return TaggedTrapdoor{ std::move(R_H), H_new, td_base.A };
}

/* ══════════════════════════════════════════════════
   §4  SamplePre（带 tag 的原像采样）
   ══════════════════════════════════════════════════ */

/**
 * sample_pre_tagged(p, td_tagged, u, seed)
 *
 * 使用带 tag H 的陷门，采样短向量 x ∈ Z^m，满足 A·x = u (mod q)
 *
 * 与基础 SamplePre 的区别：
 *   基础版：G·w = (u-v)              →  A·T·w = G·w
 *   Tagged：G·w = H^{-1}·(u-v)      →  A·T_H·w = H·G·w = H·H^{-1}·(u-v) = (u-v) ✓
 *
 * 步骤：
 *   ① p   ← D_{Z^m, s}              （高斯扰动）
 *   ② v   = A·p  mod q
 *   ③ u'  = H^{-1}·(u - v)  mod q   （用 tag 逆调整目标）
 *   ④ w   ← SampleG(u')              （G·w = u'，b 进制分解）
 *   ⑤ x   = p + [R_H; H̃]·w          （合并）
 *
 *   验证：A·x = v + A·T_H·w = v + H·G·w = v + H·u' = v + (u-v) = u  ✓
 */
inline Vec sample_pre_tagged(const Params& p,
                              const TaggedTrapdoor& td,
                              const Vec& u,
                              uint64_t seed = 0) {
    int nk = p.n * p.k;
    DGSampler dsampler(p.s, 0.0, seed);

    // ① 高斯扰动
    Vec pert = dsampler.sample_vec(p.m);

    // ② v = A·p mod q
    Vec v = mat_vec_mod(td.A, pert, p.q);

    // ③ u - v mod q
    Vec u_minus_v(p.n);
    for (int i = 0; i < p.n; i++)
        u_minus_v[i] = mod(u[i] - v[i], p.q);

    // ③ H^{-1}·(u-v) mod q
    Mat H_inv = mat_inv_mod(td.H, p.q);
    Vec adj = mat_vec_mod(H_inv, u_minus_v, p.q);

    // ④ w ← SampleG(adj)，满足 G·w = adj (mod q)
    Vec w = sample_g(p, adj);

    // ⑤ T_H·w = [R_H·w; H̃·w]
    //    H̃ = H ⊗ I_k（nk×nk）
    Mat H_tilde = kron_H_Ik(td.H, p.k);

    Vec T_w(p.m, 0);
    // 上半：R_H·w  (m_bar 行)
    for (int i = 0; i < p.m_bar; i++)
        for (int j = 0; j < nk; j++)
            T_w[i] += td.R_H[i][j] * w[j];
    // 下半：H̃·w  (nk 行)
    for (int i = 0; i < nk; i++)
        for (int j = 0; j < nk; j++)
            T_w[p.m_bar + i] += H_tilde[i][j] * w[j];

    // x = p + T_H·w
    Vec x(p.m);
    for (int i = 0; i < p.m; i++) x[i] = pert[i] + T_w[i];
    return x;
}

/* ══════════════════════════════════════════════════
   §5  SampleLeft：左扩展原像采样（ABB10 + MP12）
   ══════════════════════════════════════════════════ */

/**
 * sample_left(p, td_A, B, u, seed)
 *
 * 给定矩阵 A（有陷门 td_A）和任意矩阵 B ∈ Z_q^{n×m_B}，
 * 采样短向量 x = [x1; x2] ∈ Z^{m+m_B}，满足：
 *   [A | B] · [x1; x2] = u  (mod q)
 *
 * 应用场景（HIBE / ABE）：
 *   父节点 A 有陷门，子节点公钥含 B，
 *   可用父陷门为子节点生成密钥而无需知道 B 的结构。
 *
 * 算法：
 *   ① x2  ← D_{Z^{m_B}, s}            （对 B-部分随机高斯采样）
 *   ② u'  = u - B·x2  mod q            （消去 B 的贡献）
 *   ③ x1  ← SamplePre(A, td_A, u')    （用 A 的陷门解方程）
 *   输出 [x1; x2]
 *
 * 验证：A·x1 + B·x2 = u' + B·x2 = (u-B·x2) + B·x2 = u  ✓
 *
 * 短小性：‖x1‖ 由 SamplePre 控制（宽度 s），
 *         ‖x2‖ 由高斯采样控制（宽度 s）
 */
inline Vec sample_left(const Params& p,
                       const Trapdoor& td_A,
                       const Mat& B,
                       const Vec& u,
                       uint64_t seed = 0) {
    int m_B = (int)B[0].size();
    DGSampler dsampler(p.s, 0.0, seed);

    // ① 随机采样 x2 ∈ Z^{m_B}
    Vec x2 = dsampler.sample_vec(m_B);

    // ② u' = u - B·x2 mod q
    Vec Bx2 = mat_vec_mod(B, x2, p.q);
    Vec u_prime(p.n);
    for (int i = 0; i < p.n; i++)
        u_prime[i] = mod(u[i] - Bx2[i], p.q);

    // ③ x1 ← SamplePre(A, td_A, u')
    Vec x1 = sample_pre(p, td_A, u_prime, seed + 1);

    // 拼接 [x1; x2]
    Vec x(p.m + m_B);
    for (int i = 0; i < p.m; i++)      x[i] = x1[i];
    for (int i = 0; i < m_B; i++)  x[p.m + i] = x2[i];
    return x;
}

/* ══════════════════════════════════════════════════
   §6  SampleRight：右扩展原像采样（ABB10 + MP12）
   ══════════════════════════════════════════════════ */

/**
 * sample_right(p, A, B_R, u, seed)
 *
 * 给定矩阵 A ∈ Z_q^{n×m}，短矩阵 B_R ∈ Z^{m×nk}（"右陷门"），
 * B = A·B_R + G  (mod q)
 *
 * 采样短向量 x = [x1; x2] ∈ Z^{m+nk}，满足：
 *   [A | B] · [x1; x2] = u  (mod q)
 *
 * 应用场景（IBE 安全证明）：
 *   模拟器可构造 B = A·R + G，利用右陷门 R 生成合法密钥，
 *   而外部看来 B 是随机的（LWE 假设）。
 *
 * 算法：
 *   ① p   ← D_{Z^m, s}               （A 部分的扰动）
 *   ② v   = A·p  mod q
 *   ③ u'  = u - v  mod q
 *   ④ x2  ← SampleG(u')              （G·x2 = u'，b 进制分解）
 *   ⑤ x1  = p - B_R·x2               （整数计算，补偿 B 的贡献）
 *   输出 [x1; x2]
 *
 * 验证：
 *   A·x1 + B·x2
 *   = A·(p - B_R·x2) + (A·B_R + G)·x2
 *   = A·p - A·B_R·x2 + A·B_R·x2 + G·x2
 *   = A·p + G·x2
 *   = v + u'  = u  ✓
 *
 * 注：‖x1‖ ≈ ‖p‖ + ‖B_R‖·‖x2‖，需保证 s 足够大
 */
inline Vec sample_right(const Params& p,
                        const Mat& A,
                        const Mat& B_R,     // m × nk 短矩阵
                        const Vec& u,
                        uint64_t seed = 0) {
    int nk = p.n * p.k;
    DGSampler dsampler(p.s, 0.0, seed);

    // ① 扰动向量 p ← D_{Z^m, s}
    Vec pert = dsampler.sample_vec(p.m);

    // ② v = A·p mod q
    Vec v = mat_vec_mod(A, pert, p.q);

    // ③ u' = u - v mod q
    Vec u_prime(p.n);
    for (int i = 0; i < p.n; i++)
        u_prime[i] = mod(u[i] - v[i], p.q);

    // ④ x2 ← SampleG(u')，满足 G·x2 = u' (mod q)
    Vec x2 = sample_g(p, u_prime);

    // ⑤ B_R·x2  (整数，不取模)
    Vec BR_x2(p.m, 0);
    for (int i = 0; i < p.m; i++)
        for (int j = 0; j < nk; j++)
            BR_x2[i] += B_R[i][j] * x2[j];

    // x1 = p - B_R·x2
    Vec x1(p.m);
    for (int i = 0; i < p.m; i++) x1[i] = pert[i] - BR_x2[i];

    // 拼接 [x1; x2]
    Vec x(p.m + nk);
    for (int i = 0; i < p.m; i++)  x[i] = x1[i];
    for (int i = 0; i < nk; i++) x[p.m + i] = x2[i];
    return x;
}

/* ══════════════════════════════════════════════════
   §7  辅助验证函数
   ══════════════════════════════════════════════════ */

/**
 * 验证带 tag 的陷门关系：A · T_H = H · G  (mod q)
 * T_H = [R_H; H̃]，其中 H̃ = H ⊗ I_k
 */
inline bool verify_tagged_trapdoor(const Params& p,
                                    const TaggedTrapdoor& td) {
    int nk = p.n * p.k;
    Mat H_tilde = kron_H_Ik(td.H, p.k);

    // 构造 T_H（m × nk）
    Mat T_H = make_mat(p.m, nk);
    for (int i = 0; i < p.m_bar; i++)
        for (int j = 0; j < nk; j++)
            T_H[i][j] = td.R_H[i][j];
    for (int i = 0; i < nk; i++)
        for (int j = 0; j < nk; j++)
            T_H[p.m_bar + i][j] = H_tilde[i][j];

    // 计算 A·T_H mod q
    Mat AT_H = mat_mul_mod(td.A, T_H, p.q);

    // 计算 H·G mod q
    Mat G = gadget_matrix(p);
    Mat HG = mat_mul_mod(td.H, G, p.q);

    // 比较
    for (int i = 0; i < p.n; i++)
        for (int j = 0; j < nk; j++)
            if (AT_H[i][j] != HG[i][j]) return false;
    return true;
}

/**
 * 验证扩展矩阵的原像：[A|B]·x = u (mod q)
 */
inline bool verify_extended(const Mat& A, const Mat& B,
                              const Vec& x, const Vec& u, long q) {
    int n = (int)A.size();
    int mA = (int)A[0].size();
    int mB = (int)B[0].size();
    Vec Ax(n, 0), Bx(n, 0);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < mA; j++) Ax[i] = mod(Ax[i] + A[i][j] * x[j], q);
        for (int j = 0; j < mB; j++) Bx[i] = mod(Bx[i] + B[i][j] * x[mA + j], q);
    }
    for (int i = 0; i < n; i++)
        if (mod(Ax[i] + Bx[i], q) != mod(u[i], q)) return false;
    return true;
}

} // namespace mp12

// Entry point that runs all delegated trapdoor tests (defined at file scope in mp12deltrapgen.cpp)
void run_del_tests(const mp12::Params& p);

// Pure benchmarks for DelTrapGen and SamplePre_tagged (defined in mp12deltrapgen.cpp)
void bench_del_trap_gen(const mp12::Params& p);
void bench_sample_pre_tagged(const mp12::Params& p);
