/**
 * test_decrypt.cpp — 多身份全同态加密解密测试
 *
 * 测试流程:
 *   ① 选择参数, 构造 gadget G
 *   ② 生成 N 个私钥 t_i, 要求 Σ t_i[R-1] ≡ 1 (mod q)
 *   ③ 构造 GSW 密文 C = A·S + μ·G  (A 满足 t·A ≈ 0 的 LWE 结构)
 *   ④ 将 C 放入扩展密文 Ĉ 的对角位置 (简化的 expand)
 *   ⑤ 对每个身份执行 PartDec, 得到部分解密份额 ED_i
 *   ⑥ 执行 FinDec, 恢复明文 μ
 *   ⑦ 验证解密结果的正确性
 *
 * 编译:
 *   g++ -std=c++17 -O2 -Wall -I. test_decrypt.cpp -o test_decrypt
 */

#include "matops.h"
#include "eval.h"
#include "decrypt.h"
#include "unified_params.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>
#include <cmath>

using namespace matops;
using namespace cryptolib;

/* ═══════════════════════════════════════════════════════════
   §1  辅助: 构造满足 LWE 结构的公共矩阵 A
   ═══════════════════════════════════════════════════════════ */
/**
 * 生成 A ∈ Z_q^{R × M}, 使得对给定的秘密 t ∈ Z_q^R,
 * t · A 的每个分量都是小误差(≤ B)
 *
 * 构造方法:
 *   A 的前 (R-1) 行 ← U(Z_q)
 *   第 R-1 行 = (-(Σ_{i<R-1} t_i · A_i) + e) / t[R-1]  mod q
 *   其中 t[R-1] 取为与 q 互素的值 (默认 = 1)
 *
 * 这保证 t · A = e (小向量)
 */
static Mat generate_lwe_matrix(const Vec& t, int R, int M, long q,
                                int noise_bound, std::mt19937_64& rng)
{
    std::uniform_int_distribution<long> unif(0, q - 1);
    std::uniform_int_distribution<long> noise(-noise_bound, noise_bound);

    Mat A = make_mat(R, M);

    /* 前 R-1 行: 均匀随机 */
    for (int i = 0; i < R - 1; ++i)
        for (int j = 0; j < M; ++j)
            A[i][j] = unif(rng);

    /* 第 R-1 行: 使 t·A 的每个分量 = e_j (小噪声) */
    /* t[R-1] 应为 1 (或与 q 互素) */
    assert(t[R - 1] == 1 && "需要 t[R-1] = 1 以简化 LWE 构造");

    for (int j = 0; j < M; ++j) {
        long inner = 0;
        for (int i = 0; i < R - 1; ++i)
            inner = mod_pos(inner + t[i] * A[i][j], q);
        long e_j = noise(rng);
        /* t·A_col_j = inner + t[R-1]·A[R-1][j] = e_j */
        /* A[R-1][j] = (e_j - inner) / t[R-1] = e_j - inner  (因为 t[R-1]=1) */
        A[R - 1][j] = mod_pos(e_j - inner, q);
    }

    return A;
}

/* ═══════════════════════════════════════════════════════════
   §2  辅助: 构造 GSW 密文 C = A·S + μ·G
   ═══════════════════════════════════════════════════════════ */
/**
 * 标准 GSW 加密:
 *   S ← {0,1}^{M × M}         (短随机矩阵)
 *   C = A·S + μ·G  ∈ Z_q^{R × M}
 *
 * 解密正确性:
 *   t·C = t·A·S + μ·t·G
 *       ≈ (小误差)·S + μ·t·G
 *   t·C·G⁻¹(ŵ^T) ≈ μ·t·ŵ^T = μ·t[R-1]·⌈q/2⌉ ≈ μ·⌈q/2⌉
 */
static Mat gsw_encrypt(const Mat& A, const Mat& G, int mu,
                        long q, std::mt19937_64& rng)
{
    int R = (int)A.size();
    int M = (int)A[0].size();

    /* S ← {0,1}^{M × M} */
    std::uniform_int_distribution<int> bit(0, 1);
    Mat S = make_mat(M, M);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            S[i][j] = bit(rng);

    /* A · S  mod q */
    Mat AS = mat_mul(A, S, q);

    /* μ · G */
    Mat muG = make_mat(R, M, 0);
    if (mu != 0) {
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < M; ++j)
                muG[i][j] = mod_pos((long)mu * G[i][j], q);
    }

    /* C = AS + μG  mod q */
    Mat C = mat_add(AS, muG, q);
    return C;
}

/* ═══════════════════════════════════════════════════════════
   §3  辅助: 构造简化的扩展密文 (对角块 = C, 其余 = 0)
   ═══════════════════════════════════════════════════════════ */
/**
 * 简化的 expand: Ĉ_{a,a} = C (对角), Ĉ_{a,b} = 0 (非对角)
 *
 * 这对应 expand.h 中 i 行的非对角块为 0 的特例,
 * 适用于所有身份共享同一密文的基本场景.
 *
 * 维度: (N·R) × (N·M)
 */
static Mat simple_expand(const Mat& C, int N_id) {
    int R = (int)C.size();
    int M = (int)C[0].size();
    Mat C_hat(N_id * R, Vec(N_id * M, 0));

    for (int a = 0; a < N_id; ++a)
        for (int r = 0; r < R; ++r)
            for (int c = 0; c < M; ++c)
                C_hat[a * R + r][a * M + c] = C[r][c];

    return C_hat;
}

/* ═══════════════════════════════════════════════════════════
   §4  辅助: 生成加法秘密共享的私钥
   ═══════════════════════════════════════════════════════════ */
/**
 * 生成 N 个私钥 t_0, …, t_{N-1}, 满足 Σ t_i ≡ t_master (mod q)
 *
 * 策略:
 *   t_0, …, t_{N-2} ← 随机
 *   t_{N-1} = t_master - Σ_{i<N-1} t_i   (mod q)
 *
 * 要求: t_master[R-1] = 1, 从而 Σ t_i[R-1] ≡ 1 (mod q)
 */
static std::vector<Vec> generate_shared_keys(const Vec& t_master,
                                              int N_id, long q,
                                              std::mt19937_64& rng)
{
    int R = (int)t_master.size();
    std::uniform_int_distribution<long> unif(0, q - 1);
    std::vector<Vec> keys(N_id, Vec(R, 0));

    /* 前 N-1 个: 随机 */
    for (int i = 0; i < N_id - 1; ++i)
        for (int j = 0; j < R; ++j)
            keys[i][j] = unif(rng);

    /* 最后一个: 补齐使和 = t_master */
    for (int j = 0; j < R; ++j) {
        long partial_sum = 0;
        for (int i = 0; i < N_id - 1; ++i)
            partial_sum = mod_pos(partial_sum + keys[i][j], q);
        keys[N_id - 1][j] = mod_pos(t_master[j] - partial_sum, q);
    }

    return keys;
}

/* ═══════════════════════════════════════════════════════════
   §5  验证函数
   ═══════════════════════════════════════════════════════════ */

/** 验证 G · G⁻¹(X) ≡ X  (mod q) */
static bool verify_gadget(int R, long q, int b) {
    Mat G = build_gadget(R, q, b);
    int M = R * eval_compute_k(q, b);

    /* 测试向量 */
    Vec v = {(q + 1) / 2};
    Mat X(1, Vec(1, (q + 1) / 2));
    Mat Ginv = gadget_inverse(X, q, b);
    /* G · Ginv 应该得到 X */
    /* 这里 X 是 1×1, Ginv 是 k×1, G 是 1×k (当 R=1 时) */
    /* 一般情况: 用 R×1 测试 */
    Mat test_col(R, Vec(1, 0));
    test_col[R - 1][0] = (q + 1) / 2;
    Mat decomp = gadget_inverse(test_col, q, b);
    /* decomp: M × 1 */
    Mat reconstructed = mat_mul(G, decomp, q);
    /* reconstructed 应等于 test_col */
    for (int i = 0; i < R; ++i)
        if (mod_pos(reconstructed[i][0], q) != mod_pos(test_col[i][0], q))
            return false;
    return true;
}

/** 验证 t · A ≈ 0 (LWE 正确性) */
static bool verify_lwe(const Vec& t, const Mat& A, long q, int bound) {
    int M = (int)A[0].size();
    Vec tA = vec_mat_mul(t, A, q);
    for (int j = 0; j < M; ++j) {
        long c = center_mod(tA[j], q);
        if (c < -bound || c > bound) return false;
    }
    return true;
}

/* ═══════════════════════════════════════════════════════════
   §6  主测试
   ═══════════════════════════════════════════════════════════ */
int main() {
    std::cout << "==========================================================\n";
    std::cout << "  多身份全同态加密 — 阈值解密测试 (GSW 类型)\n";
    std::cout << "  使用: matops.h, eval.h (gadget_inverse, build_gadget)\n";
    std::cout << "==========================================================\n\n";

    /* ─── 参数设置 (统一 128-bit) ─── */
    const int d = 1;       // 深度
    const int N_id = 3;    // 身份数量
    auto params = unified::default_midparams_128(d, N_id);

    std::cout << "[参数]\n"
              << "  n     = " << params.n     << "\n"
              << "  d     = " << params.d     << "\n"
              << "  q     = " << params.q     << "\n"
              << "  b     = " << params.b     << "\n"
              << "  k     = " << params.k     << " (⌈log₂q⌉)\n"
              << "  R     = " << params.R     << " ((d+1)n+1)\n"
              << "  M     = " << params.M     << " (R·k)\n"
              << "  N_id  = " << params.N_id  << "\n"
              << "  B_sm  = " << smudging_bound(params) << " (掩蔽噪声界)\n"
              << "\n";

    /* ─── 验证 gadget 性质 ─── */
    std::cout << "[验证] G · G⁻¹(X) ≡ X  (mod q) ... ";
    bool gok = verify_gadget(params.R, q, b);
    std::cout << (gok ? "✓ 通过" : "✗ 失败") << "\n\n";

    /* ─── 生成主密钥 t_master ─── */
    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<long> unif(0, q - 1);

    Vec t_master(params.R);
    for (int i = 0; i < params.R - 1; ++i)
        t_master[i] = unif(rng);
    t_master[params.R - 1] = 1;           // 最后分量 = 1

    std::cout << "[主密钥] t_master = (";
    for (int i = 0; i < params.R; ++i) {
        if (i) std::cout << ", ";
        std::cout << t_master[i];
    }
    std::cout << ")\n";

    /* ─── 生成秘密共享的子密钥 ─── */
    auto keys = generate_shared_keys(t_master, N_id, q, rng);
    std::cout << "[密钥共享] " << N_id << " 个子密钥已生成\n";

    /* 验证: Σ t_i ≡ t_master (mod q) */
    Vec sum_keys(params.R, 0);
    for (int i = 0; i < N_id; ++i)
        for (int j = 0; j < params.R; ++j)
            sum_keys[j] = mod_pos(sum_keys[j] + keys[i][j], q);
    bool keys_ok = true;
    for (int j = 0; j < params.R; ++j)
        if (sum_keys[j] != mod_pos(t_master[j], q)) keys_ok = false;
    std::cout << "[验证] Σ t_i ≡ t_master  ... "
              << (keys_ok ? "✓ 通过" : "✗ 失败") << "\n\n";

    /* ─── 构造 LWE 矩阵 A ─── */
    Mat A = generate_lwe_matrix(t_master, params.R, params.M, q, B_chi, rng);
    std::cout << "[LWE] A ∈ Z_q^{" << params.R << " × " << params.M << "} 已生成\n";
    bool lwe_ok = verify_lwe(t_master, A, q, B_chi);
    std::cout << "[验证] ‖t · A‖_∞ ≤ " << B_chi << "  ... "
              << (lwe_ok ? "✓ 通过" : "✗ 失败") << "\n\n";

    /* ─── 构造 gadget G ─── */
    Mat G = build_gadget(params.R, q, b);
    std::cout << "[Gadget] G ∈ Z_q^{" << params.R << " × " << params.M << "} 已构造\n\n";

    /* ─── 对 μ = 0 和 μ = 1 分别测试 ─── */
    int total_tests = 0, passed_tests = 0;

    for (int mu_test = 0; mu_test <= 1; ++mu_test) {
        std::cout << "──────────────────────────────────────────\n";
        std::cout << "  加密明文 μ = " << mu_test << "\n";
        std::cout << "──────────────────────────────────────────\n";

        /* GSW 加密 */
        Mat C = gsw_encrypt(A, G, mu_test, q, rng);
        std::cout << "  [Enc] C = A·S + " << mu_test << "·G  ∈ Z_q^{"
                  << params.R << " × " << params.M << "}\n";

        /* 验证单密钥解密正确性: t · C · G⁻¹(ŵ^T) ≈ μ · ⌈q/2⌉ */
        Vec w_hat = build_w_hat(params.R, q);
        Mat w_col(params.R, Vec(1));
        for (int i = 0; i < params.R; ++i) w_col[i][0] = w_hat[i];
        Mat u_col = gadget_inverse(w_col, q, b);
        Vec u(params.M);
        for (int i = 0; i < params.M; ++i) u[i] = u_col[i][0];

        Vec tC = vec_mat_mul(t_master, C, q);
        long single_dec = vec_dot_mod(tC, u, q);
        long single_centered = center_mod(single_dec, q);
        long expected = mu_test * ((q + 1) / 2);
        std::cout << "  [单密钥验证] t·C·G⁻¹(ŵ^T) = " << single_dec
                  << " (中心化: " << single_centered << ")"
                  << ", 期望 ≈ " << (mu_test ? (long)((q + 1) / 2) : 0L) << "\n";

        /* 扩展密文 */
        Mat C_hat = simple_expand(C, N_id);
        std::cout << "  [Expand] Ĉ ∈ Z_q^{" << C_hat.size()
                  << " × " << C_hat[0].size() << "} (对角块 = C)\n";

        /* 带追踪的解密 */
        DecryptTrace trace = decrypt_with_trace(C_hat, keys, params, 42 + mu_test);

        std::cout << "  [PartDec]\n";
        for (int i = 0; i < N_id; ++i) {
            long gamma_c = center_mod(trace.gamma[i], q);
            long ed_c    = center_mod(trace.ED[i], q);
            std::cout << "    身份 " << i << ":  γ = " << std::setw(6) << gamma_c
                      << ",  ED = " << std::setw(6) << ed_c << "\n";
        }

        std::cout << "  [FinDec]\n"
                  << "    p = Σ ED_i = " << trace.p_sum
                  << " (中心化: " << trace.p_centered << ")\n"
                  << "    |p| = " << std::abs(trace.p_centered)
                  << ",  q/4 = " << q / 4 << "\n"
                  << "    μ_D = " << trace.mu << "\n";

        bool correct = (trace.mu == mu_test);
        std::cout << "  >> 解密结果: μ = " << trace.mu
                  << " (期望 " << mu_test << ")  "
                  << (correct ? "✓ 正确" : "✗ 错误") << "\n\n";

        total_tests++;
        if (correct) passed_tests++;
    }

    /* ─── 多轮随机测试 ─── */
    std::cout << "──────────────────────────────────────────\n";
    std::cout << "  多轮随机测试 (50 轮)\n";
    std::cout << "──────────────────────────────────────────\n";

    int extra_pass = 0;
    for (int trial = 0; trial < 50; ++trial) {
        int mu_trial = trial % 2;

        /* 每轮重新生成共享密钥 */
        auto trial_keys = generate_shared_keys(t_master, N_id, q, rng);

        /* 加密 + 扩展 */
        Mat C_trial = gsw_encrypt(A, G, mu_trial, q, rng);
        Mat C_hat_trial = simple_expand(C_trial, N_id);

        /* 解密 */
        int result = decrypt_full(C_hat_trial, trial_keys, params,
                                  (uint64_t)(1000 + trial));

        total_tests++;
        if (result == mu_trial) { passed_tests++; extra_pass++; }
    }
    std::cout << "  50 轮中通过: " << extra_pass << " / 50\n\n";

    /* ─── 总结 ─── */
    std::cout << "==========================================================\n";
    std::cout << "  总测试: " << total_tests
              << ",  通过: " << passed_tests
              << ",  失败: " << (total_tests - passed_tests) << "\n";
    std::cout << "  " << (passed_tests == total_tests ? "✓ 全部通过!" : "✗ 有失败") << "\n";
    std::cout << "==========================================================\n";

    return (passed_tests == total_tests) ? 0 : 1;
}

