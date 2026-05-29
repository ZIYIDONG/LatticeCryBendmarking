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
 * 编译: 由顶层 CMakeLists.txt 管理，通过 make 构建
 */

#include "matops_plain-LWE.h"
#include "eval_plain-LWE.h"
#include "decrypt_plain-LWE.h"
#include "unified_params_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cassert>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace matops;

/* ───── 文件输出辅助 ───── */
#include "bench_utils_plain-LWE.h"
static std::ostringstream dec_oss;

using namespace cryptolib;
#include "lwe_test_helpers_plain-LWE.h"

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
void run_test_decrypt() {
    std::cout << "==========================================================\n";
    std::cout << "  多身份全同态加密 — 阈值解密测试 (GSW 类型)\n";
    std::cout << "  使用: matops.h, eval.h (gadget_inverse, build_gadget)\n";
    std::cout << "==========================================================\n\n";

    /* ─── 参数设置 (统一 128-bit) ─── */
    const int d = 1;       // 深度
    const int N_id = 3;    // 身份数量
    auto params = unified::default_midparams_128(d, N_id);
    const long q = params.q;
    const int b = params.b;
    const int B_chi = params.B_chi;

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
    std::cout << "[LWE] 正在构造 A ∈ Z_q^{" << params.R << " × " << params.M << "} ... " << std::flush;
    Mat A = generate_lwe_matrix(t_master, params.R, params.M, q, B_chi, rng);
    std::cout << "完成\n" << std::flush;
    std::cout << "[验证] ‖t · A‖_∞ ≤ " << B_chi << " ... " << std::flush;
    bool lwe_ok = verify_lwe(t_master, A, q, B_chi);
    std::cout << (lwe_ok ? "✓ 通过" : "✗ 失败") << "\n\n";

    /* ─── 构造 gadget G ─── */
    Mat G = build_gadget(params.R, q, b);
    std::cout << "[Gadget] G ∈ Z_q^{" << params.R << " × " << params.M << "} 已构造\n\n";

    /* ─── 对 μ = 0 和 μ = 1 分别测试 ─── */
    int total_tests = 0, passed_tests = 0;

    for (int mu_test = 0; mu_test <= 1; ++mu_test) {
        std::cout << "──────────────────────────────────────────\n";
        std::cout << "  加密明文 μ = " << mu_test << "\n";
        std::cout << "──────────────────────────────────────────\n";

        /* GSW 加密 (L1 下极慢) */
        std::cout << "  [Enc] 正在加密 ... " << std::flush;
        Mat C = gsw_encrypt(A, G, mu_test, q, rng);
        std::cout << "完成\n" << std::flush;
        std::cout << "  [Enc] C = A·S + " << mu_test << "·G  ∈ Z_q^{"
                  << params.R << " × " << params.M << "}\n";

        /* 验证单密钥解密正确性 */
        std::cout << "  [单密钥验证] 正在计算 ... " << std::flush;
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
        std::cout << "完成\n" << std::flush;
        std::cout << "  [结果] t·C·G⁻¹(ŵ^T) = " << single_dec
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

    int random_trials = 50;
    std::string hist_path;
    unsigned long long provided_seed = 0;

    std::cout << "  多轮随机测试 (" << random_trials << " 轮)\n";
    std::cout << "──────────────────────────────────────────\n";

    /* 可选: 使用提供的种子或随机设备初始化 RNG 以获得非确定性测试 */
    if (provided_seed != 0) rng.seed(provided_seed);
    else {
        std::random_device rd;
        rng.seed(rd());
    }

    std::ofstream hist_file;
    if (!hist_path.empty()) {
        hist_file.open(hist_path);
        if (hist_file.is_open()) hist_file << "trial,mu,p_centered,pass\n";
    }

    int extra_pass = 0;
    long long sum_abs_p = 0;
    long long max_abs_p = 0;

    std::cout << "  " << random_trials << " 轮随机测试: " << std::flush;
    for (int trial = 0; trial < random_trials; ++trial) {
        std::cout << "." << std::flush;
        int mu_trial = trial % 2;

        /* 每轮重新生成共享密钥 */
        auto trial_keys = generate_shared_keys(t_master, N_id, q, rng);

        /* 加密 + 扩展 */
        Mat C_trial = gsw_encrypt(A, G, mu_trial, q, rng);
        Mat C_hat_trial = simple_expand(C_trial, N_id);

        /* 使用带追踪的解密以获得 p/ED/gamma */
        DecryptTrace trace = decrypt_with_trace(C_hat_trial, trial_keys, params, (uint64_t)(1000 + trial));
        int result = trace.mu;

        total_tests++;
        if (result == mu_trial) { passed_tests++; extra_pass++; }

        long long p_cent = std::llabs(trace.p_centered);
        sum_abs_p += p_cent;
        if (p_cent > max_abs_p) max_abs_p = p_cent;

        if (hist_file.is_open()) {
            hist_file << trial << "," << mu_trial << "," << trace.p_centered << "," << (result==mu_trial ? 1 : 0) << "\n";
        }
    }

    double mean_abs_p = (random_trials > 0) ? ((double)sum_abs_p / random_trials) : 0.0;
    std::cout << "  " << random_trials << " 轮中通过: " << extra_pass << " / " << random_trials << "\n";
    std::cout << "  p (|centered|) mean = " << mean_abs_p << ", max = " << max_abs_p << "\n\n";

    std::cout << " 完成\n" << std::flush;
    if (hist_file.is_open()) hist_file.close();

    /* ─── 总结 ─── */
    std::cout << "==========================================================\n";
    std::cout << "  总测试: " << total_tests
              << ",  通过: " << passed_tests
              << ",  失败: " << (total_tests - passed_tests) << "\n";
    std::cout << "  " << (passed_tests == total_tests ? "✓ 全部通过!" : "✗ 有失败") << "\n";
    std::cout << "==========================================================\n";

    /* 写入文件 */
    dec_oss << "\n=== Test: Multi-ID Threshold Decryption ===\n"
            << "  Total tests: " << total_tests
            << ", Passed: " << passed_tests
            << ", Failed: " << (total_tests - passed_tests) << "\n"
            << "  Result: " << (passed_tests == total_tests ? "ALL PASS" : "SOME FAIL") << "\n"
            << "  p (|centered|) mean = " << std::fixed << std::setprecision(0) << mean_abs_p
            << ", max = " << max_abs_p << "\n";
    bench_write(dec_oss.str());

}


