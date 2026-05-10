/**
 * bench_decrypt.cpp — 多身份 FHE 解密 逐步耗时拆解
 *
 * 将 PartDec / FinDec 拆成最细粒度的基本操作,
 * 对每个操作多轮计时取平均, 输出:
 *   ① 单步绝对耗时 (μs)
 *   ② 占 PartDec 总时间的百分比
 *   ③ 随参数 (n, N_id, q) 变化的缩放趋势
 *
 * 编译:
 *   g++ -std=c++17 -O2 -Wall -I. bench_decrypt.cpp -o bench_decrypt
 */

#include "matops.h"
#include "eval.h"
#include "decrypt.h"
#include "unified_params.h"

#include <iostream>
#include <iomanip>
#include <chrono>
#include <random>
#include <vector>
#include <string>
#include <functional>
#include <cassert>
#include <cmath>
#include <numeric>

using namespace matops;
using namespace cryptolib;

/* ══════════════════════════════════════════════════
   计时工具
   ══════════════════════════════════════════════════ */
using Clock = std::chrono::high_resolution_clock;
using us_t  = std::chrono::microseconds;
using ns_t  = std::chrono::nanoseconds;

/**
 * 运行 func() 共 repeats 次, 返回总耗时 (μs)、
 * 单次平均 (μs) 和标准差 (μs)
 */
struct TimingResult {
    double total_us;
    double avg_us;
    double stddev_us;
    int    repeats;
};

static TimingResult bench(std::function<void()> func, int repeats) {
    std::vector<double> times(repeats);
    for (int r = 0; r < repeats; ++r) {
        auto t0 = Clock::now();
        func();
        auto t1 = Clock::now();
        times[r] = std::chrono::duration_cast<ns_t>(t1 - t0).count() / 1000.0;
    }
    double sum = std::accumulate(times.begin(), times.end(), 0.0);
    double avg = sum / repeats;
    double var = 0;
    for (double t : times) var += (t - avg) * (t - avg);
    var /= (repeats > 1) ? (repeats - 1) : 1;
    return { sum, avg, std::sqrt(var), repeats };
}

/* ══════════════════════════════════════════════════
   LWE 矩阵 / GSW 密文生成 (同 test_decrypt.cpp)
   ══════════════════════════════════════════════════ */
static Mat generate_lwe_matrix(const Vec& t, int R, int M, long q,
                                int noise_bound, std::mt19937_64& rng)
{
    std::uniform_int_distribution<long> unif(0, q - 1);
    std::uniform_int_distribution<long> noise(-noise_bound, noise_bound);
    Mat A = make_mat(R, M);
    for (int i = 0; i < R - 1; ++i)
        for (int j = 0; j < M; ++j)
            A[i][j] = unif(rng);
    for (int j = 0; j < M; ++j) {
        long inner = 0;
        for (int i = 0; i < R - 1; ++i)
            inner = mod_pos(inner + t[i] * A[i][j], q);
        A[R - 1][j] = mod_pos(noise(rng) - inner, q);
    }
    return A;
}

static Mat gsw_encrypt(const Mat& A, const Mat& G, int mu,
                        long q, std::mt19937_64& rng)
{
    int R = (int)A.size(), M = (int)A[0].size();
    std::uniform_int_distribution<int> bit(0, 1);
    Mat S = make_mat(M, M);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            S[i][j] = bit(rng);
    Mat AS = mat_mul(A, S, q);
    Mat muG = make_mat(R, M, 0);
    if (mu) for (int i = 0; i < R; ++i)
                for (int j = 0; j < M; ++j)
                    muG[i][j] = mod_pos((long)mu * G[i][j], q);
    return mat_add(AS, muG, q);
}

static Mat simple_expand(const Mat& C, int N_id) {
    int R = (int)C.size(), M = (int)C[0].size();
    Mat Ch(N_id * R, Vec(N_id * M, 0));
    for (int a = 0; a < N_id; ++a)
        for (int r = 0; r < R; ++r)
            for (int c = 0; c < M; ++c)
                Ch[a * R + r][a * M + c] = C[r][c];
    return Ch;
}

static std::vector<Vec> generate_shared_keys(const Vec& t_master,
                                              int N_id, long q,
                                              std::mt19937_64& rng)
{
    int R = (int)t_master.size();
    std::uniform_int_distribution<long> unif(0, q - 1);
    std::vector<Vec> keys(N_id, Vec(R, 0));
    for (int i = 0; i < N_id - 1; ++i)
        for (int j = 0; j < R; ++j)
            keys[i][j] = unif(rng);
    for (int j = 0; j < R; ++j) {
        long s = 0;
        for (int i = 0; i < N_id - 1; ++i)
            s = mod_pos(s + keys[i][j], q);
        keys[N_id - 1][j] = mod_pos(t_master[j] - s, q);
    }
    return keys;
}

/* ══════════════════════════════════════════════════
   格式化打印工具
   ══════════════════════════════════════════════════ */
struct StepEntry {
    std::string name;
    std::string formula;
    TimingResult timing;
};

static void print_table(const std::string& title,
                        const std::vector<StepEntry>& entries,
                        double total_partdec_us)
{
    std::cout << "\n┌─────────────────────────────────────────────"
                 "──────────────────────────────────────────────┐\n";
    std::cout << "│  " << std::left << std::setw(90) << title << "│\n";
    std::cout << "├────┬──────────────────────────┬─────────────────"
                 "────────┬────────────┬──────────┬───────┤\n";
    std::cout << "│ #  │ 操作名称                 │ 公式 / 维度"
                 "              │  平均 (μs) │ 标准差   │ 占比  │\n";
    std::cout << "├────┼──────────────────────────┼─────────────────"
                 "────────┼────────────┼──────────┼───────┤\n";

    for (size_t i = 0; i < entries.size(); ++i) {
        const auto& e = entries[i];
        double pct = (total_partdec_us > 0)
                     ? (e.timing.avg_us / total_partdec_us * 100.0) : 0;
        std::cout << "│ " << std::setw(2) << (i + 1) << " │ "
                  << std::left  << std::setw(24) << e.name << " │ "
                  << std::left  << std::setw(23) << e.formula << " │ "
                  << std::right << std::setw(10) << std::fixed
                  << std::setprecision(2) << e.timing.avg_us << " │ "
                  << std::setw(8) << std::setprecision(2) << e.timing.stddev_us << " │ "
                  << std::setw(5) << std::setprecision(1) << pct << "│\n";
    }
    std::cout << "└────┴──────────────────────────┴─────────────────"
                 "────────┴────────────┴──────────┴───────┘\n";
}

/* ══════════════════════════════════════════════════
   主函数: 逐步拆解 PartDec + FinDec 并计时
   ══════════════════════════════════════════════════ */
static void run_benchmark(int n, int d, long q, int N_id, int REPS)
{
    const int b     = 2;
    const int B_chi = 1;

    MIDParams params = MIDParams::make(n, d, q, N_id, b, /*lambda=*/4, B_chi);

    std::cout << "\n══════════════════════════════════════════════════\n";
    std::cout << "  参数: n=" << n << ", d=" << d << ", q=" << q
              << ", N=" << N_id
              << ", R=" << params.R << ", M=" << params.M
              << ", k=" << params.k
              << "\n  每步重复 " << REPS << " 次取平均\n";
    std::cout << "══════════════════════════════════════════════════\n";

    /* ─── 预处理: 生成密钥、密文 ─── */
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<long> unif(0, q - 1);

    Vec t_master(params.R);
    for (int i = 0; i < params.R - 1; ++i)
        t_master[i] = unif(rng);
    t_master[params.R - 1] = 1;

    auto keys = generate_shared_keys(t_master, N_id, q, rng);
    Mat A  = generate_lwe_matrix(t_master, params.R, params.M, q, B_chi, rng);
    Mat G  = build_gadget(params.R, q, b);
    Mat C  = gsw_encrypt(A, G, 1, q, rng);
    Mat Ch = simple_expand(C, N_id);

    const int k_idx = 0;                    // 被测的身份索引
    const Vec& t_k  = keys[k_idx];

    /* ─── 预计算一些中间结果供后续步骤使用 ─── */
    Vec w_hat;
    Mat w_col;
    Mat u_col;
    Vec u;
    Mat Ckj;
    Vec vj;
    long dot_val = 0;
    long gamma_k = 0;
    long e_sm    = 0;
    long ED_k    = 0;

    std::vector<StepEntry> steps;

    /* ================================================================
       PartDec 步骤拆解
       ================================================================ */

    /* ── Step 1: 构造 ŵ ── */
    steps.push_back({
        "构造 ŵ",
        "ŵ=[0,..,0,⌈q/2⌉] (R)",
        bench([&](){
            w_hat = build_w_hat(params.R, q);
        }, REPS)
    });

    /* ── Step 2: ŵ → 列矩阵 (R×1) ── */
    w_hat = build_w_hat(params.R, q);
    steps.push_back({
        "ŵ→列矩阵",
        "Vec→Mat (R×1)",
        bench([&](){
            w_col = Mat(params.R, Vec(1));
            for (int i = 0; i < params.R; ++i)
                w_col[i][0] = w_hat[i];
        }, REPS)
    });

    /* ── Step 3: G⁻¹(ŵ^T) — gadget_inverse ── */
    w_col = Mat(params.R, Vec(1));
    for (int i = 0; i < params.R; ++i) w_col[i][0] = w_hat[i];
    steps.push_back({
        "G⁻¹(ŵ^T)",
        "gadget_inv (R×1→M×1)",
        bench([&](){
            u_col = gadget_inverse(w_col, q, b);
        }, REPS)
    });

    /* ── Step 4: 列矩阵 → Vec ── */
    u_col = gadget_inverse(w_col, q, b);
    steps.push_back({
        "u_col→Vec",
        "Mat(M×1)→Vec(M)",
        bench([&](){
            u = Vec(params.M);
            for (int i = 0; i < params.M; ++i)
                u[i] = u_col[i][0];
        }, REPS)
    });
    u = Vec(params.M);
    for (int i = 0; i < params.M; ++i) u[i] = u_col[i][0];

    /* ── Step 5: 提取子块 Ĉ_{k,j} (单个块) ── */
    steps.push_back({
        "提取子块 Ĉ_{k,j}",
        "extract_block (R×M)",
        bench([&](){
            Ckj = extract_block(Ch, k_idx, 0, params.R, params.M);
        }, REPS)
    });

    /* ── Step 6: 向量-矩阵乘 t_k · Ĉ_{k,j} ── */
    Ckj = extract_block(Ch, k_idx, 0, params.R, params.M);
    steps.push_back({
        "向量-矩阵乘 t·Ĉ_{k,j}",
        "vec_mat_mul (1×R·R×M)",
        bench([&](){
            vj = vec_mat_mul(t_k, Ckj, q);
        }, REPS)
    });

    /* ── Step 7: 向量内积 v_j · u ── */
    vj = vec_mat_mul(t_k, Ckj, q);
    steps.push_back({
        "内积 v_j·u",
        "dot (M)",
        bench([&](){
            dot_val = vec_dot_mod(vj, u, q);
        }, REPS)
    });

    /* ── Step 8: γ_k 累加 (N个块的总循环) ── */
    steps.push_back({
        "γ_k 累加 (N块循环)",
        "Σ_j (t·Ĉ_{k,j})·u",
        bench([&](){
            gamma_k = 0;
            for (int j = 0; j < params.N_id; ++j) {
                Mat blk = extract_block(Ch, k_idx, j, params.R, params.M);
                Vec v   = vec_mat_mul(t_k, blk, q);
                gamma_k = mod_pos(gamma_k + vec_dot_mod(v, u, q), q);
            }
        }, REPS)
    });

    /* ── Step 9: 掩蔽噪声采样 ── */
    std::mt19937_64 noise_rng(123);
    long B_sm = smudging_bound(params);
    steps.push_back({
        "掩蔽噪声 e^{sm}",
        "uniform [-B_sm,B_sm]",
        bench([&](){
            e_sm = sample_symmetric(B_sm, noise_rng);
        }, REPS)
    });

    /* ── Step 10: ED_k = γ_k + e^{sm} ── */
    steps.push_back({
        "ED_k = γ_k + e^{sm}",
        "mod_pos (标量加法)",
        bench([&](){
            ED_k = mod_pos(gamma_k + e_sm, q);
        }, REPS)
    });

    /* ── 完整 PartDec (单个身份) ── */
    TimingResult full_partdec = bench([&](){
        part_dec(Ch, k_idx, t_k, params, 42);
    }, REPS);

    steps.push_back({
        "★ PartDec 完整 (单身份)",
        "steps 1–10 combined",
        full_partdec
    });

    print_table("ξ.PartDec 基本操作拆解  (身份 k=" + std::to_string(k_idx) + ")",
                steps, full_partdec.avg_us);

    /* ================================================================
       FinDec 步骤拆解
       ================================================================ */
    /* 先执行所有 PartDec 得到 ED 向量 */
    std::vector<long> ED_all(N_id);
    for (int i = 0; i < N_id; ++i)
        ED_all[i] = part_dec(Ch, i, keys[i], params, 100 + i);

    std::vector<StepEntry> fin_steps;

    /* ── FD-1: 累加 Σ ED_i ── */
    long p_sum = 0;
    fin_steps.push_back({
        "累加 Σ ED_i",
        "N 次 mod_pos 加法",
        bench([&](){
            p_sum = 0;
            for (int i = 0; i < N_id; ++i)
                p_sum = mod_pos(p_sum + ED_all[i], q);
        }, REPS)
    });

    /* ── FD-2: 中心化 mod ── */
    fin_steps.push_back({
        "中心化 center_mod",
        "(-q/2, q/2] 映射",
        bench([&](){
            (void)center_mod(p_sum, q);
        }, REPS)
    });

    /* ── FD-3: |p| 与 q/4 比较 ── */
    long p_c = center_mod(p_sum, q);
    fin_steps.push_back({
        "舍入判定 μ_D",
        "|p| vs q/4 → {0,1}",
        bench([&](){
            long ap = (p_c >= 0) ? p_c : -p_c;
            (void)(ap > q / 4 ? 1 : 0);
        }, REPS)
    });

    /* ── FinDec 完整 ── */
    TimingResult full_findec = bench([&](){
        fin_dec(ED_all, q);
    }, REPS);
    fin_steps.push_back({
        "★ FinDec 完整",
        "steps 1–3 combined",
        full_findec
    });

    print_table("ξ.FinDec 基本操作拆解", fin_steps, full_findec.avg_us);

    /* ================================================================
       全流程汇总 (N 个 PartDec + 1 个 FinDec)
       ================================================================ */
    TimingResult full_all = bench([&](){
        std::vector<long> ed(N_id);
        for (int i = 0; i < N_id; ++i)
            ed[i] = part_dec(Ch, i, keys[i], params, 200 + i);
        fin_dec(ed, q);
    }, REPS);

    std::cout << "\n┌─────────────────────────────────────────────────────┐\n";
    std::cout << "│  全流程汇总 (N=" << N_id << " × PartDec + FinDec)"
              << std::string(24 - std::to_string(N_id).size(), ' ') << "│\n";
    std::cout << "├──────────────────────────┬────────────┬──────────────┤\n";
    std::cout << "│ 阶段                     │  平均 (μs) │  占总时间    │\n";
    std::cout << "├──────────────────────────┼────────────┼──────────────┤\n";

    double all_partdec_us = full_partdec.avg_us * N_id;
    std::cout << "│ N × PartDec              │ "
              << std::setw(10) << std::fixed << std::setprecision(2)
              << all_partdec_us << " │ "
              << std::setw(10) << std::setprecision(1)
              << (all_partdec_us / full_all.avg_us * 100) << "% │\n";
    std::cout << "│ FinDec                   │ "
              << std::setw(10) << std::setprecision(2)
              << full_findec.avg_us << " │ "
              << std::setw(10) << std::setprecision(1)
              << (full_findec.avg_us / full_all.avg_us * 100) << "% │\n";
    std::cout << "├──────────────────────────┼────────────┼──────────────┤\n";
    std::cout << "│ 总计                     │ "
              << std::setw(10) << std::setprecision(2)
              << full_all.avg_us << " │    100.0% │\n";
    std::cout << "└──────────────────────────┴────────────┴──────────────┘\n";

    /* ─── 验证解密正确性 ─── */
    int mu_dec = decrypt_full(Ch, keys, params, 999);
    std::cout << "  解密验证: μ = " << mu_dec << " (期望 1)  "
              << (mu_dec == 1 ? "✓" : "✗") << "\n";
}

/* ══════════════════════════════════════════════════
   参数缩放测试
   ══════════════════════════════════════════════════ */
static void scaling_test()
{
    std::cout << "\n\n";
    std::cout << "╔═══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  参数缩放趋势: PartDec 瓶颈分析  (各操作平均耗时, μs)        ║\n";
    std::cout << "╠═════╤═════╤═══════╤═════╤══════╤═══════════╤════════╤════════╣\n";
    std::cout << "║  n  │  d  │   q   │  N  │   R  │ G⁻¹(ŵ^T) │ t·Ĉ·u │ Total  ║\n";
    std::cout << "╠═════╪═════╪═══════╪═════╪══════╪═══════════╪════════╪════════╣\n";

    struct Config { int n, d, N; long q; };
    std::vector<Config> configs = {
        {2, 1, 2, 65537},
        {2, 1, 3, 65537},
        {2, 1, 5, 65537},
        {3, 1, 3, 65537},
        {4, 1, 3, 65537},
        {2, 2, 3, 65537},
        {2, 3, 3, 65537},
        {4, 2, 3, 65537},
        {4, 2, 5, 65537},
    };

    const int REPS = 200;
    const int b = 2;

    for (auto& cfg : configs) {
        MIDParams p = MIDParams::make(cfg.n, cfg.d, cfg.q, cfg.N, b, 4, 1);

        std::mt19937_64 rng(42);
        std::uniform_int_distribution<long> unif(0, cfg.q - 1);

        Vec t_master(p.R);
        for (int i = 0; i < p.R - 1; ++i) t_master[i] = unif(rng);
        t_master[p.R - 1] = 1;

        auto keys = generate_shared_keys(t_master, cfg.N, cfg.q, rng);
        Mat A  = generate_lwe_matrix(t_master, p.R, p.M, cfg.q, 1, rng);
        Mat G  = build_gadget(p.R, cfg.q, b);
        Mat C  = gsw_encrypt(A, G, 1, cfg.q, rng);
        Mat Ch = simple_expand(C, cfg.N);

        /* G⁻¹ 耗时 */
        Vec wh = build_w_hat(p.R, cfg.q);
        Mat wc(p.R, Vec(1));
        for (int i = 0; i < p.R; ++i) wc[i][0] = wh[i];
        TimingResult t_ginv = bench([&](){ gadget_inverse(wc, cfg.q, b); }, REPS);

        /* γ_k 计算 (extract + vec_mat_mul + dot) × N 块 */
        Mat uc = gadget_inverse(wc, cfg.q, b);
        Vec u(p.M);
        for (int i = 0; i < p.M; ++i) u[i] = uc[i][0];

        TimingResult t_gamma = bench([&](){
            long gamma = 0;
            for (int j = 0; j < p.N_id; ++j) {
                Mat blk = extract_block(Ch, 0, j, p.R, p.M);
                Vec v   = vec_mat_mul(keys[0], blk, cfg.q);
                gamma = mod_pos(gamma + vec_dot_mod(v, u, cfg.q), cfg.q);
            }
        }, REPS);

        /* 完整 PartDec */
        TimingResult t_full = bench([&](){
            part_dec(Ch, 0, keys[0], p, 42);
        }, REPS);

        std::cout << "║ " << std::setw(3) << cfg.n
                  << " │ " << std::setw(3) << cfg.d
                  << " │ " << std::setw(5) << cfg.q
                  << " │ " << std::setw(3) << cfg.N
                  << " │ " << std::setw(4) << p.R
                  << " │ " << std::setw(9) << std::fixed << std::setprecision(1)
                  << t_ginv.avg_us
                  << " │ " << std::setw(6) << std::setprecision(1)
                  << t_gamma.avg_us
                  << " │ " << std::setw(6) << std::setprecision(1)
                  << t_full.avg_us
                  << " ║\n";
    }

    std::cout << "╚═════╧═════╧═══════╧═════╧══════╧═══════════╧════════╧════════╝\n";
    std::cout << "\n";
    std::cout << "  说明:\n";
    std::cout << "    R = (d+1)·n + 1,  M = R·⌈log₂q⌉\n";
    std::cout << "    G⁻¹(ŵ^T): gadget_inverse 的开销, 与 R·k 成正比\n";
    std::cout << "    t·Ĉ·u:    N 个块的 (提取+向量矩阵乘+内积), 是主要瓶颈\n";
    std::cout << "    Total:     完整 PartDec (含噪声采样等)\n";
    std::cout << "    瓶颈: 向量-矩阵乘 O(R·M) = O(R²·k), 随 n, d 二次增长\n";
}

/* ══════════════════════════════════════════════════
   main
   ══════════════════════════════════════════════════ */
int main()
{
    std::cout << "╔══════════════════════════════════════════════════════╗\n";
    std::cout << "║  多身份 FHE 解密 — 基本操作耗时拆解与缩放分析       ║\n";
    std::cout << "║  使用: matops.h, eval.h (gadget_inverse/build_gadget)║\n";
    std::cout << "╚══════════════════════════════════════════════════════╝\n";

    /* 详细拆解 (统一 128-bit 参数) */
    {
        auto mp = unified::default_midparams_128(1, 3);
        run_benchmark(mp.n, 1, mp.q, mp.N_id, 500);
    }

    /* 较大参数 (仍使用统一 modulus q; 调整深度与身份数) */
    {
        auto mp2 = unified::default_midparams_128(2, 5);
        run_benchmark(mp2.n, 2, mp2.q, mp2.N_id, 100);
    }

    /* 缩放趋势对比 */
    scaling_test();

    return 0;
}

