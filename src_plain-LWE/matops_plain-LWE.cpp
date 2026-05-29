/**
 * matops.cpp
 *
 * 完整复现图片里的 HIBE 第 ℓ 层委托步骤:
 *   H_{id_ℓ} = FRD(id_ℓ)
 *   T  = H_{id_ℓ} · G                       (mat_mul)
 *   B  = A_ℓ + T                            (mat_add)
 *   A_{id_ℓ} = [A_{id_{ℓ-1}} || B]          (mat_hcat)
 *   target = H_ℓ · G - B                    (mat_mul + mat_sub)
 *   ── 验证 ── A_{id_ℓ} · R_ℓ = target      (mat_mul)
 *
 * 并对每个矩阵操作做独立的微观基准测试,
 * 在 HIBE 典型维度下比较各操作的耗时。
 */

#include "matops_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include "unified_params_plain-LWE.h"
#include "bench_utils_plain-LWE.h"

using namespace matops;
using clk = std::chrono::high_resolution_clock;
using ms_t = std::chrono::duration<double, std::milli>;

static void hr() { std::cout << std::string(72, '-') << "\n"; }

static std::ostringstream bench_oss;  // 全局收集器

/* ════════════════════════════════════════════════════
   §1  正确性验证: 所有操作的代数性质
   ════════════════════════════════════════════════════ */
void test_correctness() {
    std::cout << "\n========== 正确性测试 ==========\n";
    auto u_p = unified::default_mp12_params(); long q = u_p.q;

    Mat A = random_mat(8, 6, q, 1);
    Mat B = random_mat(8, 6, q, 2);
    Mat C = random_mat(6, 5, q, 3);

    // 1. 加减法互逆: (A + B) - B = A
    Mat AB  = mat_add(A, B, q);
    Mat ABB = mat_sub(AB, B, q);
    std::cout << "  (A+B) - B == A           "
              << (mat_eq(ABB, A, q) ? "PASS" : "FAIL") << "\n";

    // 2. 乘法分配律: A·(B+B') = A·B + A·B'  (用 D 代替 B' 因为维度)
    Mat D = random_mat(8, 6, q, 4);
    // (注意: A·(B+D) 这里 A 不能左乘 (B+D), 我们换成 (A+B)·C
    Mat AplusB_C = mat_mul(mat_add(A, B, q), C, q);
    Mat AC_plus_BC = mat_add(mat_mul(A, C, q), mat_mul(B, C, q), q);
    std::cout << "  (A+B)·C == A·C + B·C     "
              << (mat_eq(AplusB_C, AC_plus_BC, q) ? "PASS" : "FAIL") << "\n";

    // 3. 横向拼接: [A|B] · [C; C'] = A·C + B·C' (按块)
    //    用 [A|B] · [C; D'] 验证
    Mat AB_h = mat_hcat(A, B);          // 8 × 12
    Mat C2   = random_mat(6, 5, q, 5);  // C 部分
    Mat D2   = random_mat(6, 5, q, 6);  // D' 部分
    Mat CD_v = mat_vcat(C2, D2);        // 12 × 5
    Mat lhs  = mat_mul(AB_h, CD_v, q);  // 8 × 5
    Mat rhs  = mat_add(mat_mul(A, C2, q), mat_mul(B, D2, q), q);
    std::cout << "  [A|B]·[C;D] == A·C + B·D "
              << (mat_eq(lhs, rhs, q) ? "PASS" : "FAIL") << "\n";

    // 4. 维度检查
    auto [r, c] = dim(AB_h);
    std::cout << "  hcat 维度: 8×6 + 8×6 → " << r << "×" << c << "\n";
    auto [r2, c2] = dim(CD_v);
    std::cout << "  vcat 维度: 6×5 + 6×5 → " << r2 << "×" << c2 << "\n";
}

/* ════════════════════════════════════════════════════
   §2  完整复现图片里的 HIBE 委托步骤
   ════════════════════════════════════════════════════ */
void simulate_hibe_delegation_step() {
    std::cout << "\n========== HIBE 第 ℓ 层委托步骤 ==========\n";

    // HIBE 参数
    auto u_p = unified::default_mp12_params(); long q = u_p.q;
    int  n  = 8;
    int  k  = 14;          // ⌈log_2 q⌉
    int  nk = n * k;       // = 112
    int  level = 3;        // 当前层级 ℓ
    int  m_prev = n + (level - 1) * nk;   // 上层 A 的列数
    int  m_curr = m_prev + nk;            // 本层 A 的列数

    std::cout << "  参数: n=" << n << ", q=" << q << ", k=" << k
              << ", nk=" << nk << ", ℓ=" << level << "\n";
    std::cout << "  m_{ℓ-1} = " << m_prev << ", m_ℓ = " << m_curr << "\n\n";

    // ─── 准备输入 ───
    Mat A_prev = random_mat(n, m_prev, q, 100);   // A_{id_{ℓ-1}}
    Mat A_l    = random_mat(n, nk,     q, 101);   // A_ℓ
    Mat G      = random_mat(n, nk,     q, 102);   // gadget (用随机矩阵代替)
    Mat H_id   = random_mat(n, n,      q, 103);   // FRD(id_ℓ)
    Mat H_l    = random_mat(n, n,      q, 104);   // 另一个 tag

    std::cout << "  A_{ℓ-1}: " << n << "×" << m_prev << "\n";
    std::cout << "  A_ℓ    : " << n << "×" << nk << "\n";
    std::cout << "  G      : " << n << "×" << nk << "\n";
    std::cout << "  H_{id_ℓ}: " << n << "×" << n << "\n";

    // ─── Step 1: T = H_{id_ℓ} · G ───
    auto t0 = clk::now();
    Mat T = mat_mul(H_id, G, q);
    auto t1 = clk::now();
    double ms_T = ms_t(t1 - t0).count();

    // ─── Step 2: B = A_ℓ + T ───
    auto t2 = clk::now();
    Mat B = mat_add(A_l, T, q);
    auto t3 = clk::now();
    double ms_B = ms_t(t3 - t2).count();

    // ─── Step 3: A_{id_ℓ} = [A_{ℓ-1} | B] ───
    auto t4 = clk::now();
    Mat A_curr = mat_hcat(A_prev, B);
    auto t5 = clk::now();
    double ms_cat = ms_t(t5 - t4).count();

    // ─── Step 4: target = H_ℓ · G - B ───
    auto t6 = clk::now();
    Mat HlG = mat_mul(H_l, G, q);
    Mat target = mat_sub(HlG, B, q);
    auto t7 = clk::now();
    double ms_target = ms_t(t7 - t6).count();

    auto [rT, cT] = dim(T);
    auto [rA, cA] = dim(A_curr);
    auto [rTg, cTg] = dim(target);

    std::ostringstream local_oss;
    local_oss << "\n--- HIBE Delegation Step Breakdown ---\n"
              << "  Step 1  T = H_{id}·G       : " << rT << "×" << cT
              << "    " << std::fixed << std::setprecision(3) << ms_T << " ms\n"
              << "  Step 2  B = A_ℓ + T        : " << rT << "×" << cT
              << "    " << ms_B << " ms\n"
              << "  Step 3  [A_{ℓ-1} || B]    : " << rA << "×" << cA
              << "  " << ms_cat << " ms\n"
              << "  Step 4  H_ℓ·G - B         : " << rTg << "×" << cTg
              << "    " << ms_target << " ms\n";

    std::cout << "\n  Step 1  T = H_{id}·G       : " << rT << "×" << cT
              << "    " << std::fixed << std::setprecision(3) << ms_T << " ms\n";
    std::cout << "  Step 2  B = A_ℓ + T        : " << rT << "×" << cT
              << "    " << ms_B << " ms\n";
    std::cout << "  Step 3  [A_{ℓ-1} || B]    : " << rA << "×" << cA
              << "  " << ms_cat << " ms\n";
    std::cout << "  Step 4  H_ℓ·G - B         : " << rTg << "×" << cTg
              << "    " << ms_target << " ms\n";

    double total = ms_T + ms_B + ms_cat + ms_target;
    std::cout << "  ─────────────────────────────────────────\n";
    std::cout << "  总耗时(单次委托步骤)   : "
              << std::setprecision(3) << total << " ms\n";
    local_oss << "  ─────────────────────────────────────────\n";
    local_oss << "  总耗时(单次委托步骤)   : "
              << std::setprecision(3) << total << " ms\n";

    bench_oss << local_oss.str() << "\n";

    // 防止编译器优化掉
    long sink = A_curr[0][0] + target[0][0];
    (void)sink;  // prevent compiler optimization
}

/* ════════════════════════════════════════════════════
   §3  四个操作的微观基准测试
   ════════════════════════════════════════════════════ */

void run_benchmarks() {
    std::cout << "\n========== 微观基准测试 ==========\n\n";
    auto u_p = unified::default_mp12_params(); long q = u_p.q;

    constexpr int TOTAL_ITERS = 1000;
    constexpr int WARMUP = 50;

    struct Dim { int r, c; const char* name; };
    Dim sizes[] = {
        { 16,  16, "tiny    " },
        { 64,  64, "small   " },
        {128, 128, "medium  " },
        {256, 256, "large   " },
        { 32, 256, "wide    " },
        {256,  32, "tall    " },
    };

    /* ───── 加法 ───── */
    std::cout << ">>> mat_add: 加法 (A+B mod q)\n";
    std::cout << "  shape          time(ms)    throughput      median(us)  stddev(us)\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.r, s.c, q, 2);
        auto stats = run_benchmark([&]{ Mat C = mat_add(A, B, q); (void)C; }, TOTAL_ITERS, WARMUP);
        double t = stats.avg_us / 1000.0;
        double mops = (double)s.r * s.c / (t * 1000.0);
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "    " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(1) << mops << " M ops/s"
                  << "  " << std::setw(8) << std::setprecision(1) << stats.median_us
                  << "  " << std::setw(8) << std::setprecision(1) << stats.stddev_us << "\n";
        bench_oss << "  " << s.name << std::setw(4) << s.r << "x" << std::setw(4) << s.c
                  << "  " << std::setw(10) << std::setprecision(4) << t
                  << " ms  " << std::setw(8) << std::setprecision(1) << mops << " M ops/s"
                  << "  median=" << std::setprecision(1) << stats.median_us << "us"
                  << "  (total=" << TOTAL_ITERS << " warmup=" << WARMUP
                  << " active=" << stats.active_samples << ")\n";
    }

    /* ───── 减法 ───── */
    std::cout << "\n>>> mat_sub: 减法 (A-B mod q)\n";
    std::cout << "  shape          time(ms)    throughput      median(us)  stddev(us)\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.r, s.c, q, 2);
        auto stats = run_benchmark([&]{ Mat C = mat_sub(A, B, q); (void)C; }, TOTAL_ITERS, WARMUP);
        double t = stats.avg_us / 1000.0;
        double mops = (double)s.r * s.c / (t * 1000.0);
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "    " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(1) << mops << " M ops/s"
                  << "  " << std::setw(8) << std::setprecision(1) << stats.median_us
                  << "  " << std::setw(8) << std::setprecision(1) << stats.stddev_us << "\n";
        bench_oss << "  " << s.name << std::setw(4) << s.r << "x" << std::setw(4) << s.c
                  << "  " << std::setw(10) << std::setprecision(4) << t
                  << " ms  " << std::setw(8) << std::setprecision(1) << mops << " M ops/s"
                  << "  median=" << std::setprecision(1) << stats.median_us << "us"
                  << "  (total=" << TOTAL_ITERS << " warmup=" << WARMUP
                  << " active=" << stats.active_samples << ")\n";
    }

    /* ───── 横向拼接 ───── */
    std::cout << "\n>>> mat_hcat: 横向拼接 [A|B]\n";
    std::cout << "  shape          time(ms)    throughput      median(us)  stddev(us)\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.r, s.c, q, 2);
        auto stats = run_benchmark([&]{ Mat C = mat_hcat(A, B); (void)C; }, TOTAL_ITERS, WARMUP);
        double t = stats.avg_us / 1000.0;
        double mops = (double)s.r * (s.c * 2) / (t * 1000.0);
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "    " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(1) << mops << " M cells/s"
                  << "  " << std::setw(8) << std::setprecision(1) << stats.median_us
                  << "  " << std::setw(8) << std::setprecision(1) << stats.stddev_us << "\n";
        bench_oss << "  " << s.name << std::setw(4) << s.r << "x" << std::setw(4) << s.c
                  << "  " << std::setw(10) << std::setprecision(4) << t
                  << " ms  " << std::setw(8) << std::setprecision(1) << mops << " M cells/s"
                  << "  median=" << std::setprecision(1) << stats.median_us << "us"
                  << "  (total=" << TOTAL_ITERS << " warmup=" << WARMUP
                  << " active=" << stats.active_samples << ")\n";
    }

    /* ───── 乘法 ───── */
    std::cout << "\n>>> mat_mul: 乘法 (A·B mod q)\n";
    std::cout << "  A shape    B shape       time(ms)    throughput      median(us)  stddev(us)\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.c, s.r, q, 2);
        auto stats = run_benchmark([&]{ Mat C = mat_mul(A, B, q); (void)C; }, TOTAL_ITERS, WARMUP);
        double t = stats.avg_us / 1000.0;
        long flops_per_call = 2L * s.r * s.c * s.r;
        double gflops = (double)flops_per_call / (t * 1e6);
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "  " << std::setw(4) << s.c << "×" << std::setw(4) << s.r
                  << "   " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(2) << gflops << " GOPS"
                  << "  " << std::setw(8) << std::setprecision(1) << stats.median_us
                  << "  " << std::setw(8) << std::setprecision(1) << stats.stddev_us << "\n";
    }

    /* ───── 操作间相对耗时 ───── */
    std::cout << "\n>>> 同一尺寸下不同操作的相对耗时 (128×128)\n";
    int N = 128;
    Mat A = random_mat(N, N, q, 1);
    Mat B = random_mat(N, N, q, 2);
    auto s_add = run_benchmark([&]{ Mat C = mat_add(A, B, q); (void)C; }, TOTAL_ITERS, WARMUP);
    auto s_sub = run_benchmark([&]{ Mat C = mat_sub(A, B, q); (void)C; }, TOTAL_ITERS, WARMUP);
    auto s_cat = run_benchmark([&]{ Mat C = mat_hcat(A, B); (void)C; }, TOTAL_ITERS, WARMUP);
    auto s_mul = run_benchmark([&]{ Mat C = mat_mul(A, B, q); (void)C; }, TOTAL_ITERS, WARMUP);
    double t_add = s_add.avg_us / 1000.0;
    double t_sub = s_sub.avg_us / 1000.0;
    double t_cat = s_cat.avg_us / 1000.0;
    double t_mul = s_mul.avg_us / 1000.0;
    std::cout << "  add :  " << std::setw(8) << std::fixed << std::setprecision(4) << t_add << " ms  (基准  1×)\n";
    std::cout << "  sub :  " << std::setw(8) << t_sub << " ms  ("
              << std::setprecision(2) << t_sub/t_add << "×)\n";
    std::cout << "  hcat:  " << std::setw(8) << std::setprecision(4) << t_cat << " ms  ("
              << std::setprecision(2) << t_cat/t_add << "×)\n";
    std::cout << "  mul :  " << std::setw(8) << std::setprecision(4) << t_mul << " ms  ("
              << std::setprecision(2) << t_mul/t_add << "×)\n";

    bench_oss << "\n=== Benchmark: Matrix Ops Relative Timing (128×128) ===\n"
              << "  Parameters: q=" << q << "\n"
              << "  total_iters=" << TOTAL_ITERS << " warmup=" << WARMUP << "\n"
              << "  add :  " << std::setw(8) << std::fixed << std::setprecision(4)
              << t_add << " ms  (基准  1×)\n"
              << "  sub :  " << std::setw(8) << t_sub << " ms  ("
              << std::setprecision(2) << t_sub/t_add << "×)\n"
              << "  hcat:  " << std::setw(8) << std::setprecision(4) << t_cat << " ms  ("
              << std::setprecision(2) << t_cat/t_add << "×)\n"
              << "  mul :  " << std::setw(8) << std::setprecision(4) << t_mul << " ms  ("
              << std::setprecision(2) << t_mul/t_add << "×)\n";
}

void run_bench_matops() {
    // 重置收集器
    bench_oss.str("");
    bench_oss.clear();

    bench_oss << "=== Benchmark: Matrix Operations (matops.h) ===\n"
              << "  API: mat_add / mat_sub / mat_hcat / mat_mul\n\n";

    std::cout << "════════════════════════════════════════════════════════════════\n";
    std::cout << "  HIBE Matrix Operations Benchmark\n";
    std::cout << "  (mat_add / mat_sub / mat_mul / mat_hcat)\n";
    std::cout << "════════════════════════════════════════════════════════════════\n";

    test_correctness();
    simulate_hibe_delegation_step();
    run_benchmarks();

    // 统一写入文件
    bench_write(bench_oss.str());

    std::cout << "\nDone.\n";
}
