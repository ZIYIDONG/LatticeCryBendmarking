/**
 * bench_matops.cpp
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

#include "../include/matops.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <vector>
#include <string>
#include "unified_params.h"

using namespace matops;
using clk = std::chrono::high_resolution_clock;
using ms_t = std::chrono::duration<double, std::milli>;

/* ════════════════════════════════════════════════════
   计时器 — 重复 N 次取均值,首次预热
   ════════════════════════════════════════════════════ */
template <typename F>
double time_ms(int repeats, F&& f) {
    f();  // warmup
    auto t0 = clk::now();
    for (int i = 0; i < repeats; i++) f();
    auto t1 = clk::now();
    return ms_t(t1 - t0).count() / repeats;
}

static void hr() { std::cout << std::string(72, '-') << "\n"; }

/* ════════════════════════════════════════════════════
   §1  正确性验证: 所有操作的代数性质
   ════════════════════════════════════════════════════ */
void test_correctness() {
    std::cout << "\n========== 正确性测试 ==========\n";
    auto __u_p = unified::default_mp12_params_128(); long q = __u_p.q;

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
    auto __u_p = unified::default_mp12_params_128(); long q = __u_p.q;
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

    // 防止编译器优化掉
    volatile long sink = A_curr[0][0] + target[0][0];
    (void)sink;
}

/* ════════════════════════════════════════════════════
   §3  四个操作的微观基准测试
   ════════════════════════════════════════════════════ */
struct BenchResult {
    std::string op;
    int dim_a_r, dim_a_c, dim_b_r, dim_b_c;
    int dim_out_r, dim_out_c;
    double time_ms;
    double throughput_mops;   // 百万次基本运算/秒
};

void run_benchmarks() {
    std::cout << "\n========== 微观基准测试 ==========\n\n";
    long q = 8209;

    // 测试维度: 模拟 MP12 / HIBE 中常见的矩阵尺寸
    struct Dim { int r, c; const char* name; };
    Dim sizes[] = {
        { 16,  16, "tiny    " },
        { 64,  64, "small   " },
        {128, 128, "medium  " },
        {256, 256, "large   " },
        { 32, 256, "wide    " },     // n × nk 形状
        {256,  32, "tall    " },
    };

    std::vector<BenchResult> results;

    /* ───── 加法 ───── */
    std::cout << ">>> mat_add: 加法 (A+B mod q)\n";
    std::cout << "  shape          time(ms)    throughput\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.r, s.c, q, 2);
        int reps = std::max(10, 5000000 / (s.r * s.c + 1));
        double t = time_ms(reps, [&]{ Mat C = mat_add(A, B, q); (void)C; });
        double mops = (double)s.r * s.c / (t * 1000.0);  // 百万 ops/sec
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "    " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(1) << mops << " M ops/s\n";
        results.push_back({"add", s.r, s.c, s.r, s.c, s.r, s.c, t, mops});
    }

    /* ───── 减法 ───── */
    std::cout << "\n>>> mat_sub: 减法 (A-B mod q)\n";
    std::cout << "  shape          time(ms)    throughput\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.r, s.c, q, 2);
        int reps = std::max(10, 5000000 / (s.r * s.c + 1));
        double t = time_ms(reps, [&]{ Mat C = mat_sub(A, B, q); (void)C; });
        double mops = (double)s.r * s.c / (t * 1000.0);
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "    " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(1) << mops << " M ops/s\n";
        results.push_back({"sub", s.r, s.c, s.r, s.c, s.r, s.c, t, mops});
    }

    /* ───── 横向拼接 ───── */
    std::cout << "\n>>> mat_hcat: 横向拼接 [A|B]\n";
    std::cout << "  shape          time(ms)    throughput\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.r, s.c, q, 2);
        int reps = std::max(10, 5000000 / (s.r * s.c + 1));
        double t = time_ms(reps, [&]{ Mat C = mat_hcat(A, B); (void)C; });
        double mops = (double)s.r * (s.c * 2) / (t * 1000.0);
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "    " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(1) << mops << " M cells/s\n";
        results.push_back({"hcat", s.r, s.c, s.r, s.c, s.r, s.c*2, t, mops});
    }

    /* ───── 乘法 ───── */
    std::cout << "\n>>> mat_mul: 乘法 (A·B mod q)\n";
    std::cout << "  A shape    B shape       time(ms)    throughput\n";
    for (auto& s : sizes) {
        Mat A = random_mat(s.r, s.c, q, 1);
        Mat B = random_mat(s.c, s.r, q, 2);  // 互换让乘积合法
        // 总浮点运算量: 2·r·c·r (mul-add)
        long flops_per_call = 2L * s.r * s.c * s.r;
        int reps = std::max(3, (int)(2000000000L / (flops_per_call + 1)));
        if (reps > 200) reps = 200;
        double t = time_ms(reps, [&]{ Mat C = mat_mul(A, B, q); (void)C; });
        double gflops = (double)flops_per_call / (t * 1e6);  // GFLOPS
        std::cout << "  " << std::setw(4) << s.r << "×" << std::setw(4) << s.c
                  << "  " << std::setw(4) << s.c << "×" << std::setw(4) << s.r
                  << "   " << std::setw(10) << std::fixed << std::setprecision(4) << t
                  << "    " << std::setw(8) << std::setprecision(2) << gflops << " GOPS\n";
        results.push_back({"mul", s.r, s.c, s.c, s.r, s.r, s.r, t, gflops*1000});
    }

    /* ───── 操作间相对耗时 ───── */
    std::cout << "\n>>> 同一尺寸下不同操作的相对耗时 (128×128)\n";
    int N = 128;
    Mat A = random_mat(N, N, q, 1);
    Mat B = random_mat(N, N, q, 2);
    double t_add = time_ms(50, [&]{ Mat C = mat_add(A, B, q); (void)C; });
    double t_sub = time_ms(50, [&]{ Mat C = mat_sub(A, B, q); (void)C; });
    double t_cat = time_ms(50, [&]{ Mat C = mat_hcat(A, B); (void)C; });
    double t_mul = time_ms(20, [&]{ Mat C = mat_mul(A, B, q); (void)C; });
    std::cout << "  add :  " << std::setw(8) << std::fixed << std::setprecision(4) << t_add << " ms  (基准  1×)\n";
    std::cout << "  sub :  " << std::setw(8) << t_sub << " ms  ("
              << std::setprecision(2) << t_sub/t_add << "×)\n";
    std::cout << "  hcat:  " << std::setw(8) << std::setprecision(4) << t_cat << " ms  ("
              << std::setprecision(2) << t_cat/t_add << "×)\n";
    std::cout << "  mul :  " << std::setw(8) << std::setprecision(4) << t_mul << " ms  ("
              << std::setprecision(2) << t_mul/t_add << "×)\n";
}

int main() {
    std::cout << "════════════════════════════════════════════════════════════════\n";
    std::cout << "  HIBE Matrix Operations Benchmark\n";
    std::cout << "  (mat_add / mat_sub / mat_mul / mat_hcat)\n";
    std::cout << "════════════════════════════════════════════════════════════════\n";

    test_correctness();
    simulate_hibe_delegation_step();
    run_benchmarks();

    std::cout << "\nDone.\n";
    return 0;
}
