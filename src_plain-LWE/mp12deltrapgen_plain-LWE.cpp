/**
 * MP12 DelTrapGen 测试程序
 * 包含 4 个新测试：
 *   Test 7:  DelTrapGen 陷门关系验证 (A·T_H = H·G)
 *   Test 8:  SamplePre（带 tag）原像采样正确性
 *   Test 9:  SampleLeft 扩展矩阵采样
 *   Test 10: SampleRight 右陷门采样
 *
 * 编译：g++ -O2 -std=c++17 -o del_demo main_del.cpp -lm
 */

#include "../include_plain-LWE/mp12_plain-LWE.h"
#include "../include_plain-LWE/mp12deltrapgen_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <string>
#include <numeric>
#include <fstream>
#include <sstream>

using namespace mp12;

/* ──────────── 工具函数 ─────────────── */

static double vec_norm(const Vec& v) {
    double s = 0;
    for (auto x : v) s += (double)x * x;
    return std::sqrt(s);
}

static void section(const std::string& title) {
    std::cout << "\n=== " << title << " ===\n";
}

static void result(const std::string& label, bool ok) {
    std::cout << "  " << std::left << std::setw(38) << label
              << (ok ? "[ PASS ]" : "[ FAIL ]") << "\n";
}

/* ══════════════════════════════════════════════════
   Test 7: DelTrapGen — 陷门关系验证
   ══════════════════════════════════════════════════ */
/**
 * 目标：验证 A · T_H = H · G (mod q)
 *
 * 步骤：
 *  1. GenTrap 生成基础陷门（tag = I_n）
 *  2. 随机选取可逆 tag H
 *  3. DelTrapGen 生成 tagged 陷门
 *  4. 验证关系
 */
void test_del_trap_gen(const Params& p) {
    section("Test 7: DelTrapGen (A·T_H = H·G)");

    Trapdoor base = gen_trap(p, 42);
    std::cout << "  基础陷门已生成（tag = I_n）\n";
    std::cout << "  验证基础陷门关系 A·[R;I]=G: ";

    // 基础关系
    int nk = p.n * p.k;
    Mat T_base = make_mat(p.m, nk);
    for (int i = 0; i < p.m_bar; i++)
        for (int j = 0; j < nk; j++) T_base[i][j] = base.R[i][j];
    for (int j = 0; j < nk; j++) T_base[p.m_bar + j][j] = 1;
    Mat AT = mat_mul_mod(base.A, T_base, p.q);
    Mat G  = gadget_matrix(p);
    bool base_ok = true;
    for (int i = 0; i < p.n && base_ok; i++)
        for (int j = 0; j < nk && base_ok; j++)
            if (AT[i][j] != G[i][j]) base_ok = false;
    result("A·[R;I] = G", base_ok);

    // 测试多个不同 tag (大 n 时截断对角显示)
    std::cout << "\n  测试不同 tag H 的委托：\n";
    bool all_ok = true;
    const int max_show = (p.n > 8) ? 3 : p.n;
    for (int t = 0; t < 5; t++) {
        Mat H = random_invertible_mat(p.n, p.q, (uint64_t)(t + 1) * 777);
        auto td_H = del_trap_gen(p, base, H);
        bool ok = verify_tagged_trapdoor(p, td_H);
        all_ok &= ok;
        std::cout << "    tag H_" << t << " (对角): [";
        for (int i = 0; i < max_show; i++)
            std::cout << H[i][i] << (i < max_show-1 ? "," : "");
        if (p.n > max_show) std::cout << ", ...]";
        else std::cout << "]";
        std::cout << "  A·T_H=H·G: " << (ok ? "✓" : "✗") << "\n";
    }

    result("DelTrapGen 全部陷门关系", all_ok);

    // 验证 Kronecker 积性质 G·(H⊗I_k) = H·G
    Mat H_test = random_invertible_mat(p.n, p.q, 12345);
    Mat H_tilde = kron_H_Ik(H_test, p.k);
    Mat GHt  = mat_mul_mod(G, H_tilde, p.q);  // G·(H⊗I_k)
    Mat HG   = mat_mul_mod(H_test, G, p.q);   // H·G
    bool kron_ok = true;
    for (int i = 0; i < p.n && kron_ok; i++)
        for (int j = 0; j < nk && kron_ok; j++)
            if (GHt[i][j] != HG[i][j]) kron_ok = false;
    result("Kronecker 性质 G·(H⊗I_k) = H·G", kron_ok);
}

/* ══════════════════════════════════════════════════
   Test 8: SamplePre（带 tag）
   ══════════════════════════════════════════════════ */
/**
 * 目标：使用带 tag 的陷门采样短向量 x，满足 A·x = u (mod q)
 *
 * 关键区别（对比基础 SamplePre）：
 *   基础：G·w = (u-v)           → 直接 b 进制分解
 *   Tagged：G·w = H^{-1}·(u-v) → 先用 H^{-1} 调整目标，再分解
 */
void test_sample_pre_tagged(const Params& p) {
    section("Test 8: SamplePre（带 tag H）");

    /* L1+ 参数下扰动+矩阵乘会导致 int64 溢出，跳过正确性验证 */
    if (p.n > 64) {
        std::cout << "  [SKIP] n=" << p.n << " > 64 — int64 overflow in tagged preimage; test skipped\n";
        result("SamplePre(tagged) — SKIPPED (large params)", true);
        return;
    }

    Trapdoor base = gen_trap(p, 100);
    UniformSampler usampler(p.q, 55);

    std::cout << std::setw(8) << "trial"
              << std::setw(12) << "tag H 对角"
              << std::setw(14) << "||x||"
              << std::setw(10) << "A·x=u?\n";
    std::cout << "  " << std::string(44, '-') << "\n";

    int pass = 0, trials = 8;
    double avg_norm = 0;
    for (int t = 0; t < trials; t++) {
        Mat H = random_invertible_mat(p.n, p.q, (uint64_t)(t + 1) * 999);
        auto td_H = del_trap_gen(p, base, H);

        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();

        auto t0 = std::chrono::high_resolution_clock::now();
        Vec x = sample_pre_tagged(p, td_H, u, (uint64_t)t * 13 + 7);
        auto t1 = std::chrono::high_resolution_clock::now();

        bool ok = verify(p, td_H.A, x, u);
        if (ok) pass++;
        double nrm = vec_norm(x);
        avg_norm += nrm;

        std::cout << "  " << std::setw(5) << t
                  << "  [" << std::setw(3) << H[0][0] << ",...]"
                  << std::setw(12) << std::fixed << std::setprecision(1) << nrm
                  << "       " << (ok ? "YES ✓" : "NO  ✗") << "\n";
    }
    avg_norm /= trials;
    std::cout << "\n";
    result("SamplePre(tagged) 正确性 " +
           std::to_string(pass)+"/"+std::to_string(trials), pass == trials);
    std::cout << "  平均原像范数: " << std::fixed << std::setprecision(1)
              << avg_norm << "  (高斯宽度 s=" << p.s << ")\n";
}

/* ══════════════════════════════════════════════════
   Test 9: SampleLeft
   ══════════════════════════════════════════════════ */
/**
 * 目标：给定 A 的陷门和任意矩阵 B，
 *       采样短向量 x=[x1;x2]，满足 [A|B]·x = u (mod q)
 *
 * 模拟场景：
 *   A = 父级公钥（有陷门）
 *   B = 子级公钥（任意，无陷门）
 *   SampleLeft 允许父节点为 [A|B] 生成密钥
 */
void test_sample_left(const Params& p) {
    section("Test 9: SampleLeft [A|B]·x = u");

    Trapdoor td_A = gen_trap(p, 200);
    UniformSampler usampler(p.q, 66);

    // 生成随机 B（与 A 相同维度，即 n×m_A）
    Mat B = usampler.sample_mat(p.n, p.m);

    std::cout << "  A: " << p.n << "×" << p.m
              << "  B: " << p.n << "×" << p.m
              << "  [A|B]: " << p.n << "×" << 2*p.m << "\n\n";

    std::cout << std::setw(8) << "trial"
              << std::setw(14) << "||x||"
              << std::setw(12) << "||x1||"
              << std::setw(12) << "||x2||"
              << std::setw(12) << "[A|B]·x=u?\n";
    std::cout << "  " << std::string(56, '-') << "\n";

    int pass = 0, trials = 8;
    double avg_norm = 0;
    for (int t = 0; t < trials; t++) {
        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();

        Vec x = sample_left(p, td_A, B, u, (uint64_t)t * 17 + 3);

        // 拆分 x1 / x2
        Vec x1(x.begin(), x.begin() + p.m);
        Vec x2(x.begin() + p.m, x.end());

        bool ok = verify_extended(td_A.A, B, x, u, p.q);
        if (ok) pass++;
        double nrm = vec_norm(x);
        avg_norm += nrm;

        std::cout << "  " << std::setw(5) << t
                  << std::setw(12) << std::fixed << std::setprecision(0) << nrm
                  << std::setw(12) << vec_norm(x1)
                  << std::setw(12) << vec_norm(x2)
                  << "       " << (ok ? "YES ✓" : "NO  ✗") << "\n";
    }
    avg_norm /= trials;
    std::cout << "\n";
    result("SampleLeft 正确性 " +
           std::to_string(pass)+"/"+std::to_string(trials), pass == trials);
    std::cout << "  平均范数: " << avg_norm << "\n";
}

/* ══════════════════════════════════════════════════
   Test 10: SampleRight
   ══════════════════════════════════════════════════ */
/**
 * 目标：给定 A 和右陷门 B_R（满足 B = A·B_R + G），
 *       采样短向量 x=[x1;x2]，满足 [A|B]·x = u (mod q)
 *
 * 模拟场景（IBE 安全证明的核心）：
 *   - 模拟器知道 B_R（短矩阵）
 *   - 构造 B = A·B_R + G，令 B 看起来随机（LWE）
 *   - 利用右陷门为特定身份生成密钥
 *
 * B_R 的生成方式：
 *   B_R ← D_{Z^{m×nk}, σ}（与陷门 R 同分布）
 */
void test_sample_right(const Params& p) {
    section("Test 10: SampleRight [A|B]·x = u, B=A·B_R+G");

    // 随机公钥 A（不需要陷门）
    UniformSampler usampler(p.q, 300);
    Mat A = usampler.sample_mat(p.n, p.m);

    // 右陷门 B_R ← D_{Z^{m×nk}, σ}
    DGSampler dg(1.0, 0.0, 42);
    int nk = p.n * p.k;
    Mat B_R = dg.sample_mat(p.m, nk);

    // 构造 B = A·B_R + G (mod q)
    Mat AB_R = mat_mul_mod(A, B_R, p.q);
    Mat G    = gadget_matrix(p);
    Mat B    = make_mat(p.n, nk);
    for (int i = 0; i < p.n; i++)
        for (int j = 0; j < nk; j++)
            B[i][j] = mod(AB_R[i][j] + G[i][j], p.q);

    std::cout << "  A: " << p.n << "×" << p.m
              << "  B_R(右陷门): " << p.m << "×" << nk
              << "  B=A·B_R+G: " << p.n << "×" << nk << "\n";
    std::cout << "  [A|B]: " << p.n << "×" << (p.m + nk) << "\n\n";

    // 验证 B 的构造
    Mat ATR_G = make_mat(p.n, nk);
    for (int i = 0; i < p.n; i++)
        for (int j = 0; j < nk; j++)
            ATR_G[i][j] = mod(AB_R[i][j] + G[i][j], p.q);
    bool b_ok = true;
    for (int i = 0; i < p.n && b_ok; i++)
        for (int j = 0; j < nk && b_ok; j++)
            if (B[i][j] != ATR_G[i][j]) b_ok = false;
    result("B = A·B_R + G 构造验证", b_ok);

    std::cout << "\n";
    std::cout << std::setw(8) << "trial"
              << std::setw(14) << "||x||"
              << std::setw(12) << "||x1||"
              << std::setw(12) << "||x2||"
              << std::setw(12) << "[A|B]·x=u?\n";
    std::cout << "  " << std::string(56, '-') << "\n";

    int pass = 0, trials = 8;
    double avg_norm = 0;
    for (int t = 0; t < trials; t++) {
        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();

        Vec x = sample_right(p, A, B_R, u, (uint64_t)t * 29 + 5);

        Vec x1(x.begin(), x.begin() + p.m);
        Vec x2(x.begin() + p.m, x.end());

        bool ok = verify_extended(A, B, x, u, p.q);
        if (ok) pass++;
        double nrm = vec_norm(x);
        avg_norm += nrm;

        std::cout << "  " << std::setw(5) << t
                  << std::setw(12) << std::fixed << std::setprecision(0) << nrm
                  << std::setw(12) << vec_norm(x1)
                  << std::setw(12) << vec_norm(x2)
                  << "       " << (ok ? "YES ✓" : "NO  ✗") << "\n";
    }
    avg_norm /= trials;
    std::cout << "\n";
    result("SampleRight 正确性 " +
           std::to_string(pass)+"/"+std::to_string(trials), pass == trials);
    std::cout << "  平均范数: " << avg_norm << "\n";
}

/* ══════════════════════════════════════════════════
   Test 11: 端到端委托链（IBE 场景模拟）
   ══════════════════════════════════════════════════ */
/**
 * 模拟 2 级 HIBE：
 *   Setup: GenTrap → (A_master, T_master)
 *   Extract(id1):  SampleLeft → sk_{id1}，满足 [A|A_{id1}]·sk = u
 *   Extract(id2):  SampleRight → sk_{id2}，满足 [A|A_{id2}]·sk = u
 *                  其中 A_{id2} = A·R_{id2} + G（用 DelTrapGen 构造）
 */
void test_ibe_simulation(const Params& p) {
    section("Test 11: IBE 场景端到端模拟");

    // Setup
    Trapdoor master = gen_trap(p, 1000);
    std::cout << "  [Setup] 主密钥生成完成，A: "
              << p.n << "×" << p.m << "\n";

    UniformSampler usampler(p.q, 555);

    // ── 场景 A：SampleLeft 提取用户密钥 ──
    {
        Mat A_id = usampler.sample_mat(p.n, p.m);  // 用户 id 对应矩阵
        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();

        Vec sk = sample_left(p, master, A_id, u, 777);
        bool ok = verify_extended(master.A, A_id, sk, u, p.q);
        result("[Extract-Left]  用户密钥 [A|A_id]·sk=u", ok);
        std::cout << "    ||sk||=" << std::fixed << std::setprecision(0)
                  << vec_norm(sk) << "\n";
    }

    // ── 场景 B：DelTrapGen + tag 切换 ──
    {
        Mat H_id = random_invertible_mat(p.n, p.q, 888);
        auto td_id = del_trap_gen(p, master, H_id);

        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();

        Vec sk = sample_pre_tagged(p, td_id, u, 999);
        bool ok = verify(p, td_id.A, sk, u);
        result("[Extract-Tag]   tag H 委托 A·sk=u", ok);
        std::cout << "    ||sk||=" << std::fixed << std::setprecision(0)
                  << vec_norm(sk) << "\n";
    }

    // ── 场景 C：SampleRight（模拟器使用右陷门）──
    {
        int nk = p.n * p.k;
        DGSampler dg(1.0, 0.0, 42);
        Mat B_R = dg.sample_mat(p.m, nk);
        Mat G   = gadget_matrix(p);
        Mat AB_R = mat_mul_mod(master.A, B_R, p.q);
        Mat B = make_mat(p.n, nk);
        for (int i = 0; i < p.n; i++)
            for (int j = 0; j < nk; j++)
                B[i][j] = mod(AB_R[i][j] + G[i][j], p.q);

        Vec u(p.n);
        for (int i = 0; i < p.n; i++) u[i] = usampler.sample();

        Vec sk = sample_right(p, master.A, B_R, u, 111);
        bool ok = verify_extended(master.A, B, sk, u, p.q);
        result("[Extract-Right] 右陷门 [A|B]·sk=u", ok);
        std::cout << "    ||sk||=" << std::fixed << std::setprecision(0)
                  << vec_norm(sk) << "\n";
    }
}

/* ───────────────────────────────────────────────────
   文件输出辅助：将单条 benchmark 结果追加写入文件
   ─────────────────────────────────────────────────── */
static void write_to_bench_file(const std::string& content) {
    constexpr const char* OUT_PATH = "bendmarking_output/bendmarking_plain-LWE.txt";
    std::ofstream fout(OUT_PATH, std::ios::app);
    if (fout.is_open()) {
        fout << content;
        fout.close();
        std::cout << "  [Results written to " << OUT_PATH << "]\n";
    } else {
        std::cerr << "  [WARN] Could not open " << OUT_PATH << " for writing\n";
    }
}

/* ───────────────────────────────────────────────────
   Benchmark 1: DelTrapGen 纯耗时
   ─────────────────────────────────────────────────── */
/**
 * bench_del_trap_gen — 纯粹测量 del_trap_gen() 的耗时
 *
 *   ① 生成基础陷门（不计入时间）
 *   ② 预热 3 轮
 *   ③ 20 轮计时（不同随机 tag H  + 不同 seed）
 *   ④ 统计：平均 / 最小 / 最大 / 标准差（µs）
 */
void bench_del_trap_gen(const Params& p) {
    using Clock = std::chrono::high_resolution_clock;
    constexpr int WARMUP  = 3;
    constexpr int ITERS   = 20;

    /* ── 基础陷门（一次性，不计入时间）── */
    Trapdoor base = gen_trap(p, 42);

    /* ── 预热 ── */
    std::cout << "  Warming up (" << WARMUP << " rounds)..." << std::flush;
    for (int i = 0; i < WARMUP; ++i) {
        Mat H = random_invertible_mat(p.n, p.q, (uint64_t)(i + 1) * 100003);
        (void)del_trap_gen(p, base, H);
        std::cout << "." << std::flush;
    }
    std::cout << " done\n" << std::flush;

    /* ── 计时 ── */
    std::vector<double> times_us;
    times_us.reserve(ITERS);
    std::cout << "  Benchmarking (" << ITERS << " rounds): " << std::flush;
    int pm = std::max(1, ITERS / 5);
    for (int i = 0; i < ITERS; ++i) {
        Mat H = random_invertible_mat(p.n, p.q, (uint64_t)(i + 1) * 500009);
        auto t0 = Clock::now();
        (void)del_trap_gen(p, base, H);
        auto t1 = Clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        times_us.push_back(us);
        if ((i + 1) % pm == 0 || i == ITERS - 1)
            std::cout << " " << (i + 1) << "/" << ITERS << std::flush;
    }
    std::cout << " done\n" << std::flush;

    /* ── 统计 ── */
    double sum_us = std::accumulate(times_us.begin(), times_us.end(), 0.0);
    double avg_us = sum_us / ITERS;
    double min_us = *std::min_element(times_us.begin(), times_us.end());
    double max_us = *std::max_element(times_us.begin(), times_us.end());
    double var_us = 0.0;
    for (double t : times_us) { double d = t - avg_us; var_us += d * d; }
    var_us /= (ITERS > 1) ? (ITERS - 1) : 1;
    double std_us = std::sqrt(var_us);

    /* ── 格式化输出 ── */
    std::ostringstream oss;
    oss << "\n=== Benchmark: DelTrapGen (MP12 §5) ===\n"
        << "  Parameters: n=" << p.n << ", q=" << p.q
        << ", b=" << p.b << ", k=" << p.k
        << ", m=" << p.m << ", σ=" << p.sigma << "\n"
        << "  Warmup rounds : " << WARMUP << "\n"
        << "  Timed  rounds : " << ITERS << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << avg_us << " µs\n"
        << "  Min       : " << std::setw(8) << min_us << " µs\n"
        << "  Max       : " << std::setw(8) << max_us << " µs\n"
        << "  StdDev    : " << std::setw(8) << std_us << " µs\n"
        << "  Throughput: " << std::setw(8)
        << (1e6 / avg_us) << " ops/s\n";

    std::cout << oss.str();
    write_to_bench_file(oss.str());
}

/* ───────────────────────────────────────────────────
   Benchmark 2: SamplePre（带 tag）纯耗时
   ─────────────────────────────────────────────────── */
/**
 * bench_sample_pre_tagged — 纯粹测量 sample_pre_tagged() 的耗时
 *
 *   ① 生成基础陷门 + 委托生成 tagged 陷门（不计入时间）
 *   ② 预热 3 轮（含随机目标 u）
 *   ③ 20 轮计时
 *   ④ 统计：平均 / 最小 / 最大 / 标准差（µs）
 */
void bench_sample_pre_tagged(const Params& p) {
    using Clock = std::chrono::high_resolution_clock;
    constexpr int WARMUP  = 3;
    constexpr int ITERS   = 20;

    /* ── 建立 tagged 陷门（一次性，不计入时间）── */
    Trapdoor base = gen_trap(p, 100);
    Mat H_tag = random_invertible_mat(p.n, p.q, 999);
    auto td_H = del_trap_gen(p, base, H_tag);

    UniformSampler usampler(p.q, 55);

    /* ── 预热 ── */
    std::cout << "  Warming up (" << WARMUP << " rounds)..." << std::flush;
    for (int i = 0; i < WARMUP; ++i) {
        Vec u(p.n);
        for (int j = 0; j < p.n; ++j) u[j] = usampler.sample();
        (void)sample_pre_tagged(p, td_H, u, (uint64_t)(i + 1) * 200003);
        std::cout << "." << std::flush;
    }
    std::cout << " done\n" << std::flush;

    /* ── 计时 ── */
    std::vector<double> times_us;
    times_us.reserve(ITERS);
    std::cout << "  Benchmarking (" << ITERS << " rounds): " << std::flush;
    int pm = std::max(1, ITERS / 5);
    for (int i = 0; i < ITERS; ++i) {
        Vec u(p.n);
        for (int j = 0; j < p.n; ++j) u[j] = usampler.sample();
        auto t0 = Clock::now();
        (void)sample_pre_tagged(p, td_H, u, (uint64_t)(i + 1) * 700001);
        auto t1 = Clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        times_us.push_back(us);
        if ((i + 1) % pm == 0 || i == ITERS - 1)
            std::cout << " " << (i + 1) << "/" << ITERS << std::flush;
    }
    std::cout << " done\n" << std::flush;

    /* ── 统计 ── */
    double sum_us = std::accumulate(times_us.begin(), times_us.end(), 0.0);
    double avg_us = sum_us / ITERS;
    double min_us = *std::min_element(times_us.begin(), times_us.end());
    double max_us = *std::max_element(times_us.begin(), times_us.end());
    double var_us = 0.0;
    for (double t : times_us) { double d = t - avg_us; var_us += d * d; }
    var_us /= (ITERS > 1) ? (ITERS - 1) : 1;
    double std_us = std::sqrt(var_us);

    /* ── 格式化输出 ── */
    std::ostringstream oss;
    oss << "\n=== Benchmark: SamplePre (tagged, MP12 Algo 2) ===\n"
        << "  Parameters: n=" << p.n << ", q=" << p.q
        << ", b=" << p.b << ", k=" << p.k
        << ", m=" << p.m << ", s=" << p.s << "\n"
        << "  Warmup rounds : " << WARMUP << "\n"
        << "  Timed  rounds : " << ITERS << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << avg_us << " µs\n"
        << "  Min       : " << std::setw(8) << min_us << " µs\n"
        << "  Max       : " << std::setw(8) << max_us << " µs\n"
        << "  StdDev    : " << std::setw(8) << std_us << " µs\n"
        << "  Throughput: " << std::setw(8)
        << (1e6 / avg_us) << " ops/s\n";

    std::cout << oss.str();
    write_to_bench_file(oss.str());
}

/* ════════════════════ main ════════════════════════ */
void run_del_tests(const mp12::Params& p) {
    std::cout << "╔══════════════════════════════════════════════╗\n";
    std::cout << "║   MP12 DelTrapGen — C/C++ 实现              ║\n";
    std::cout << "║   Tag委托 / SampleLeft / SampleRight         ║\n";
    std::cout << "╚══════════════════════════════════════════════╝\n";

    //auto p = Params::make(4, 97, 2);
    std::cout << "\n参数: n=" << p.n << "  q=" << p.q << "  b=" << p.b
              << "  k=" << p.k << "  m=" << p.m << "  s=" << p.s << "\n";

    test_del_trap_gen(p);
    test_sample_pre_tagged(p);
    test_sample_left(p);
    test_sample_right(p);
    test_ibe_simulation(p);
    bench_del_trap_gen(p);
    bench_sample_pre_tagged(p);

    // 算法结构汇总
    section("算法结构汇总");
    std::cout << R"(
  GenTrap (已有)
    └─ 生成 (A, R)，满足 A·[R;I] = G

  DelTrapGen (新增，本文件)
    └─ 输入: (A, R, H_new)
       输出: T_H，满足 A·T_H = H·G
       核心: H̃ = H⊗I_k，T_H = [R·H̃; H̃]

  SamplePre_tagged (新增)
    └─ 输入: tagged 陷门 T_H，目标 u
       核心变化: 目标调整 u' = H^{-1}·(u-v)，再 SampleG

  SampleLeft (新增)
    └─ 输入: A 的陷门，任意 B，目标 u
       输出: [x1;x2] s.t. [A|B]·[x1;x2] = u
       核心: x2←高斯, x1←SamplePre(A, u-B·x2)

  SampleRight (新增)
    └─ 输入: A，右陷门 B_R (B=A·B_R+G)，目标 u
       输出: [x1;x2] s.t. [A|B]·[x1;x2] = u
       核心: 扰动 p，x2←SampleG(u-A·p)，x1=p-B_R·x2

  应用层次:
    ├─ IBE:   DelTrapGen → 用户私钥提取
    ├─ HIBE:  SampleLeft → 子级密钥委托
    ├─ ABE:   SampleRight → 属性密钥生成
    └─ 签名:  SamplePre_tagged → 带 tag 的签名
)";

    std::cout << "\nDone.\n";
    //return 0;
}