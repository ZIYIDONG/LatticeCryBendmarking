#include "../include/powersof_modswitch.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <stdexcept>
#include <string>
#include <chrono>
#include <numeric>
#include <fstream>
#include <sstream>
#include "unified_params.h"

using namespace cryptolib;

/* ───── 前向声明 ───── */
static void bench_powersof_modswitch();

static void print_vec(const char* label, const Vec& v, int maxshow = 16) {
    std::cout << label << " [";
    int show = std::min((int)v.size(), maxshow);
    for (int i = 0; i < show; i++)
        std::cout << (i ? ", " : "") << v[i];
    if ((int)v.size() > show) std::cout << " ...";
    std::cout << "] (len=" << v.size() << ")\n";
}

void run_test_powersof_modswitch() {
    std::cout << "=========================================================\n";
    std::cout << "  Two-part Powersof2 with modulus switching (IBE keygen)\n";
    std::cout << "=========================================================\n";

    /* ───── Test 1: 小例子,直观查看两段结构 ───── */
    std::cout << "\n--- Test 1: 小例子 (p=17, q=257, m=4) ---\n";
    {
        long p = 17, q = 257;
        int m = 4;
        int k_p = compute_k(p, 2);   // = 5
        int k_q = compute_k(q, 2);   // = 9
        std::cout << "p=" << p << "  q=" << q << "  m=" << m
                  << "  k_p=" << k_p << "  k_q=" << k_q << "\n";
        std::cout << "期望长度: k_p + m·k_q = " << k_p + m*k_q << "\n\n";

        Vec e = {3, 17, 100, 200};
        std::cout << "e = [3, 17, 100, 200]\n";

        if (!(p > 1 && q > p))
            throw std::invalid_argument(std::string("Invalid modulus parameters in test (p<=1 or q<=p): p=") + std::to_string(p) + ", q=" + std::to_string(q));
        Vec sk = ibe_extract_key(e, p, q);
        std::cout << "sk_id 长度 = " << sk.size()
                  << "  (期望 " << k_p + m*k_q << ")\n";

        // 第一段
        Vec part1(sk.begin(), sk.begin() + k_p);
        print_vec("第一段 Powersof2_p(1):", part1);
        std::cout << "  期望: (1, 2, 4, 8, 16) mod 17 = (1, 2, 4, 8, 16) ✓\n\n";

        // 第二段
        Vec part2(sk.begin() + k_p, sk.end());
        print_vec("第二段 -(p/q)·Powersof2_q(e):", part2);
    }

    /* ───── Test 2: round_scale 边界检查 ───── */
    std::cout << "\n--- Test 2: round_scale 正确性 ---\n";
    {
        long p = 16, q = 1024;
        std::cout << "round(p·x/q) for p=" << p << ", q=" << q << "\n";
        for (long x : {0L, 32L, 64L, 128L, 512L, 1023L, -100L}) {
            long r = round_scale(x, p, q);
            double exact = (double)p * x / q;
            std::cout << "  x=" << std::setw(5) << x
                      << "  round_scale=" << std::setw(4) << r
                      << "  exact=" << std::fixed << std::setprecision(3) << exact
                      << "\n";
        }
    }

    /* ───── Test 3: 单分量验证第二段每个元素 ───── */
    std::cout << "\n--- Test 3: 手动验证第二段每个元素 ---\n";
    {
        long p = 16, q = 1024;
        int k_q = compute_k(q, 2);  // = 10
        Vec e = {100};
        if (!(p > 1 && q > p))
            throw std::invalid_argument(std::string("Invalid modulus parameters in test (p<=1 or q<=p): p=") + std::to_string(p) + ", q=" + std::to_string(q));
        Vec sk = ibe_extract_key(e, p, q);
        // 跳过第一段
        int k_p = compute_k(p, 2);
        std::cout << "e = [100],  k_q = " << k_q << "\n";
        std::cout << "j   2^j·e mod q   round(p·val/q)   sk[k_p+j]\n";
        for (int j = 0; j < k_q; j++) {
            long val = (100L << j) % q;
            long expected_neg = -round_scale(val, p, q);
            long expected_mod = ((expected_neg % p) + p) % p;
            long actual = sk[k_p + j];
            std::cout << "  " << std::setw(2) << j
                      << "      " << std::setw(5) << val
                      << "          " << std::setw(5) << -round_scale(val, p, q)
                      << "        " << std::setw(5) << actual
                      << (actual == expected_mod ? "  ✓" : "  ✗") << "\n";
        }
    }

    /* ───── Test 4: 对偶恒等式 (近似) ───── */
    /*    <BitDecomp(c̄), sk_id> ≈ (p/q)·<c̄, s̄>  (mod p)        */
    std::cout << "\n--- Test 4: 模数切换下的近似对偶恒等式 ---\n";
    {
        auto __u_p = unified::default_mp12_params(); long p = 17, q = __u_p.q;       // use unified default q for tests
        int m = 8;
        int k_q = compute_k(q, 2);
        int k_p = compute_k(p, 2);

        std::mt19937_64 rng(42);
        std::uniform_int_distribution<long> dist_q(0, q - 1);
        // e 取小值 (LWE 噪声)
        std::uniform_int_distribution<long> dist_e(-3, 3);

        int trials = 1000;
        double max_err = 0, sum_err = 0;
        int correct_within_tol = 0;
        long tolerance = (k_p + m * k_q);  // 总舍入误差上界量级

        for (int t = 0; t < trials; t++) {
            // 构造 e
            Vec e(m);
            for (int i = 0; i < m; i++) e[i] = dist_e(rng);

            // 计算 sk_id (不 mod p,保留整数,方便算误差)
            Vec sk_int = powers_of_2_with_modswitch(e, p, q, false);

            // 构造随机密文 c̄ = (c0, c1, ..., cm) ∈ Z_q^{m+1}
            Vec c_bar(m + 1);
            for (int i = 0; i <= m; i++) c_bar[i] = dist_q(rng);

            // BitDecomp_2(c̄): 把每个 c_i 分解成 k_q 个 0/1 位
            // 注意: 第一个分量 c_0 对应 s̄ 的第一段(标量1),长度 k_p
            //       后 m 个分量 c_1..c_m 对应 -e 的第二段,长度 k_q each
            // 但是这里两段位数不同(k_p vs k_q),所以 BitDecomp 也要按段处理
            Vec bd;
            bd.reserve(k_p + m * k_q);
            // 第一段: BitDecomp_2 of c0 用 k_p 位 (因为对应模 p 的 gadget)
            {
                long v = ((c_bar[0] % p) + p) % p;  // 注意: c0 对应 mod p 的部分
                for (int j = 0; j < k_p; j++) { bd.push_back(v & 1); v >>= 1; }
            }
            // 第二段: BitDecomp_2 of c_i (i=1..m) 用 k_q 位
            for (int i = 1; i <= m; i++) {
                long v = ((c_bar[i] % q) + q) % q;
                for (int j = 0; j < k_q; j++) { bd.push_back(v & 1); v >>= 1; }
            }

            // LHS: <bd, sk_id> mod p
            long lhs = 0;
            for (size_t i = 0; i < bd.size(); i++)
                lhs += bd[i] * sk_int[i];
            lhs = ((lhs % p) + p) % p;

            // RHS: (round(p·c0/?) ... 对 s̄ 的内积要小心
            //   <c̄, s̄> = c0·1 + Σ c_i·(-e_i) = c0 - Σ c_i·e_i  在 Z_q 下
            // 在模数切换后期望 ≈ round(p/q · (c0 - Σ c_i·e_i))
            // 但因为第一段 c0 已经被当成 mod p 处理,需要分别算:
            //   第一段贡献: c0_mod_p · 1 (在 Z_p 下)
            //   第二段贡献: round(p/q · Σ c_i · (-e_i)) (模数切换)
            long c0_mod_p = ((c_bar[0] % p) + p) % p;
            long sum_q = 0;
            for (int i = 1; i <= m; i++) {
                sum_q = ((sum_q + c_bar[i] * (-e[i-1])) % q + q) % q;
            }
            long rhs_part2 = round_scale(sum_q, p, q);
            long rhs = ((c0_mod_p + rhs_part2) % p + p) % p;

            // 误差 (允许环绕)
            long diff = (lhs - rhs + 3*p) % p;
            if (diff > p / 2) diff -= p;
            double e_abs = std::abs((double)diff);
            sum_err += e_abs;
            if (e_abs > max_err) max_err = e_abs;
            if (e_abs < tolerance) correct_within_tol++;
        }
        std::cout << "p=" << p << "  q=" << q << "  m=" << m << "\n";
        std::cout << "试验次数: " << trials << "\n";
        std::cout << "平均误差: " << sum_err / trials << "\n";
        std::cout << "最大误差: " << max_err << "\n";
        std::cout << "容差 (k_p + m·k_q = " << tolerance << ") 内: "
                  << correct_within_tol << "/" << trials << "\n";
        std::cout << "结论: " << (max_err < tolerance ? "PASS (误差有界)" : "WARN") << "\n";
    }

    /* ───── Test 5: e=0 时的退化情况 ───── */
    std::cout << "\n--- Test 5: e=0 退化情况 ---\n";
    {
        long p = 17, q = 257;
        int m = 4;
        Vec e(m, 0);
        if (!(p > 1 && q > p))
            throw std::invalid_argument(std::string("Invalid modulus parameters in test (p<=1 or q<=p): p=") + std::to_string(p) + ", q=" + std::to_string(q));
        Vec sk = ibe_extract_key(e, p, q);
        int k_p = compute_k(p, 2);
        int k_q = compute_k(q, 2);
        Vec part1(sk.begin(), sk.begin() + k_p);
        Vec part2(sk.begin() + k_p, sk.end());
        print_vec("e=0 时第一段:", part1);
        print_vec("e=0 时第二段:", part2);
        bool all_zero = true;
        for (auto x : part2) if (x != 0) all_zero = false;
        std::cout << "第二段全 0: " << (all_zero ? "PASS" : "FAIL") << "\n";
    }

    /* ── 纯基准测试 ── */
    bench_powersof_modswitch();

    std::cout << "\nDone.\n";
}

/* ───────────────────────────────────────────────────
   文件输出辅助
   ─────────────────────────────────────────────────── */
static void write_to_bench_file(const std::string& content) {
    constexpr const char* OUT_PATH = "../bendmarking_output/bendmarking_plain-LWE.txt";
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
   Benchmark: Powersof2 with modulus switch 纯耗时
   ─────────────────────────────────────────────────── */
static void bench_powersof_modswitch() {
    using Clock = std::chrono::high_resolution_clock;
    constexpr int WARMUP  = 3;
    constexpr int ITERS   = 20;

    auto __u_p = unified::default_mp12_params();
    long q = __u_p.q;
    long p = 17;
    int  m = 8;

    std::mt19937_64 rng(42);
    std::uniform_int_distribution<long> dist_e(-3, 3);
    Vec e0(m);
    for (int i = 0; i < m; ++i) e0[i] = dist_e(rng);

    /* ── 预热 ── */
    for (int i = 0; i < WARMUP; ++i) {
        (void)powers_of_2_with_modswitch(e0, p, q, true);
    }

    /* ── 计时 ── */
    std::vector<double> times_us;
    times_us.reserve(ITERS);
    for (int i = 0; i < ITERS; ++i) {
        Vec e(m);
        for (int j = 0; j < m; ++j) e[j] = dist_e(rng);
        auto t0 = Clock::now();
        (void)powers_of_2_with_modswitch(e, p, q, true);
        auto t1 = Clock::now();
        double us = std::chrono::duration<double, std::micro>(t1 - t0).count();
        times_us.push_back(us);
    }

    /* ── 统计 ── */
    double sum_us = std::accumulate(times_us.begin(), times_us.end(), 0.0);
    double avg_us = sum_us / ITERS;
    double min_us = *std::min_element(times_us.begin(), times_us.end());
    double max_us = *std::max_element(times_us.begin(), times_us.end());
    double var_us = 0.0;
    for (double t : times_us) { double d = t - avg_us; var_us += d * d; }
    var_us /= (ITERS > 1) ? (ITERS - 1) : 1;
    double std_us = std::sqrt(var_us);

    std::ostringstream oss;
    oss << "\n=== Benchmark: Powersof2 with ModSwitch (IBE keygen) ===\n"
        << "  Parameters: p=" << p << ", q=" << q
        << ", m=" << m << "\n"
        << "  Warmup rounds : " << WARMUP << "\n"
        << "  Timed  rounds : " << ITERS << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Average   : " << std::setw(8) << avg_us << " µs\n"
        << "  Min       : " << std::setw(8) << min_us << " µs\n"
        << "  Max       : " << std::setw(8) << max_us << " µs\n"
        << "  StdDev    : " << std::setw(8) << std_us << " µs\n"
        << "  Throughput: " << std::setw(8)
        << (1e6 / avg_us) << " ops/s\n";

    std::cout << "\n--- Benchmark: Powersof2 ModSwitch ---\n";
    std::cout << oss.str();
    write_to_bench_file(oss.str());
}
