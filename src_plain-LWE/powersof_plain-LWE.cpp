#include "../include_plain-LWE/powersof_plain-LWE.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include <numeric>
#include <fstream>
#include <sstream>
#include "unified_params_plain-LWE.h"

using namespace cryptolib;

/* ───── 前向声明 ───── */
static void bench_powersof();

static void print_vec(const char* label, const Vec& v) {
    std::cout << label << " [";
    for (size_t i = 0; i < v.size(); i++)
        std::cout << (i ? ", " : "") << v[i];
    std::cout << "]\n";
}

void run_test_powersof() {
    std::cout << "================================================\n";
    std::cout << "  Powersof_b / BitDecomp_b  测试\n";
    std::cout << "================================================\n";

    /* ───── Test 1: 标量 Powersof2,小例子 ───── */
    std::cout << "\n--- Test 1: Powersof2 标量 (q=17, b=2) ---\n";
    {
        auto __u_p = unified::default_mp12_params(); long q = __u_p.q;
        int b = 2;
        int k = compute_k(q, b);
        std::cout << "q=" << q << "  b=" << b << "  k=" << k << "\n";
        for (long a : {1L, 3L, 5L, 7L}) {
            Vec p = powers_of_b_scalar(a, b, q);
            std::cout << "  Powersof2(" << a << ") = ";
            print_vec("", p);
        }
    }

    /* ───── Test 2: 标量 Powers-of-3 ───── */
    {
        auto __u_p = unified::default_mp12_params(); long q = __u_p.q;
        int b = 3;
        Vec p = powers_of_b_scalar(5, b, q);
        std::cout << "\n--- Test 2: Powersof3 标量 (q=" << q << ", b=3) ---\n";
        std::cout << "  Powersof3(5)  = ";
        print_vec("", p);
        // 期望: (5, 15, 45, 38, 17, 51)  (因为 5·81 = 405 = 4·97+17 = 17)
    }

    /* ───── Test 3: BitDecomp 是 Powersof 的逆向操作 ───── */
    std::cout << "\n--- Test 3: BitDecomp 重构原值 ---\n";
    {
        auto __u_p = unified::default_mp12_params(); long q = __u_p.q;
        int b = 2;
        int k = compute_k(q, b);
        std::cout << "q=" << q << "  k=" << k << "\n";
        long check_count = (q > 10000) ? 10000 : q;
        bool all_ok = true;
        for (long a = 0; a < check_count; a++) {
            Vec d = bit_decomp_scalar(a, b, q);
            long recon = 0, pw = 1;
            for (int j = 0; j < k; j++) { recon += d[j] * pw; pw *= b; }
            if (mod_q(recon, q) != a) { all_ok = false; break; }
        }
        if (all_ok && q > check_count) {
            std::mt19937_64 rng2(7);
            std::uniform_int_distribution<long> dist2(0, q - 1);
            for (int t = 0; t < 1000 && all_ok; t++) {
                long a = dist2(rng2);
                Vec d = bit_decomp_scalar(a, b, q);
                long recon = 0, pw = 1;
                for (int j = 0; j < k; j++) { recon += d[j] * pw; pw *= b; }
                if (mod_q(recon, q) != a) all_ok = false;
            }
        }
        std::cout << "  " << (q > check_count ? std::to_string(check_count) + " + 1000 spot"
                                              : std::to_string(q))
                  << " values: " << (all_ok ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 4: 关键对偶恒等式 ───── */
    /*    <BitDecomp(x), Powersof(y)> == x·y (mod q)         */
    std::cout << "\n--- Test 4: 对偶恒等式 (核心性质) ---\n";
    {
        auto __u_p = unified::default_mp12_params(); long q = __u_p.q;          // 素数
        int b = 2;
        std::mt19937_64 rng(42);
        std::uniform_int_distribution<long> dist(0, q - 1);
        int trials = 10000, pass = 0;
        for (int t = 0; t < trials; t++) {
            long x = dist(rng), y = dist(rng);
            Vec dx = bit_decomp_scalar(x, b, q);
            Vec py = powers_of_b_scalar(y, b, q);
            long lhs = inner_prod_mod(dx, py, q);
            long rhs = mod_q(x * y, q);
            if (lhs == rhs) pass++;
        }
        std::cout << "  <BitDecomp(x), Powersof(y)> = x·y (mod q)\n";
        std::cout << "  随机测试 " << pass << "/" << trials
                  << "  " << (pass == trials ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 5: 多种 (q, b) 组合 ───── */
    std::cout << "\n--- Test 5: 不同基 b 的对偶恒等式 ---\n";
    {
        std::mt19937_64 rng(123);
        struct Case { long q; int b; };
        auto __u_p = unified::default_mp12_params();
        Case cases[] = {{__u_p.q, 2}, {__u_p.q, 3}, {__u_p.q, 5}, {__u_p.q, 2},
                        {__u_p.q, 4}, {__u_p.q, 8}, {__u_p.q, 16}};
        for (auto c : cases) {
            std::uniform_int_distribution<long> dist(0, c.q - 1);
            int pass = 0, trials = 1000;
            for (int t = 0; t < trials; t++) {
                long x = dist(rng), y = dist(rng);
                Vec dx = bit_decomp_scalar(x, c.b, c.q);
                Vec py = powers_of_b_scalar(y, c.b, c.q);
                if (inner_prod_mod(dx, py, c.q) == mod_q(x * y, c.q))
                    pass++;
            }
            std::cout << "  q=" << std::setw(6) << c.q
                      << "  b=" << std::setw(2) << c.b
                      << "  k=" << std::setw(2) << compute_k(c.q, c.b)
                      << "  " << pass << "/" << trials
                      << "  " << (pass == trials ? "PASS" : "FAIL") << "\n";
        }
    }

    /* ───── Test 6: 向量版 + 平衡分解 ───── */
    std::cout << "\n--- Test 6: 向量版对偶恒等式 ---\n";
    {
        auto __u_p = unified::default_mp12_params(); long q = __u_p.q;
        int b = 2;
        int n = 8;
        std::mt19937_64 rng(7);
        std::uniform_int_distribution<long> dist(0, q - 1);

        Vec v1(n), v2(n);
        for (int i = 0; i < n; i++) { v1[i] = dist(rng); v2[i] = dist(rng); }

        Vec p1 = powers_of_b_vec(v1, b, q);
        Vec d2 = bit_decomp_vec(v2, b, q);

        // <BitDecomp(v2), Powersof(v1)> = <v2, v1>
        long lhs = inner_prod_mod(d2, p1, q);
        long rhs = 0;
        for (int i = 0; i < n; i++) rhs = mod_q(rhs + v1[i] * v2[i], q);
        std::cout << "  <BitDecomp(v2), Powersof(v1)> = " << lhs << "\n";
        std::cout << "  <v1, v2> mod q             = " << rhs << "\n";
        std::cout << "  " << (lhs == rhs ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 7: 平衡分解的范数优势 ───── */
    std::cout << "\n--- Test 7: 平衡 vs 非平衡分解 (b=8) ---\n";
    {
        auto __u_p = unified::default_mp12_params(); long q = __u_p.q;
        int b = 8;
        std::mt19937_64 rng(99);
        std::uniform_int_distribution<long> dist(0, q - 1);
        double sum_unbalanced = 0, sum_balanced = 0;
        int trials = 1000;
        for (int t = 0; t < trials; t++) {
            long a = dist(rng);
            Vec u = bit_decomp_scalar(a, b, q);
            Vec ba = bit_decomp_balanced(a, b, q);
            for (auto x : u)  sum_unbalanced += (double)x * x;
            for (auto x : ba) sum_balanced   += (double)x * x;
        }
        std::cout << "  非平衡 (∈[0,b))   平均 ‖·‖²: "
                  << sum_unbalanced / trials << "\n";
        std::cout << "  平衡   (∈[-b/2,b/2)) 平均 ‖·‖²: "
                  << sum_balanced / trials << "\n";
        std::cout << "  减小比例: " << std::fixed << std::setprecision(1)
                  << (1 - sum_balanced / sum_unbalanced) * 100 << "%\n";
    }

    /* ── 纯基准测试 ── */
    bench_powersof();

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
   Benchmark: Powersof / BitDecomp 纯耗时
   ─────────────────────────────────────────────────── */
static void bench_powersof() {
    using Clock = std::chrono::high_resolution_clock;
    constexpr int WARMUP  = 3;
    constexpr int ITERS   = 20;

    auto __u_p = unified::default_mp12_params();
    long q = __u_p.q;
    int  b = __u_p.b;
    int  k = compute_k(q, b);
    int  n = 8;

    std::mt19937_64 rng(0);
    std::uniform_int_distribution<long> dist(0, q - 1);

    Vec v(n);
    for (int i = 0; i < n; ++i) v[i] = dist(rng);

    /* ── 预热 ── */
    for (int i = 0; i < WARMUP; ++i) {
        (void)powers_of_b_vec(v, b, q);
        (void)bit_decomp_vec(v, b, q);
    }

    /* ── Powersof 计时 ── */
    std::vector<double> times_po, times_bd;
    times_po.reserve(ITERS);
    times_bd.reserve(ITERS);
    for (int i = 0; i < ITERS; ++i) {
        Vec vi(n);
        for (int j = 0; j < n; ++j) vi[j] = dist(rng);

        auto t0 = Clock::now();
        (void)powers_of_b_vec(vi, b, q);
        auto t1 = Clock::now();
        times_po.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());

        auto t2 = Clock::now();
        (void)bit_decomp_vec(vi, b, q);
        auto t3 = Clock::now();
        times_bd.push_back(std::chrono::duration<double, std::micro>(t3 - t2).count());
    }

    auto stats = [&](const std::vector<double>& tv) {
        double s = std::accumulate(tv.begin(), tv.end(), 0.0);
        double a = s / ITERS;
        double mn = *std::min_element(tv.begin(), tv.end());
        double mx = *std::max_element(tv.begin(), tv.end());
        double vv = 0;
        for (double t : tv) { double d = t - a; vv += d * d; }
        vv /= (ITERS > 1) ? (ITERS - 1) : 1;
        return std::make_tuple(a, mn, mx, std::sqrt(vv));
    };

    auto [avg_po, min_po, max_po, std_po] = stats(times_po);
    auto [avg_bd, min_bd, max_bd, std_bd] = stats(times_bd);

    std::ostringstream oss;
    oss << "\n=== Benchmark: Powersof / BitDecomp ===\n"
        << "  Parameters: n=" << n << ", q=" << q
        << ", b=" << b << ", k=" << k << "\n"
        << "  Warmup rounds : " << WARMUP << "\n"
        << "  Timed  rounds : " << ITERS << "\n\n"
        << std::fixed << std::setprecision(1)
        << "  Powersof_b     Avg=" << std::setw(8) << avg_po
        << " µs  Min=" << std::setw(8) << min_po
        << " µs  Max=" << std::setw(8) << max_po
        << " µs  σ=" << std::setw(8) << std_po << " µs\n"
        << "  BitDecomp_b   Avg=" << std::setw(8) << avg_bd
        << " µs  Min=" << std::setw(8) << min_bd
        << " µs  Max=" << std::setw(8) << max_bd
        << " µs  σ=" << std::setw(8) << std_bd << " µs\n";

    std::cout << "\n--- Benchmark: Powersof / BitDecomp ---\n";
    std::cout << oss.str();
    write_to_bench_file(oss.str());
}
