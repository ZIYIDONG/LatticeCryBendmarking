#include "../include/powersof.h"
#include <iostream>
#include <iomanip>
#include <random>
#include "unified_params.h"

using namespace cryptolib;

static void print_vec(const char* label, const Vec& v) {
    std::cout << label << " [";
    for (size_t i = 0; i < v.size(); i++)
        std::cout << (i ? ", " : "") << v[i];
    std::cout << "]\n";
}

int main() {
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
    std::cout << "\n--- Test 2: Powersof3 标量 (q=97, b=3) ---\n";
    {
        auto __u_p = unified::default_mp12_params(); long q = __u_p.q;
        int b = 3;
        Vec p = powers_of_b_scalar(5, b, q);
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
        bool all_ok = true;
        for (long a = 0; a < q; a++) {
            Vec d = bit_decomp_scalar(a, b, q);
            // 重构: Σ d_j · b^j
            long recon = 0, pw = 1;
            for (int j = 0; j < k; j++) {
                recon += d[j] * pw;
                pw *= b;
            }
            if (mod_q(recon, q) != a) { all_ok = false; break; }
        }
        std::cout << "  全部 " << q << " 个值重构: "
                  << (all_ok ? "PASS" : "FAIL") << "\n";
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
        Case cases[] = {{97, 2}, {97, 3}, {97, 5}, {1031, 2},
                        {1031, 4}, {1031, 8}, {65537, 16}};
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

    std::cout << "\nDone.\n";
    return 0;
}
