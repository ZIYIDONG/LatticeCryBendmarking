#include "../include/frd.h"
#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include "unified_params.h"

using namespace cryptolib;

/* ───── 矩阵打印 ───── */
static void print_mat(const char* label, const Mat& M, int q = 0) {
    std::cout << label << "  (" << M.size() << "×" << M[0].size() << ")\n";
    for (size_t i = 0; i < M.size(); i++) {
        std::cout << "  [";
        for (size_t j = 0; j < M[i].size(); j++)
            std::cout << std::setw(5) << M[i][j];
        std::cout << " ]\n";
    }
}

/* ───── 高斯消元求矩阵秩 (mod q) ───── */
static int matrix_rank_mod(Mat M, long q) {
    int rows = (int)M.size(), cols = (int)M[0].size();
    int rank = 0;
    int col = 0;
    for (int r = 0; r < rows && col < cols; ) {
        // 找列 col 中第 r 行起的非零主元
        int pivot = -1;
        for (int i = r; i < rows; i++)
            if (M[i][col] != 0) { pivot = i; break; }
        if (pivot < 0) { col++; continue; }
        std::swap(M[r], M[pivot]);
        long inv = mod_inv_q(M[r][col], q);
        for (int j = col; j < cols; j++)
            M[r][j] = mod_pos(M[r][j] * inv, q);
        for (int i = 0; i < rows; i++) {
            if (i == r) continue;
            long fac = M[i][col];
            if (fac == 0) continue;
            for (int j = col; j < cols; j++)
                M[i][j] = mod_pos(M[i][j] - fac * M[r][j], q);
        }
        rank++;
        r++;
        col++;
    }
    return rank;
}

/* ───── 矩阵相减 mod q ───── */
static Mat mat_sub(const Mat& A, const Mat& B, long q) {
    int r = (int)A.size(), c = (int)A[0].size();
    Mat R(r, Vec(c, 0));
    for (int i = 0; i < r; i++)
        for (int j = 0; j < c; j++)
            R[i][j] = mod_pos(A[i][j] - B[i][j], q);
    return R;
}

int main() {
    std::cout << "================================================\n";
    std::cout << "  FRD: Full-Rank Difference Encoding (ABB10)\n";
    std::cout << "================================================\n";

    /* ───── Test 1: 找不可约多项式 ───── */
    std::cout << "\n--- Test 1: 寻找 F_q 上 n 次不可约多项式 ---\n";
    {
        struct Case { int n; long q; };
        Case cases[] = {{2, 7}, {3, 11}, {4, 17}, {5, 23}, {6, 31}, {4, 97}};
        for (auto c : cases) {
            Vec f = find_irreducible(c.n, c.q, 42);
            std::cout << "  n=" << c.n << "  q=" << std::setw(3) << c.q
                      << "  f(x) = ";
            for (int i = (int)f.size() - 1; i >= 0; i--) {
                if (f[i] == 0) continue;
                if (i < (int)f.size() - 1) std::cout << " + ";
                if (i == 0 || f[i] != 1) std::cout << f[i];
                if (i > 0) std::cout << "x";
                if (i > 1) std::cout << "^" << i;
            }
            bool ok = is_irreducible(f, c.q);
            std::cout << "   " << (ok ? "✓ irreducible" : "✗ FAIL") << "\n";
        }
    }

    /* ───── Test 2: FRD 输出维度和示例 ───── */
    std::cout << "\n--- Test 2: FRD 输出示例 (n=4, q=17) ---\n";
    {
        auto ctx = FRDContext::setup(4, 17, 7);
        std::cout << "  使用不可约多项式 f(x) = ";
        for (int i = (int)ctx.f.size() - 1; i >= 0; i--) {
            if (ctx.f[i] == 0) continue;
            if (i < (int)ctx.f.size() - 1) std::cout << " + ";
            if (i == 0 || ctx.f[i] != 1) std::cout << ctx.f[i];
            if (i > 0) std::cout << "x";
            if (i > 1) std::cout << "^" << i;
        }
        std::cout << "\n\n";

        Vec id1 = {3, 1, 4, 1};
        Mat H1 = frd_encode(ctx, id1);
        print_mat("FRD([3,1,4,1])", H1, 17);

        Vec id2 = {0, 0, 0, 0};
        Mat H0 = frd_encode(ctx, id2);
        std::cout << "\nFRD(0) 应该是零矩阵:\n";
        print_mat("FRD([0,0,0,0])", H0, 17);
    }

    /* ───── Test 3: 全秩差分核心性质 ───── */
    std::cout << "\n--- Test 3: 全秩差分性质 (核心) ---\n";
    {
        auto __u_p = unified::default_mp12_params_128(); long q = __u_p.q;
        int n = 4;
        auto ctx = FRDContext::setup(n, q, 13);

        std::mt19937_64 rng(2025);
        std::uniform_int_distribution<long> dist(0, q - 1);

        int trials = 500;
        int full_rank = 0;
        for (int t = 0; t < trials; t++) {
            Vec id1(n), id2(n);
            for (int i = 0; i < n; i++) {
                id1[i] = dist(rng);
                id2[i] = dist(rng);
            }
            // 确保 id1 ≠ id2
            if (id1 == id2) { id2[0] = (id2[0] + 1) % q; }

            Mat H1 = frd_encode(ctx, id1);
            Mat H2 = frd_encode(ctx, id2);
            Mat D  = mat_sub(H1, H2, q);
            int r = matrix_rank_mod(D, q);
            if (r == n) full_rank++;
        }
        std::cout << "  q=" << q << "  n=" << n << "  trials=" << trials << "\n";
        std::cout << "  H_id1 - H_id2 满秩: " << full_rank << "/" << trials
                  << "   " << (full_rank == trials ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 4: 线性性 FRD(a) - FRD(b) = FRD(a-b) ───── */
    std::cout << "\n--- Test 4: 线性性 FRD(a) - FRD(b) = FRD(a-b) ---\n";
    {
        auto __u_p = unified::default_mp12_params_128(); long q = __u_p.q;
        int n = 5;
        auto ctx = FRDContext::setup(n, q, 99);
        std::mt19937_64 rng(7);
        std::uniform_int_distribution<long> dist(0, q - 1);

        int pass = 0, trials = 200;
        for (int t = 0; t < trials; t++) {
            Vec a(n), b(n);
            for (int i = 0; i < n; i++) { a[i] = dist(rng); b[i] = dist(rng); }
            Vec ab(n);
            for (int i = 0; i < n; i++) ab[i] = mod_pos(a[i] - b[i], q);

            Mat lhs = mat_sub(frd_encode(ctx, a), frd_encode(ctx, b), q);
            Mat rhs = frd_encode(ctx, ab);
            bool eq = true;
            for (int i = 0; i < n && eq; i++)
                for (int j = 0; j < n && eq; j++)
                    if (lhs[i][j] != rhs[i][j]) eq = false;
            if (eq) pass++;
        }
        std::cout << "  线性性测试: " << pass << "/" << trials
                  << "  " << (pass == trials ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 5: 零身份 → 零矩阵 ───── */
    std::cout << "\n--- Test 5: FRD(0) = 0 ---\n";
    {
        auto __u_p = unified::default_mp12_params_128(); long q = __u_p.q;
        int n = 6;
        auto ctx = FRDContext::setup(n, q, 1);
        Vec zero(n, 0);
        Mat H = frd_encode(ctx, zero);
        bool all_zero = true;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (H[i][j] != 0) all_zero = false;
        std::cout << "  " << (all_zero ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 6: FRD(单位元 e_0) = I_n ───── */
    std::cout << "\n--- Test 6: FRD([1,0,0,...]) = I_n ---\n";
    {
        auto __u_p = unified::default_mp12_params_128(); long q = __u_p.q;
        int n = 5;
        auto ctx = FRDContext::setup(n, q, 5);
        Vec e0(n, 0); e0[0] = 1;     // 多项式 a(x) = 1
        Mat H = frd_encode(ctx, e0);
        bool is_id = true;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                long expected = (i == j) ? 1 : 0;
                if (H[i][j] != expected) is_id = false;
            }
        std::cout << "  乘以 1 应得单位矩阵: "
                  << (is_id ? "PASS" : "FAIL") << "\n";
    }

    /* ───── Test 7: 性能测试 ───── */
    std::cout << "\n--- Test 7: 性能 (n=8, q=8209) ---\n";
    {
        auto __u_p = unified::default_mp12_params_128(); long q = __u_p.q;
        int n = 8;
        auto ctx = FRDContext::setup(n, q, 11);
        std::mt19937_64 rng(0);
        std::uniform_int_distribution<long> dist(0, q - 1);

        auto t0 = std::chrono::high_resolution_clock::now();
        int N = 1000;
        long acc = 0;
        for (int t = 0; t < N; t++) {
            Vec id(n);
            for (int i = 0; i < n; i++) id[i] = dist(rng);
            Mat H = frd_encode(ctx, id);
            acc += H[0][0];
        }
        auto t1 = std::chrono::high_resolution_clock::now();
        double ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        std::cout << "  " << N << " 次 FRD 编码: " << std::fixed
                  << std::setprecision(2) << ms << " ms\n";
        std::cout << "  平均每次: " << (ms / N) * 1000 << " µs\n";
    }

    std::cout << "\nDone.\n";
    return 0;
}
