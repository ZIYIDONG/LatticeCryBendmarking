/**
 * LatticeCryBenchmarking — 主入口
 *
 * 统一运行所有 MP12 陷门、委托陷门、以及各模块的测试/基准。
 *
 * 构建：在项目根目录执行
 *   mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && cmake --build .
 *
 * 运行:
 *   build_<level>/LatticeCryBenchmarking               # 完整运行
 *   build_<level>/LatticeCryBenchmarking --no-decrypt   # 跳过解密（L1 加速）
 *
 * @author Ziyi Dong, 2026
 */

#include "../include_plain-LWE/mp12_plain-LWE.h"
#include "../include_plain-LWE/mp12deltrapgen_plain-LWE.h"
#include "../include_plain-LWE/mp12trapgen_plain-LWE.h"
#include "../include_plain-LWE/unified_params_plain-LWE.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <cstring>
#include <filesystem>

int main(int argc, char** argv) {
    // ── 可跳过解密/基准测试（L1 下极慢）──
    bool skip_decrypt = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "--no-decrypt") == 0)
            skip_decrypt = true;
    }

    // ── 确保输出目录存在（避免后续 append 时 WARN）──
    {
        std::error_code ec;
        std::filesystem::create_directories("bendmarking_output", ec);
    }

    // ── 计时起点 ──
    auto t_start = std::chrono::high_resolution_clock::now();

    // ── 横幅 ──
    std::cout << "╔════════════════════════════════════════════════╗\n"
              << "║  LatticeCryBendmarking — C/C++ Implementation  ║\n"
              << "║  @Author: Ziyi Dong, 2026                      ║\n"
              << "╚════════════════════════════════════════════════╝\n";

    // ── 参数 ──
    const auto p = unified::default_mp12_params();

    std::cout << "\nParameters:\n"
              << "  n = " << p.n << "  (lattice dimension)\n"
              << "  q = " << p.q << "  (modulus)\n"
              << "  b = " << p.b << "  (gadget base)\n"
              << "  k = " << p.k << "  (k = ceil(log_b q))\n"
              << "  m_bar = " << p.m_bar << "\n"
              << "  m = " << p.m << "  (total columns of A)\n"
              << "  sigma = " << p.sigma << "  (trapdoor Gaussian width)\n"
              << "  s = " << p.s << "  (preimage Gaussian width)\n";

    // ── 将参数写入基准测试输出文件 ──
    {
        constexpr const char* OUT_PATH = "bendmarking_output/bendmarking_plain-LWE.txt";
        std::ofstream fout(OUT_PATH, std::ios::out);
        if (fout.is_open()) {
            fout << "\n==========================================================\n"
                 << "  LatticeCryBenchmarking — Run Parameters\n"
                 << "==========================================================\n"
                 << "  n     = " << p.n   << "  (lattice dimension)\n"
                 << "  q     = " << p.q   << "  (modulus)\n"
                 << "  b     = " << p.b   << "  (gadget base)\n"
                 << "  k     = " << p.k   << "  (k = ceil(log_b q))\n"
                 << "  m_bar = " << p.m_bar << "\n"
                 << "  m     = " << p.m   << "  (total columns of A)\n"
                 << "  sigma = " << p.sigma << "  (trapdoor Gaussian width)\n"
                 << "  s     = " << p.s   << "  (preimage Gaussian width)\n";
            fout.close();
            std::cout << "  [Parameters written to " << OUT_PATH << "]\n";
        }
    }

    // ── 基础 MP12 陷门测试（Test 1–7）──
    run_mp12_trap_tests(p);

    // ── 委托陷门测试（Test 7–11）──
    run_del_tests(p);

    // ── 其他模块测试 / 基准 ──
    std::cout << "\n--- Running consolidated demos from other translation units ---\n";
    run_test_encrypt();
    run_test_powersof();
    run_test_powersof_modswitch();
    run_test_frd();
    run_bench_matops();
    run_test_expand();
    run_test_eval();

    if (!skip_decrypt) {
        run_test_decrypt();
        run_bench_decrypt();
    } else {
        std::cout << "\n--- Skipping decrypt & bench_decrypt (--no-decrypt) ---\n\n";
    }

    run_test_security();

    // ── 总耗时 ──
    auto t_end = std::chrono::high_resolution_clock::now();
    double total_s  = std::chrono::duration<double>(t_end - t_start).count();
    double total_ms = total_s * 1000.0;

    std::cout << "\nAll tests completed.\n"
              << "Total elapsed time: " << std::fixed << std::setprecision(2)
              << total_s << " s  (" << total_ms << " ms)\n";

    /* 写入文件 */
    {
        constexpr const char* OUT_PATH = "bendmarking_output/bendmarking_plain-LWE.txt";
        std::ofstream fout(OUT_PATH, std::ios::app);
        if (fout.is_open()) {
            fout << "\n----------------------------------------------------------\n"
                 << "  Total elapsed time: " << std::fixed << std::setprecision(2)
                 << total_s << " s  (" << total_ms << " ms)\n"
                 << "----------------------------------------------------------\n";
            fout.close();
            std::cout << "  [Total time written to " << OUT_PATH << "]\n";
        }
    }

    return 0;
}
