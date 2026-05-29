/**
 * @file test_ibags_all.cpp
 * @brief IBAGS 统一测试入口 — 串联 core / poly / ntru_trapgen 三个套件
 *
 * 构建:
 *   cmake --build . --target test_ibags_all -j$(nproc)
 * 运行:
 *   ./test_ibags_all                 # 详细模式（全部 PASSED 行，彩色）
 *   ./test_ibags_all -q              # 静默模式（仅失败 + 汇总）
 *   ./test_ibags_all -q --no-color   # 静默 + 禁用 ANSI 颜色
 *   ./test_ibags_all --no-color      # 详细 + 禁用 ANSI 颜色
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "../include_RLWEorNTRU/test_colors.h"

// 全局标志 — 各测试文件通过 extern 引用
bool g_ibags_quiet   = false;
bool g_ibags_nocolor = false;

// 前向声明 — 各翻译单元导出的 run_* 函数
int run_ibags_core_tests();
int run_ibags_poly_tests();
int run_ntru_trapgen_tests();

int main(int argc, char** argv) {
    // 解析命令行标志
    for (int i = 1; i < argc; ++i) {
        if (std::strcmp(argv[i], "-q") == 0 ||
            std::strcmp(argv[i], "--quiet") == 0) {
            g_ibags_quiet = true;
        } else if (std::strcmp(argv[i], "--no-color") == 0) {
            g_ibags_nocolor = true;
        }
    }

    // 颜色包装宏 (仅在 --no-color 时禁用)
    #define C(suffix) (g_ibags_nocolor ? "" : COLOR_##suffix)

    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  %sIBAGS — Unified Test Suite%s                  ║\n",
           C(CYAN), COLOR_RESET);
    printf("║  Core, Poly & NTT, NTRU TrapGen              ║\n");
    printf("║  Mode: %s%-36s%s ║\n",
           C(BOLD),
           g_ibags_quiet ? "QUIET (failures only)" : "VERBOSE (all results)",
           COLOR_RESET);
    if (g_ibags_nocolor)
        printf("║  Color: OFF                                  ║\n");
    printf("╚══════════════════════════════════════════════╝\n\n");

    int total_passed = 0;
    int total_failed = 0;

    // ── Suite 1: Core (errors, params, secure_memory) ──
    printf("┌─ Suite 1/3: Core ────────────────────────────┐\n\n");
    int ret1 = run_ibags_core_tests();
    if (ret1 != 0) {
        printf("\n└─ Suite 1/3: Core — %sONE OR MORE FAILURES%s\n\n",
               C(RED), COLOR_RESET);
        ++total_failed;
    } else {
        printf("\n└─ Suite 1/3: Core — %sALL PASSED%s\n\n",
               C(GREEN), COLOR_RESET);
        ++total_passed;
    }

    // ── Suite 2: Poly & NTT ──
    printf("┌─ Suite 2/3: Poly & NTT ──────────────────────┐\n\n");
    int ret2 = run_ibags_poly_tests();
    if (ret2 != 0) {
        printf("\n└─ Suite 2/3: Poly & NTT — %sONE OR MORE FAILURES%s\n\n",
               C(RED), COLOR_RESET);
        ++total_failed;
    } else {
        printf("\n└─ Suite 2/3: Poly & NTT — %sALL PASSED%s\n\n",
               C(GREEN), COLOR_RESET);
        ++total_passed;
    }

    // ── Suite 3: NTRU TrapGen ──
    printf("┌─ Suite 3/3: NTRU TrapGen ────────────────────┐\n\n");
    int ret3 = run_ntru_trapgen_tests();
    if (ret3 != 0) {
        printf("\n└─ Suite 3/3: NTRU TrapGen — %sONE OR MORE FAILURES%s\n\n",
               C(RED), COLOR_RESET);
        ++total_failed;
    } else {
        printf("\n└─ Suite 3/3: NTRU TrapGen — %sALL PASSED%s\n\n",
               C(GREEN), COLOR_RESET);
        ++total_passed;
    }

    // ── Grand total ──
    printf("╔══════════════════════════════════════════════╗\n");
    printf("║  %sIBAGS Unified Suite — Grand Total%s           ║\n",
           C(CYAN), COLOR_RESET);
    printf("║  Suites passed : %s%d / %d%s                      ║\n",
           (total_passed == total_passed + total_failed) ? C(GREEN) : C(YELLOW),
           total_passed, total_passed + total_failed,
           COLOR_RESET);
    printf("║  Suites failed : %s%d%s                          ║\n",
           total_failed ? C(RED) : C(GREEN),
           total_failed,
           COLOR_RESET);
    printf("╚══════════════════════════════════════════════╝\n");

    #undef C
    return (total_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
