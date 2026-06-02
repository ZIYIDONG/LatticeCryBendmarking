/**
 * bench_utils.cpp — 向后兼容 stub
 *
 * 实际 benchmark 逻辑在 bench_utils.h 的模板函数 run_benchmark_t 中。
 * 此文件仅为链接器提供 run_benchmark(std::function) 的实例化目标。
 */

#include "bench_utils.h"
// run_benchmark 的 std::function 版本在 bench_utils.h 中 inline 定义
