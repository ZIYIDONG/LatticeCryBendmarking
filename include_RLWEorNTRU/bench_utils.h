#pragma once
/**
 * @file bench_utils.h
 * @brief 计时工具 — 支持 μs 延迟、P99/P999 尾延迟、CPU 周期计数 (x86)
 *
 * 设计原则:
 *  - 模板函数消除 std::function 虚调用开销，适合 μs 级微基准
 *  - P99/P999 尾延迟用于评估实时/恒定时间特性
 *  - __rdtsc 周期计数避免 CPU 频率缩放对 μs 测量的干扰
 *
 * 用法:
 *   auto stats = run_benchmark([&]{ do_work(); }, 1000, 50);
 *   printf("avg=%.1f us, p99=%.1f us\n", stats.avg_us, stats.p99_us);
 */

#include <vector>
#include <functional>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <algorithm>

/// 平台相关的 CPU 时间戳计数器
/// x86-64: __rdtsc; 其他平台返回 0，仅依赖 chrono
inline uint64_t read_cpu_cycles() {
#if defined(__x86_64__) || defined(_M_X64)
    unsigned int lo, hi;
    __asm__ volatile("rdtsc" : "=a"(lo), "=d"(hi));
    return (static_cast<uint64_t>(hi) << 32) | lo;
#else
    return 0;
#endif
}

struct BenchStats {
    // ── μs 级指标 ──
    double avg_us;
    double median_us;
    double stddev_us;
    double min_us;
    double max_us;

    // ── 尾延迟 (百分位) ──
    double p99_us;    ///< 第 99 百分位延迟
    double p999_us;   ///< 第 99.9 百分位延迟

    // ── CPU 周期 (x86-64 only, 其他平台为 0) ──
    double avg_cycles; ///< 每次调用平均 CPU 周期数

    // ── 元信息 ──
    int active_samples;  ///< 有效样本数 (total_iters - warmup)
};

// ============================================================================
// 模板版本 — 零抽象开销，适合 μs 级微基准
// ============================================================================

/**
 * @brief 对可调用对象 func 运行 benchmark，返回统计结果。
 *
 * @tparam F  可调用类型 (lambda / 函数指针)，无 std::function 包装开销
 * @param func          被测函数
 * @param total_iters   总迭代次数 (含 warmup)
 * @param warmup        丢弃的前 N 次样本数
 * @return BenchStats   统计结果
 */
template <typename F>
BenchStats run_benchmark_t(F&& func, int total_iters, int warmup) {
    using Clock = std::chrono::high_resolution_clock;

    std::vector<double> times_us;
    std::vector<double> cycles;
    times_us.reserve(total_iters);

    const bool has_tsc = (read_cpu_cycles() != 0);
    if (has_tsc) {
        cycles.reserve(total_iters);
    }

    for (int i = 0; i < total_iters; ++i) {
        uint64_t c0 = has_tsc ? read_cpu_cycles() : 0;
        auto t0 = Clock::now();

        func();

        auto t1 = Clock::now();
        uint64_t c1 = has_tsc ? read_cpu_cycles() : 0;

        times_us.push_back(
            std::chrono::duration<double, std::micro>(t1 - t0).count());
        if (has_tsc) {
            cycles.push_back(static_cast<double>(c1 - c0));
        }
    }

    // 丢弃 warmup 样本
    std::vector<double> active(times_us.begin() + warmup, times_us.end());
    std::sort(active.begin(), active.end());

    const int n = static_cast<int>(active.size());

    // ── 基础统计量 ──
    double sum = 0.0;
    for (double t : active) sum += t;
    const double avg = sum / n;

    double median;
    if (n % 2 == 0) {
        median = (active[n / 2 - 1] + active[n / 2]) / 2.0;
    } else {
        median = active[n / 2];
    }

    double var = 0.0;
    for (double t : active) { double d = t - avg; var += d * d; }
    var /= (n > 1) ? (n - 1) : 1;

    // ── 百分位 (P99, P999) ──
    auto percentile = [&](double pct) -> double {
        double idx = pct / 100.0 * static_cast<double>(n - 1);
        int lo = static_cast<int>(idx);
        int hi = std::min(lo + 1, n - 1);
        double frac = idx - lo;
        return active[lo] * (1.0 - frac) + active[hi] * frac;
    };

    // ── 周期数均值 ──
    double avg_cycles = 0.0;
    if (has_tsc) {
        double sum_cyc = 0.0;
        for (int i = warmup; i < total_iters; ++i) sum_cyc += cycles[i];
        avg_cycles = sum_cyc / n;
    }

    return {avg, median, std::sqrt(var),
            active.front(), active.back(),
            percentile(99.0), percentile(99.9),
            avg_cycles, n};
}

// ============================================================================
// 便捷包装 (保持向后兼容 + std::function 接口)
// ============================================================================

/**
 * @brief 便捷版本 — 通过 std::function 调用，适合非热点路径。
 *        热点基准请使用 run_benchmark_t 模板版本。
 */
inline BenchStats run_benchmark(const std::function<void()>& func,
                                int total_iters = 1000,
                                int warmup = 50) {
    return run_benchmark_t(func, total_iters, warmup);
}
