#include "bench_utils.h"
#include <chrono>
#include <algorithm>
#include <numeric>
#include <cmath>

BenchStats run_benchmark(const std::function<void()>& func,
                          int total_iters, int warmup) {
    using Clock = std::chrono::high_resolution_clock;

    std::vector<double> times;
    times.reserve(total_iters);
    for (int i = 0; i < total_iters; ++i) {
        auto t0 = Clock::now();
        func();
        auto t1 = Clock::now();
        times.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
    }

    std::vector<double> active(times.begin() + warmup, times.end());
    std::sort(active.begin(), active.end());

    double sum = std::accumulate(active.begin(), active.end(), 0.0);
    int n = (int)active.size();
    double avg = sum / n;

    double median;
    if (n % 2 == 0) {
        median = (active[n / 2 - 1] + active[n / 2]) / 2.0;
    } else {
        median = active[n / 2];
    }

    double var = 0;
    for (double t : active) { double d = t - avg; var += d * d; }
    var /= (n > 1) ? (n - 1) : 1;

    return { avg, median, std::sqrt(var),
             active.front(), active.back(), n, std::move(active) };
}
