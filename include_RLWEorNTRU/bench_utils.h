#pragma once
#include <vector>
#include <functional>

struct BenchStats {
    double avg_us;
    double median_us;
    double stddev_us;
    double min_us;
    double max_us;
    int    active_samples;
    std::vector<double> all_us;
};

BenchStats run_benchmark(const std::function<void()>& func,
                          int total_iters = 1000,
                          int warmup = 50);
