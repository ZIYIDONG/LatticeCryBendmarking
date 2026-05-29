#include "bench_utils_plain-LWE.h"
#include <fstream>
#include <iostream>

void bench_write(const std::string& content) {
    constexpr const char* OUT_PATH = "benchmarking_output/benchmarking_plain-LWE.txt";
    std::ofstream fout(OUT_PATH, std::ios::app);
    if (fout.is_open()) {
        fout << content;
        fout.close();
    } else {
        std::cerr << "  [WARN] Could not open " << OUT_PATH << " for writing\n";
    }
}
