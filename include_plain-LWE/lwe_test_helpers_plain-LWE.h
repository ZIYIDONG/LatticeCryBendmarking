#pragma once
#include "matops_plain-LWE.h"
#include <random>
#include <vector>

matops::Mat generate_lwe_matrix(const matops::Vec& t, int R, int M, long q,
                                 int noise_bound, std::mt19937_64& rng);

matops::Mat gsw_encrypt(const matops::Mat& A, const matops::Mat& G, int mu,
                         long q, std::mt19937_64& rng);

matops::Mat simple_expand(const matops::Mat& C, int N_id);

std::vector<matops::Vec> generate_shared_keys(const matops::Vec& t_master,
                                               int N_id, long q,
                                               std::mt19937_64& rng);
