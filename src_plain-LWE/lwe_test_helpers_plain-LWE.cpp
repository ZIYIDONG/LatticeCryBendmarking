#include "lwe_test_helpers_plain-LWE.h"
#include <cassert>

using namespace matops;

Mat generate_lwe_matrix(const Vec& t, int R, int M, long q,
                         int noise_bound, std::mt19937_64& rng)
{
    std::uniform_int_distribution<long> unif(0, q - 1);
    std::uniform_int_distribution<long> noise(-noise_bound, noise_bound);

    Mat A = matops::make_mat(R, M);

    for (int i = 0; i < R - 1; ++i)
        for (int j = 0; j < M; ++j)
            A[i][j] = unif(rng);

    assert(t[R - 1] == 1 && "需要 t[R-1] = 1 以简化 LWE 构造");

    for (int j = 0; j < M; ++j) {
        long inner = 0;
        for (int i = 0; i < R - 1; ++i)
            inner = matops::mod_pos(inner + t[i] * A[i][j], q);
        long e_j = noise(rng);
        A[R - 1][j] = matops::mod_pos(e_j - inner, q);
    }

    return A;
}

Mat gsw_encrypt(const Mat& A, const Mat& G, int mu,
                 long q, std::mt19937_64& rng)
{
    int R = (int)A.size();
    int M = (int)A[0].size();

    std::uniform_int_distribution<int> bit(0, 1);
    Mat S = matops::make_mat(M, M);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < M; ++j)
            S[i][j] = bit(rng);

    Mat AS = matops::mat_mul(A, S, q);

    Mat muG = matops::make_mat(R, M, 0);
    if (mu != 0) {
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < M; ++j)
                muG[i][j] = matops::mod_pos((long)mu * G[i][j], q);
    }

    return matops::mat_add(AS, muG, q);
}

Mat simple_expand(const Mat& C, int N_id) {
    int R = (int)C.size();
    int M = (int)C[0].size();
    Mat Ch(N_id * R, Vec(N_id * M, 0));

    for (int a = 0; a < N_id; ++a)
        for (int r = 0; r < R; ++r)
            for (int c = 0; c < M; ++c)
                Ch[a * R + r][a * M + c] = C[r][c];

    return Ch;
}

std::vector<Vec> generate_shared_keys(const Vec& t_master,
                                       int N_id, long q,
                                       std::mt19937_64& rng)
{
    int R = (int)t_master.size();
    std::uniform_int_distribution<long> unif(0, q - 1);
    std::vector<Vec> keys(N_id, Vec(R, 0));

    for (int i = 0; i < N_id - 1; ++i)
        for (int j = 0; j < R; ++j)
            keys[i][j] = unif(rng);

    for (int j = 0; j < R; ++j) {
        long partial_sum = 0;
        for (int i = 0; i < N_id - 1; ++i)
            partial_sum = matops::mod_pos(partial_sum + keys[i][j], q);
        keys[N_id - 1][j] = matops::mod_pos(t_master[j] - partial_sum, q);
    }

    return keys;
}
