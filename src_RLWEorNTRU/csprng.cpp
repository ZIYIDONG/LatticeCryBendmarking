/**
 * @file csprng.cpp
 * @brief CSPRNG 系统实现
 */

#include "../include_RLWEorNTRU/csprng.h"

#include <cstring>
#include <cstdio>
#include <stdexcept>
#include <limits>
#include <thread>
#include <chrono>

#ifdef __linux__
#include <sys/random.h>
#endif

namespace ibags {

// ============================================================================
// SystemCSPRNG
// ============================================================================

void SystemCSPRNG::randombytes(uint8_t* out, size_t len)
{
    if (len == 0) return;

#ifdef __linux__
    // 首选 getrandom()
    ssize_t ret = getrandom(out, len, 0);
    if (static_cast<size_t>(ret) == len) {
        return;
    }

    // getrandom 失败时回退到 /dev/urandom
    // (在 Linux 3.17+ 内核上 getrandom 总是可用，此回退用于兼容)
#endif

    // 回退: /dev/urandom
    FILE* f = fopen("/dev/urandom", "rb");
    if (!f) {
        throw std::runtime_error(
            "SystemCSPRNG: cannot open /dev/urandom");
    }

    size_t nread = fread(out, 1, len, f);
    fclose(f);

    if (nread != len) {
        throw std::runtime_error(
            "SystemCSPRNG: short read from /dev/urandom");
    }
}

uint64_t SystemCSPRNG::rand_u64()
{
    uint64_t value;
    randombytes(reinterpret_cast<uint8_t*>(&value), sizeof(value));
    return value;
}

uint64_t SystemCSPRNG::rand_uint64_bound(uint64_t bound)
{
    if (bound == 0) return 0;
    if (bound == 1) return 0;

    // 拒绝采样: 确保均匀分布
    uint64_t limit = std::numeric_limits<uint64_t>::max() -
                     (std::numeric_limits<uint64_t>::max() % bound);
    for (;;) {
        uint64_t r = rand_u64();
        if (r < limit) {
            return r % bound;
        }
    }
}

// ============================================================================
// SeedCSPRNG
// ============================================================================

SeedCSPRNG::SeedCSPRNG(const uint8_t* seed,
                       size_t         seed_len,
                       std::string_view domain_label)
    : seed_(seed, seed + seed_len)
    , counter_(0)
    , domain_label_(domain_label)
{
    // 确保种子至少 32 字节 (256 位)
    if (seed_len < CSPRNG_SEED_BYTES) {
        // 不足则补零
        seed_.resize(CSPRNG_SEED_BYTES, 0);
    }
}

void SeedCSPRNG::randombytes(uint8_t* out, size_t len)
{
    if (len == 0) return;

    size_t offset = 0;
    while (offset < len) {
        // 构建 XOF: domain || counter || seed
        Xof xof(domain_label_);
        xof.absorb_u64(counter_);
        xof.absorb(ByteSpan(seed_.data(), seed_.size()));
        xof.finalize();

        size_t chunk = len - offset;
        if (chunk > CSPRNG_MAX_OUTPUT) {
            chunk = CSPRNG_MAX_OUTPUT;
        }

        xof.squeeze(out + offset, chunk);

        offset += chunk;
        counter_++;

        // 防止计数器溢出
        if (counter_ == 0) {
            throw std::runtime_error(
                "SeedCSPRNG: counter overflow — seed exhausted");
        }
    }
}

uint64_t SeedCSPRNG::rand_u64()
{
    uint64_t value;
    randombytes(reinterpret_cast<uint8_t*>(&value), sizeof(value));
    return value;
}

uint64_t SeedCSPRNG::rand_uint64_bound(uint64_t bound)
{
    if (bound == 0) return 0;
    if (bound == 1) return 0;

    uint64_t limit = std::numeric_limits<uint64_t>::max() -
                     (std::numeric_limits<uint64_t>::max() % bound);
    for (;;) {
        uint64_t r = rand_u64();
        if (r < limit) {
            return r % bound;
        }
    }
}

} // namespace ibags