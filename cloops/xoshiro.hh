#ifndef XOSHIRO_H_
#define XOSHIRO_H_

#include <cstdint>
#include <limits>

// xoshiro256++ 1.0 (Blackman & Vigna, 2018) - public domain reference
// algorithm: https://prng.di.unimi.it/xoshiro256plusplus.c
//
// A modern, fast, well-tested 64-bit PRNG. Roughly 2-4x faster per call than
// std::mt19937.
//
// Satisfies the standard's UniformRandomBitGenerator concept (operator(),
// min(), max(), result_type), so it works as a drop-in for std::mt19937 in
// std::uniform_real_distribution, std::exponential_distribution,
// std::discrete_distribution, etc.

struct xoshiro256pp {
    using result_type = std::uint64_t;

    std::uint64_t s[4];

    explicit xoshiro256pp(std::uint64_t seed) {
        auto splitmix64 = [state = seed]() mutable -> std::uint64_t {
            std::uint64_t z = (state += 0x9E3779B97F4A7C15ULL);
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            return z ^ (z >> 31);
        };
        for (auto& word : s) { word = splitmix64(); }
    }

    static constexpr std::uint64_t rotl(std::uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    std::uint64_t operator()() {
        const std::uint64_t result = rotl(s[0] + s[3], 23) + s[0];
        const std::uint64_t t = s[1] << 17;

        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = rotl(s[3], 45);

        return result;
    }

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return std::numeric_limits<std::uint64_t>::max(); }
};

#endif // XOSHIRO_H_
