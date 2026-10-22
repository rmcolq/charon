#ifndef CHARON_HYPERLOGLOG_H
#define CHARON_HYPERLOGLOG_H

#pragma once

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

/**
 * @brief HyperLogLog sketch estimating the number of distinct 64-bit values (minimiser hashes) in a stream.
 *
 * Used to size IBF bins without holding every distinct minimiser in memory. With the default precision of 18
 * (2^18 one-byte registers, 256 KB) the standard error is 1.04 / 2^9, about 0.2 %. Cardinalities are estimated
 * with Ertl's improved estimator ("New cardinality estimation algorithms for HyperLogLog sketches", 2017), which
 * needs no empirical bias correction and is accurate from a handful of values up to billions.
 */
class HyperLogLog {
    uint8_t precision_;
    std::vector<uint8_t> registers_;

    // splitmix64 finaliser: minimiser hashes are k-mer codes XOR a seed, so mix them before use
    [[nodiscard]] static constexpr uint64_t mix(uint64_t x) {
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ULL;
        x ^= x >> 27;
        x *= 0x94d049bb133111ebULL;
        x ^= x >> 31;
        return x;
    }

    [[nodiscard]] static double sigma(double x) {
        if (x == 1.0)
            return std::numeric_limits<double>::infinity();
        double y = 1.0;
        double z = x;
        double z_old;
        do {
            x *= x;
            z_old = z;
            z += x * y;
            y += y;
        } while (z != z_old);
        return z;
    }

    [[nodiscard]] static double tau(double x) {
        if (x == 0.0 or x == 1.0)
            return 0.0;
        double y = 1.0;
        double z = 1.0 - x;
        double z_old;
        do {
            x = std::sqrt(x);
            z_old = z;
            y *= 0.5;
            z -= (1.0 - x) * (1.0 - x) * y;
        } while (z != z_old);
        return z / 3.0;
    }

public:
    explicit HyperLogLog(const uint8_t precision = 18) : precision_{precision}, registers_(size_t{1} << precision, 0) {}

    void add(const uint64_t value) {
        const uint64_t hash = mix(value);
        const size_t index = hash >> (64 - precision_);
        const uint64_t rest = hash << precision_;
        const auto rank = static_cast<uint8_t>(rest == 0 ? 64 - precision_ + 1 : std::countl_zero(rest) + 1);
        if (rank > registers_[index])
            registers_[index] = rank;
    }

    [[nodiscard]] uint64_t estimate() const {
        const size_t q = 64 - precision_;
        std::vector<uint64_t> counts(q + 2, 0);
        for (const auto value: registers_)
            ++counts[value];
        const auto m = static_cast<double>(registers_.size());
        double z = m * tau(1.0 - static_cast<double>(counts[q + 1]) / m);
        for (size_t k = q; k >= 1; --k)
            z = 0.5 * (z + static_cast<double>(counts[k]));
        z += m * sigma(static_cast<double>(counts[0]) / m);
        return static_cast<uint64_t>(std::llround(m * m / (2.0 * std::log(2.0) * z)));
    }
};

#endif // CHARON_HYPERLOGLOG_H
