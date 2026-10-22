#ifndef CHARON_CONCURRENT_IBF_H
#define CHARON_CONCURRENT_IBF_H

#pragma once

#include <array>
#include <bit>
#include <cstdint>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

/**
 * @brief Inserts into an uncompressed seqan3 IBF from several threads at once.
 *
 * seqan3's emplace sets bits with a plain read-modify-write of a 64-bit word that holds 64 bins, so two threads
 * filling different bins can lose each other's bits. This sets exactly the same bits as
 * interleaved_bloom_filter::emplace, but with an atomic OR. It mirrors seqan3's private hash_and_fit and seeds;
 * the test that compares it with emplace guards against a seqan3 update changing them.
 * The IBF must outlive the inserter and must not be resized while it is in use.
 */
class ConcurrentIbfInserter {
    static constexpr std::array<uint64_t, 5> hash_seeds{13'572'355'802'537'770'549ULL,
                                                        13'043'817'825'332'782'213ULL,
                                                        10'650'232'656'628'343'401ULL,
                                                        16'499'269'484'942'379'435ULL,
                                                        4'893'150'838'803'335'377ULL};

    uint64_t *words_;
    uint64_t bin_size_;
    uint64_t hash_shift_;
    uint64_t technical_bins_;
    size_t hash_funs_;

    [[nodiscard]] uint64_t hash_and_fit(uint64_t h, const uint64_t seed) const {
        h *= seed;
        h ^= h >> hash_shift_;
        h *= 11'400'714'819'323'198'485ULL;
        h = static_cast<uint64_t>((static_cast<__uint128_t>(h) * static_cast<__uint128_t>(bin_size_)) >> 64);
        return h * technical_bins_;
    }

public:
    explicit ConcurrentIbfInserter(seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> &ibf) :
            words_{ibf.raw_data().data()},
            bin_size_{ibf.bin_size()},
            hash_shift_{static_cast<uint64_t>(std::countl_zero(static_cast<uint64_t>(ibf.bin_size())))},
            technical_bins_{((ibf.bin_count() + 63) >> 6) << 6},
            hash_funs_{ibf.hash_function_count()} {}

    void emplace(const uint64_t value, const size_t bin) const noexcept {
        for (size_t i = 0; i < hash_funs_; ++i) {
            const uint64_t idx = hash_and_fit(value, hash_seeds[i]) + bin;
            uint64_t *word = words_ + (idx >> 6);
            const uint64_t bit = uint64_t{1} << (idx & 63);
            // Skip the atomic write (and taking the cache line exclusively) when the bit is already set
            if ((__atomic_load_n(word, __ATOMIC_RELAXED) & bit) == 0)
                __atomic_fetch_or(word, bit, __ATOMIC_RELAXED);
        }
    }
};

#endif // CHARON_CONCURRENT_IBF_H
