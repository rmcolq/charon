#ifndef CHARON_LAYOUT_H
#define CHARON_LAYOUT_H

#pragma once

#include <cstdint>
#include <vector>

#include <cereal/types/vector.hpp>

/**
 * @brief Maps scoring groups (the "bins" of InputSummary) onto IBF technical bins.
 *
 * A group whose minimisers exceed the bin capacity is split into num_parts technical bins;
 * minimiser x of group g is stored in (and only looked up in) technical bin
 * first_bin[g] + part(x, num_parts[g]). Because a minimiser only ever lives in one part, a
 * split group behaves like one tall bin with the same false positive rate, and all parts sit
 * in the same IBF row, so a lookup costs no extra memory accesses.
 * The identity layout (one technical bin per group) is the pre-balancing index format.
 */
struct Layout {
    std::vector<uint16_t> first_bin{};
    std::vector<uint16_t> num_parts{};
    uint32_t num_technical_bins{0};
    // Derived (not stored): for each technical bin, its group and which part of the group it is.
    std::vector<uint16_t> technical_group{};
    std::vector<uint16_t> technical_part{};

    // Fills technical_group/technical_part; call after first_bin/num_parts change or after loading.
    void finalise() {
        technical_group.assign(num_technical_bins, 0);
        technical_part.assign(num_technical_bins, 0);
        for (size_t group = 0; group < first_bin.size(); ++group) {
            for (uint16_t part = 0; part < num_parts[group]; ++part) {
                technical_group[first_bin[group] + part] = static_cast<uint16_t>(group);
                technical_part[first_bin[group] + part] = part;
            }
        }
    }

    [[nodiscard]] static Layout identity(const size_t num_groups) {
        Layout layout;
        for (size_t group = 0; group < num_groups; ++group) {
            layout.first_bin.push_back(static_cast<uint16_t>(group));
            layout.num_parts.push_back(1);
        }
        layout.num_technical_bins = static_cast<uint32_t>(num_groups);
        layout.finalise();
        return layout;
    }

    [[nodiscard]] size_t num_groups() const {
        return first_bin.size();
    }

    [[nodiscard]] bool is_identity() const {
        return num_technical_bins == num_groups();
    }

    // Which of `parts` partitions a minimiser hash belongs to. Hashes are mixed first because
    // seqan3 minimiser hashes are k-mer codes XOR a seed, whose low bits just encode the last bases.
    [[nodiscard]] static size_t part(const uint64_t hash, const size_t parts) {
        if (parts <= 1)
            return 0;
        const uint64_t mixed = hash * 0x9E3779B97F4A7C15ULL;
        return static_cast<size_t>(((mixed >> 32) * parts) >> 32);
    }

    [[nodiscard]] size_t technical_bin(const size_t group, const uint64_t hash) const {
        return first_bin[group] + part(hash, num_parts[group]);
    }

    template<class archive_t>
    void save(archive_t &archive) const {
        archive(first_bin, num_parts, num_technical_bins);
    }

    template<class archive_t>
    void load(archive_t &archive) {
        archive(first_bin, num_parts, num_technical_bins);
        finalise();
    }
};

#endif // CHARON_LAYOUT_H
