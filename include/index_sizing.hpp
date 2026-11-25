#ifndef CHARON_INDEX_SIZING_H
#define CHARON_INDEX_SIZING_H

#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>

#include "index_arguments.hpp"
#include "input_stats.hpp"
#include "input_summary.hpp"
#include "layout.hpp"
#include "utils.hpp"

// Bin heights are sized for 1 % more minimisers than estimated: the HyperLogLog estimates are about
// 0.2 % accurate (measured 0.12-0.17 % low), and a hashed split of a group into parts adds a little
// variation per part, so this keeps every bin at or below the target fpr.
constexpr double sizing_margin = 1.01;

inline uint64_t with_sizing_margin(const uint64_t estimated_minimisers) {
    return static_cast<uint64_t>(std::ceil(static_cast<double>(estimated_minimisers) * sizing_margin));
}

// Bloom filter bits per element for a target fpr with num_hash hash functions
inline double bits_per_element(const double fpr, const uint8_t num_hash) {
    return -static_cast<double>(num_hash) / std::log(1 - std::exp(std::log(fpr) / num_hash));
}

// A category's target fpr: its --category-fpr when that is lower than max_fpr, otherwise max_fpr.
inline double category_target_fpr(const IndexArguments &opt, const std::string &category) {
    const auto it = opt.category_fpr.find(category);
    if (it == opt.category_fpr.end())
        return opt.max_fpr;
    return std::min(it->second, opt.max_fpr);
}

// Bits a minimiser of this category needs, relative to one at max_fpr (1 for a category at max_fpr)
inline double category_weight(const IndexArguments &opt, const std::string &category) {
    return bits_per_element(category_target_fpr(opt, category), opt.num_hash) /
           bits_per_element(opt.max_fpr, opt.num_hash);
}

// Height of the technical bins of a balanced index: enough for every group's parts to meet its own
// category's target fpr. Sizing only the fullest technical bin at max_fpr is not enough. When no bin is
// filled to capacity at max_fpr (a host-only index, whose bins are all loaded for the lower host fpr),
// the IBF would shrink until the host bins were at max_fpr. Expects the summary and stats returned by
// balance_layout, whose bins are the groups.
inline size_t balanced_bin_bits(const IndexArguments &opt, const InputSummary &summary, const InputStats &stats,
                                const Layout &layout) {
    size_t num_bits = 1;
    IndexArguments group_opt = opt;
    for (size_t group = 0; group < layout.num_groups(); ++group) {
        // Layout::part is a hash, so a group's minimisers spread evenly over its parts. Group loads add up their
        // bins' counts, so for packed groups whose bins share minimisers they are an upper bound.
        const uint64_t parts = layout.num_parts[group];
        const uint64_t load = (stats.hashes_per_bin.at(static_cast<uint8_t>(group)) + parts - 1) / parts;
        group_opt.max_fpr = category_target_fpr(opt, summary.bin_to_category.at(static_cast<uint8_t>(group)));
        num_bits = std::max(num_bits, bin_size_in_bits(group_opt, with_sizing_margin(load)));
    }
    return num_bits;
}

#endif // CHARON_INDEX_SIZING_H
