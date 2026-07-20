#pragma once

#include <unordered_map>
#include <vector>
#include <string>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <plog/Log.h>

struct InputStats {
    uint32_t num_files{0};
    std::unordered_map<uint8_t, uint64_t> records_per_bin{};
    std::unordered_map<uint8_t, uint64_t> hashes_per_bin{};

    InputStats() = default;
    InputStats(InputStats const&) = default;
    InputStats(InputStats&&) = default;
    InputStats& operator=(InputStats const&) = default;
    InputStats& operator=(InputStats&&) = default;
    ~InputStats() = default;

    std::vector<std::pair<uint8_t, uint64_t>> bins_by_size() const {
        std::vector<std::pair<uint8_t, uint64_t>> sorted_pairs;
        for (const auto& it : hashes_per_bin) {
            sorted_pairs.push_back(it);
        }
        std::sort(sorted_pairs.begin(), sorted_pairs.end(),
            [](const auto& left, const auto& right) {
                return left.second < right.second;
            });
        return sorted_pairs;
    }

    // Fixed: was non-const despite not mutating state.
    // Fixed: was undefined behaviour on empty hashes_per_bin — now returns 0.
    uint64_t max_num_hashes() const {
        if (hashes_per_bin.empty())
            return 0;
        return bins_by_size().back().second;
    }

    template<seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t& archive) {
        try {
            archive(num_files);
            archive(records_per_bin);
            archive(hashes_per_bin);
        } catch (std::exception const& e) {
            PLOG_ERROR << "Cannot read input_stats: " + std::string{e.what()};
            exit(1);
        }
    }

    template<seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t& archive) {
        try {
            archive(num_files);
            archive(records_per_bin);
            archive(hashes_per_bin);
        } catch (std::exception const& e) {
            PLOG_ERROR << "Cannot read input_stats: " + std::string{e.what()};
            exit(1);
        }
    }
};