#ifndef INPUT_STATS_H
#define INPUT_STATS_H

#pragma once

#include <unordered_map>
#include <vector>
#include <string>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <plog/Log.h>

bool pair_cmp(std::pair<uint8_t, uint64_t>& a,
              std::pair<uint8_t, uint64_t>& b)
{
    return a.second < b.second;
}

struct InputStats
{
        //uint8_t num_bins{0};
        uint32_t num_files{0};
        std::unordered_map<uint8_t, uint64_t> records_per_bin{};
        std::unordered_map<uint8_t, uint64_t> hashes_per_bin{};

    public:
        InputStats() = default;
        InputStats(InputStats const &) = default;
        InputStats(InputStats &&) = default;
        InputStats & operator=(InputStats const &) = default;
        InputStats & operator=(InputStats &&) = default;
        ~InputStats() = default;

        std::vector<std::pair<uint8_t, uint64_t> > bins_by_size() const
        {
            std::vector<std::pair<uint8_t, uint64_t> > sorted_pairs;

            for (auto& it : hashes_per_bin) {
                sorted_pairs.push_back(it);
            }

            std::sort(sorted_pairs.begin(), sorted_pairs.end(), pair_cmp);
            return sorted_pairs;
        }

    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        try
        {
            //archive(num_bins);
            archive(num_files);
            archive(records_per_bin);
            archive(hashes_per_bin);
        }
            // GCOVR_EXCL_START
        catch (std::exception const & e)
        {
            PLOG_ERROR << "Cannot read input_stats: " + std::string{e.what()};
            exit(1);
        }

    }

    template <seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t & archive)
    {
        try
        {
            //archive(num_bins);
            archive(num_files);
            archive(records_per_bin);
            archive(hashes_per_bin);
        }
            // GCOVR_EXCL_START
        catch (std::exception const & e)
        {
            PLOG_ERROR << "Cannot read input_stats: " + std::string{e.what()};
            exit(1);
        }
    }
};

#endif // INPUT_STATS_H
