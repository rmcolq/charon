#ifndef INPUT_SUMMARY_H
#define INPUT_SUMMARY_H

#pragma once

#include <unordered_map>
#include <string>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <plog/Log.h>

struct InputSummary
{
    uint8_t num_bins{0};
    std::vector<std::string> categories;
    std::vector<std::pair<std::string, uint8_t>> filepath_to_bin;
    std::unordered_map<uint8_t, std::string> bin_to_category;

    public:
        InputSummary() = default;
        InputSummary(InputSummary const &) = default;
        InputSummary(InputSummary &&) = default;
        InputSummary & operator=(InputSummary const &) = default;
        InputSummary & operator=(InputSummary &&) = default;
        ~InputSummary() = default;

        uint8_t num_categories() const
        {
            return categories.size();
        }

    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        try
        {
            archive(num_bins);
            archive(categories);
            archive(filepath_to_bin);
            archive(bin_to_category);
        }
            // GCOVR_EXCL_START
        catch (std::exception const & e)
        {
            PLOG_ERROR << "Cannot read input_summary: " + std::string{e.what()};
            exit(1);
        }

    }

    template <seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t & archive)
    {
        try
        {
            archive(num_bins);
            archive(categories);
            archive(filepath_to_bin);
            archive(bin_to_category);
        }
            // GCOVR_EXCL_START
        catch (std::exception const & e)
        {
            PLOG_ERROR << "Cannot read input_summary: " + std::string{e.what()};
            exit(1);
        }
    }
};

#endif // INPUT_SUMMARY_H
