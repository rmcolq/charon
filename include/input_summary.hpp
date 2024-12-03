#ifndef INPUT_SUMMARY_H
#define INPUT_SUMMARY_H

#pragma once

#include <unordered_map>
#include <string>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <plog/Log.h>

struct InputSummary
{
        uint8_t num_bins{0};
        uint32_t num_files{0};
        std::unordered_map<uint8_t, uint64_t> records_per_bin{};
        std::unordered_map<uint8_t, uint64_t> hashes_per_bin{};
    public:
        InputSummary() = default;
        InputSummary(InputSummary const &) = default;
        InputSummary(InputSummary &&) = default;
        InputSummary & operator=(InputSummary const &) = default;
        InputSummary & operator=(InputSummary &&) = default;
        ~InputSummary() = default;

    template <seqan3::cereal_archive archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        try
        {
            archive(num_bins);
            archive(num_files);
            archive(records_per_bin);
            archive(hashes_per_bin);
        }
            // GCOVR_EXCL_START
        catch (std::exception const & e)
        {
            PLOG_ERROR << "Cannot read index: " + std::string{e.what()};
            exit(1);
        }

    }

    template <seqan3::cereal_input_archive archive_t>
    void load_parameters(archive_t & archive)
    {
        try
        {
            archive(num_bins);
            archive(num_files);
            archive(records_per_bin);
            archive(hashes_per_bin);
        }
            // GCOVR_EXCL_START
        catch (std::exception const & e)
        {
            PLOG_ERROR << "Cannot read index: " + std::string{e.what()};
            exit(1);
        }
    }
};

#endif // INPUT_SUMMARY_H
