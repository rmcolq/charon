#ifndef INPUT_SUMMARY_H
#define INPUT_SUMMARY_H

#pragma once

#include <unordered_map>
#include <string>

#include <cereal/types/string.hpp>
#include <cereal/types/unordered_map.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/utility.hpp>
#include <seqan3/core/concept/cereal.hpp>
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

        uint8_t category_index(const std::string category) const
        {
            for (auto i=0; i<categories.size(); ++i)
            {
                if (category == categories.at(i))
                    return i;
            }
            return std::numeric_limits<uint8_t>::max();
        }

        uint8_t host_category_index() const
        {
            auto index1 = category_index("human");
            auto index2 = category_index("host");
            auto index = std::min(index1, index2);
            if (index == std::numeric_limits<uint8_t>::max())
                PLOG_ERROR << "Neither 'human' nor 'host' appear as categories in the index";
            assert(index != std::numeric_limits<uint8_t>::max());
            return index;
        }

        std::string category_name(const uint8_t index) const
        {
            if (index > categories.size())
                return "";
            else
                return categories.at(index);
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
