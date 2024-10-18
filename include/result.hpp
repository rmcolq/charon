#ifndef SIFTER_RESULT_H
#define SIFTER_RESULT_H

#pragma once
#include <string>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <plog/Log.h>

#include <entry.hpp>

class Result
{
    private:
        uint8_t num_bins_{};
        std::unordered_map<std::string, ReadEntry> entries;

    public:
        Result() = default;
        Result(Result const &) = default;
        Result(Result &&) = default;
        Result & operator=(Result const &) = default;
        Result & operator=(Result &&) = default;
        ~Result() = default;

        Result(const uint8_t num_bins):
            num_bins_{num_bins}
        {};

        void update_entry(const std::string read_id, const auto & entry){
            if (entries.find(read_id) == entries.end()){
                entries[read_id] = ReadEntry(read_id, num_bins_);
            }
            entries[read_id].update_entry(entry);
        };
    };

#endif // SIFTER_RESULT_H
