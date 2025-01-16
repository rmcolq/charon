#ifndef CHARON_RESULT_H
#define CHARON_RESULT_H

#pragma once

#include <string>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <plog/Log.h>

#include "entry.hpp"
#include "input_summary.hpp"
#include "classify_stats.hpp"


class Result
{
    private:
        InputSummary summary_;
        std::unordered_map<std::string, ReadEntry> entries;
        StatsModel stats_model;
        std::vector<std::string> cached_read_ids;

    public:
        Result() = default;
        Result(Result const &) = default;
        Result(Result &&) = default;
        Result & operator=(Result const &) = default;
        Result & operator=(Result &&) = default;
        ~Result() = default;

        Result(const InputSummary summary):
            summary_{summary}
        {
            cached_read_ids.reserve(2000);
        };

        void update_entry(const std::string read_id, const uint16_t length, const auto & entry){
            /*std::cout << read_id << " ";
            for (const auto i : entry){
                std:: cout << +i;
            }
            std::cout << std::endl;*/
            if (entries.find(read_id) == entries.end()){
                PLOG_DEBUG << "Define entry for " << read_id << " with length " << length;
                entries[read_id] = ReadEntry(read_id, length, summary_);
            }
            //PLOG_DEBUG << "Update entry " << entry;
            entries[read_id].update_entry(entry);
        };

        /*void classify_read(const std::string read_id)
        {

        }*/

        void post_process_read(const std::string read_id) {
            if (entries.find(read_id) != entries.end()) {
                entries[read_id].post_process(summary_);
            }
        };
        void print_result(const std::string read_id) {
            if (entries.find(read_id) != entries.end()) {
                entries[read_id].print_result(summary_);
            }
        };
    };

#endif // CHARON_RESULT_H
