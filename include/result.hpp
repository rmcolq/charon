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

        Result(const ClassifyArguments& opt, const InputSummary & summary):
            summary_{summary}
        {
            stats_model = StatsModel(opt, summary);
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
                PLOG_DEBUG << "Done";
            }
            PLOG_DEBUG << "Update entry ";
            entries.at(read_id).update_entry(entry);
        };

        void classify_read(const std::string & read_id)
        {
            entries.at(read_id).classify(stats_model);
            entries.at(read_id).print_result(summary_);
        }

        void classify_cache()
        {
            for (const auto & read_id : cached_read_ids)
            {
                classify_read(read_id);
            }
            cached_read_ids.clear();
        }

        void complete()
        {
            classify_cache();
        }

        void post_process_read(const std::string read_id) {
            PLOG_DEBUG << "Post-process read " << read_id;
            entries.at(read_id).post_process(summary_);
            if (stats_model.ready()) {
                PLOG_DEBUG << "Classify read " << read_id;
                entries.at(read_id).classify(stats_model);
            } else {
                PLOG_DEBUG << "Add read to training " << read_id;
                cached_read_ids.push_back(read_id);
                const auto training_complete = stats_model.add_read_to_training_data(entries.at(read_id).unique_props());
                if (training_complete)
                    classify_cache();
            }
        };
    };

#endif // CHARON_RESULT_H
