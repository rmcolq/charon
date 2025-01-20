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
            cached_read_ids.reserve(opt.num_reads_to_fit*summary.num_categories()*2);
        };

        uint8_t category_index(const std::string& category)
        {
            return summary_.category_index(category);
        }

        uint8_t call(const std::string& read_id)
        {
            PLOG_VERBOSE << "Check call for read " << read_id;
            return entries.at(read_id).call();
        }

        void update_entry(const std::string read_id, const uint32_t length, const auto & entry){
            /*std::cout << read_id << " ";
            for (const auto i : entry){
                std:: cout << +i;
            }
            std::cout << std::endl;*/
            if (entries.find(read_id) == entries.end()){
                PLOG_VERBOSE << "Define entry for " << read_id << " with length " << length;
#pragma omp critical(add_read_id_to_entries)
                entries[read_id] = ReadEntry(read_id, length, summary_);
            }
            //PLOG_VERBOSE << "Update entry ";
            entries.at(read_id).update_entry(entry);
        };

        void classify_read(const std::string & read_id)
        {
            entries.at(read_id).classify(stats_model);
            entries.at(read_id).print_assignment_result(summary_);
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
            PLOG_VERBOSE << "Post-process read " << read_id;
            entries.at(read_id).post_process(summary_);
            if (stats_model.ready()) {
                PLOG_VERBOSE << "Classify read " << read_id;
                entries.at(read_id).classify(stats_model);
            } else {
                PLOG_VERBOSE << "Add read to training " << read_id;
#pragma omp critical
                {
                    cached_read_ids.push_back(read_id);
                    const auto training_complete = stats_model.add_read_to_training_data(entries.at(read_id).unique_props());
                    if (training_complete)
                        classify_cache();
                }
            }
        }

        void print_summary() {
            uint64_t unclassified = 0;
            std::vector<uint64_t> call_counts(summary_.num_categories(), 0);
            for (const auto & [read_id, entry] : entries){
                const auto & call = entry.call();
                if (call == std::numeric_limits<uint8_t>::max())
                    unclassified += 1;
                else
                    call_counts.at(call) += 1;
            }
            PLOG_INFO << "Results summary: ";
            for (auto i = 0; i < call_counts.size(); i++) {
                const auto & category = summary_.categories.at(i);
                PLOG_INFO << category << " :\t\t" << call_counts.at(i);
            }
            PLOG_INFO << "unclassified :\t" << unclassified;
        }
    };


#endif // CHARON_RESULT_H
