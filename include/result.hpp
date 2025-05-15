#ifndef CHARON_RESULT_H
#define CHARON_RESULT_H

#pragma once

#include <string>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <utility>
#include <plog/Log.h>

#include "read_entry.hpp"
#include "dehost_arguments.hpp"
#include "input_summary.hpp"
#include "classify_stats.hpp"

struct ResultSummary
{
    std::vector<uint64_t> classified_counts;
    uint64_t unclassified_count{0};
    uint64_t extracted_count{0};

    ResultSummary(const uint8_t size):
            classified_counts(size,0)
    {};
};

template<class record_type>
struct ReadRecord {
    bool is_paired;
    ReadEntry read;
    record_type record;
    record_type record2;
};

template<class record_type, class outfile_field_ids, class outfile_format>
class Result
{
    private:
        InputSummary input_summary_;
        ResultSummary result_summary_;
        StatsModel stats_model_;
        std::vector<ReadRecord<record_type>> cached_reads_;

        bool run_extract_;
        std::unordered_map<uint8_t, std::vector<seqan3::sequence_file_output<outfile_field_ids, outfile_format>>> extract_handles_;

    public:
        Result() = default;
        Result(Result const &) = default;
        Result(Result &&) = default;
        Result & operator=(Result const &) = default;
        Result & operator=(Result &&) = default;
        ~Result() = default;

        Result(const ClassifyArguments& opt, const InputSummary & summary):
                input_summary_{summary},
                result_summary_(summary.num_categories()),
                run_extract_(opt.run_extract)
        {
            stats_model_ = StatsModel(opt, summary);
            if (opt.run_extract){
                for (const auto [category_index,extract_files] : opt.extract_category_to_file) {
                    for (const auto extract_file : extract_files)
                        extract_handles_[category_index].push_back(seqan3::sequence_file_output{extract_file});
                    cached_reads_.reserve(opt.num_reads_to_fit*summary.num_categories()*4);
                }
            }


        };

        Result(const DehostArguments& opt, const InputSummary & summary):
                input_summary_{summary},
                result_summary_(summary.num_categories()),
                run_extract_(opt.run_extract)
        {
            stats_model_ = StatsModel(opt, summary);
            if (opt.run_extract){
                for (const auto [category_index,extract_files] : opt.extract_category_to_file) {
                    for (const auto extract_file : extract_files)
                        extract_handles_[category_index].push_back(seqan3::sequence_file_output{extract_file});
                    cached_reads_.reserve(opt.num_reads_to_fit*summary.num_categories()*4);
                }
            }
        };

        const InputSummary& input_summary() const
        {
            return input_summary_;
        }

        const uint8_t category_index(const std::string& category) const
        {
            return input_summary_.category_index(category);
        }

        uint8_t classify_read(ReadEntry& read_entry, const bool dehost=false)
        {
            PLOG_VERBOSE << "Classify read " << read_entry.read_id();
            if (dehost)
                read_entry.dehost(stats_model_, input_summary_.host_category_index());
            else
                read_entry.classify(stats_model_);
#pragma omp critical(print_result)
            read_entry.print_assignment_result(input_summary_);
#pragma omp critical(update_result_count)
            {
                const auto &call = read_entry.call();
                if (call < std::numeric_limits<uint8_t>::max())
                {
                    result_summary_.classified_counts[call] += 1;
                } else {
                    result_summary_.unclassified_count += 1;
                }
            }
            return read_entry.call();

        }

        void extract_read(const uint8_t category_index, const record_type& record)
        {
#pragma omp critical(extract_read)
            extract_handles_[category_index][0].push_back(record);
        }

        void extract_paired_read(const uint8_t category_index, const record_type& record, const record_type& record2)
        {
#pragma omp critical(extract_read)
            extract_handles_[category_index][0].push_back(record);
#pragma omp critical(extract_read2)
            extract_handles_[category_index][1].push_back(record2);
        }

        void add_read(ReadEntry& read_entry, const record_type& record, bool dehost=false){
            if (stats_model_.ready()) {
                auto category_index = classify_read(read_entry, dehost);
                if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()){
                    extract_read(category_index, record);
                }
            } else {
                PLOG_VERBOSE << "Add read " << read_entry.read_id() << " to training ";
#pragma omp critical(add_to_cache)
                {
                    bool training_complete = false;
                    if (cached_reads_.size() < cached_reads_.capacity())
                    {
                        cached_reads_.emplace_back(ReadRecord(false, read_entry, record, record));
                        training_complete = stats_model_.add_read_to_training_data(read_entry.unique_proportions());
                    } else {
                        stats_model_.force_ready();
                        training_complete = true;
                    }

                    if (training_complete)
                        classify_cache(dehost);
                }
            }
        }

        void add_paired_read(ReadEntry& read_entry, const record_type& record, const record_type& record2, const bool dehost=false){
            if (stats_model_.ready()) {
                auto category_index = classify_read(read_entry, dehost);
                if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()){
                    extract_paired_read(category_index, record, record2);
                }
            } else {
                PLOG_VERBOSE << "Add read " << read_entry.read_id() << " to training ";
#pragma omp critical(add_to_cache)
                {
                    bool training_complete = false;
                    if (cached_reads_.size() < cached_reads_.capacity())
                    {
                        cached_reads_.emplace_back(ReadRecord(true, read_entry, record, record2));
                        training_complete = stats_model_.add_read_to_training_data(read_entry.unique_proportions());
                    } else {
                        stats_model_.force_ready();
                        training_complete = true;
                    }

                    if (training_complete)
                        classify_cache(dehost);
                }
            }
        }

        void classify_cache(const bool dehost=false)
        {
            PLOG_VERBOSE << "Classify cached reads";

            for (auto read_record : cached_reads_)
            {
                auto & read_entry = read_record.read;
                const auto & record = read_record.record;
                auto category_index = classify_read(read_entry, dehost);
                if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()){
                    if (read_record.is_paired)
                    {
                        const auto & record2 = read_record.record2;
                        extract_paired_read(category_index, record, record2);
                    } else {
                        extract_read(category_index, record);
                    }
                }
            }
            cached_reads_.resize(0);
        }

        void complete(const bool dehost=false)
        {
            classify_cache(dehost);
        }



        void print_summary() const {

            PLOG_INFO << "Results summary: ";
            for (auto i = 0; i < result_summary_.classified_counts.size(); i++) {
                const auto & category = input_summary_.categories.at(i);
                PLOG_INFO << category << " :\t\t" << result_summary_.classified_counts.at(i);
            }
            PLOG_INFO << "unclassified :\t" << result_summary_.unclassified_count;
        }
    };


#endif // CHARON_RESULT_H
