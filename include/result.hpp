#ifndef CHARON_RESULT_H
#define CHARON_RESULT_H

#pragma once

#include <string>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <utility>
#include <plog/Log.h>

#include "entry.hpp"
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
        uint8_t extract_category_;
        seqan3::sequence_file_output<outfile_field_ids, outfile_format> extract_handle_;
        seqan3::sequence_file_output<outfile_field_ids, outfile_format> extract_handle2_;

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
                run_extract_(opt.run_extract),
                extract_handle_{std::cout, seqan3::format_fasta{}},
                extract_handle2_{std::cout, seqan3::format_fasta{}}
        {
            stats_model_ = StatsModel(opt, summary);
            if (opt.run_extract){
                extract_handle_ = seqan3::sequence_file_output{opt.extract_file};
                if (opt.is_paired)
                    extract_handle2_ = seqan3::sequence_file_output{opt.extract_file2};

                extract_category_ = category_index(opt.category_to_extract);
            }
            cached_reads_.reserve(opt.num_reads_to_fit*summary.num_categories()*4);

        };

        const InputSummary& input_summary() const
        {
            return input_summary_;
        }

        const uint8_t category_index(const std::string& category) const
        {
            return input_summary_.category_index(category);
        }

        bool classify_read(ReadEntry& read_entry)
        {
            PLOG_VERBOSE << "Classify read " << read_entry.read_id();
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
            return run_extract_ and read_entry.call() == extract_category_;
        }

        void extract_read(const record_type& record)
        {
#pragma omp critical(extract_read)
            extract_handle_.push_back(record);
        }

        void extract_paired_read(const record_type& record, const record_type& record2)
        {
#pragma omp critical(extract_read)
                extract_handle_.push_back(record);
#pragma omp critical(extract_read2)
                extract_handle2_.push_back(record2);
        }

        void add_read(ReadEntry& read_entry, const record_type& record){
            if (stats_model_.ready()) {
                auto read_to_extract = classify_read(read_entry);
                if (run_extract_ and read_to_extract){
                    extract_read(record);
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
                        classify_cache();
                }
            }
        }

        void add_paired_read(ReadEntry& read_entry, const record_type& record, const record_type& record2){
            if (stats_model_.ready()) {
                auto read_to_extract = classify_read(read_entry);
                if (run_extract_ and read_to_extract){
                    extract_paired_read(record, record2);
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
                        classify_cache();
                }
            }
        }

        void classify_cache()
        {
            PLOG_VERBOSE << "Classify cached reads";

            for (auto read_record : cached_reads_)
            {
                auto & read_entry = read_record.read;
                const auto & record = read_record.record;
                bool read_to_extract = classify_read(read_entry);
                if (run_extract_ and read_to_extract){
                    if (read_record.is_paired)
                    {
                        const auto & record2 = read_record.record2;
                        extract_paired_read(record, record2);
                    } else {
                        extract_read(record);
                    }
                }
            }
            cached_reads_.resize(0);
        }

        void complete()
        {
            classify_cache();
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
