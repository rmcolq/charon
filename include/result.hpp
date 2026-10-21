#ifndef CHARON_RESULT_H
#define CHARON_RESULT_H

#pragma once

#include <atomic>
#include <string>
#include <ankerl/unordered_dense.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <utility>
#include <plog/Log.h>

#include "memory_manager.hpp"
#include "read_entry.hpp"
#include "dehost_arguments.hpp"
#include "input_summary.hpp"
#include "classify_stats.hpp"

struct ResultSummary {
    std::vector<uint64_t> classified_counts;
    uint64_t unclassified_count{0};
    uint64_t extracted_count{0};

    ResultSummary(const uint8_t size) :
            classified_counts(size, 0) {};
};

template<class record_type>
struct ReadRecord {
    bool is_paired;
    ReadEntry read;
    record_type record;
    record_type record2;
};

template<class record_type, class outfile_field_ids, class outfile_format>
class Result {
private:
    InputSummary input_summary_;
    ResultSummary result_summary_;
    StatsModel stats_model_;
    std::vector<ReadRecord<record_type>> cached_reads_;
    size_t cache_bytes_{0};        // estimated bytes held in cached_reads_
    size_t cache_byte_budget_{0};  // 0 = unlimited
    // Cached reads waiting to be classified now the model is trained. They are moved out of
    // cached_reads_ inside the add_to_cache lock and classified in parallel by classify_pending,
    // outside any lock or parallel region.
    std::vector<ReadRecord<record_type>> pending_reads_;
    int num_threads_{1};

    bool run_extract_;
    // Set (release) once stats_model_ is trained; after that reads can be classified without the
    // add_read_to_results lock, because the model is only read from then on.
    std::atomic<bool> model_ready_{false};
    bool category_or_{false};
    ankerl::unordered_dense::map<uint8_t, std::vector<seqan3::sequence_file_output<outfile_field_ids, outfile_format>>> extract_handles_;

public:
    Result() = default;

    Result(Result const &) = default;

    Result(Result &&) = default;

    Result &operator=(Result const &) = default;

    Result &operator=(Result &&) = default;

    ~Result() = default;

    Result(const ClassifyArguments &opt, const InputSummary &summary) :
            input_summary_{summary},
            result_summary_(summary.num_categories()),
            cache_byte_budget_(opt.cache_byte_budget),
            num_threads_(opt.threads),
            run_extract_(opt.run_extract) {
        result_summary_.classified_counts.reserve(summary.num_categories());  // Pre-allocate
        stats_model_ = StatsModel(opt, summary);
        if (opt.run_extract) {
            for (const auto [category_index, extract_files]: opt.extract_category_to_file) {
                for (const auto extract_file: extract_files)
                    extract_handles_[category_index].push_back(seqan3::sequence_file_output{extract_file});
            }
        }
        // Reads are cached (and classified later) until the model is trained, whether or not we extract.
        cached_reads_.reserve(opt.num_reads_to_fit * summary.num_categories() * 4);
    };

    Result(const DehostArguments &opt, const InputSummary &summary) :
            input_summary_{summary},
            result_summary_(summary.num_categories()),
            cache_byte_budget_(opt.cache_byte_budget),
            num_threads_(opt.threads),
            run_extract_(opt.run_extract),
            category_or_(opt.scoring == "category_or") {
        result_summary_.classified_counts.reserve(summary.num_categories());  // Pre-allocate
        stats_model_ = StatsModel(opt, summary);
        if (opt.run_extract) {
            for (const auto [category_index, extract_files]: opt.extract_category_to_file) {
                for (const auto extract_file: extract_files)
                    extract_handles_[category_index].push_back(seqan3::sequence_file_output{extract_file});
            }
        }
        // Reads are cached (and classified later) until the model is trained, whether or not we extract.
        cached_reads_.reserve(opt.num_reads_to_fit * summary.num_categories() * 4);
    };

    [[nodiscard]] bool category_or() const {
        return category_or_;
    }

    [[nodiscard]] bool model_ready() const {
        return model_ready_.load(std::memory_order_acquire);
    }

    const InputSummary &input_summary() const {
        return input_summary_;
    }

    const uint8_t category_index(const std::string &category) const {
        return input_summary_.category_index(category);
    }

    uint8_t classify_read(ReadEntry &read_entry, const bool dehost = false) {
        PLOG_VERBOSE << "Classify read " << read_entry.read_id();
        if (dehost)
            read_entry.dehost(stats_model_, input_summary_.host_category_index());
        else
            read_entry.classify(stats_model_);
        const auto line = read_entry.assignment_line(input_summary_);
#pragma omp critical(print_result)
        std::cout << line;
#pragma omp critical(update_result_count)
        {
            const auto &call = read_entry.call();
            if (call < std::numeric_limits<uint8_t>::max()) {
                result_summary_.classified_counts[call] += 1;
            } else {
                result_summary_.unclassified_count += 1;
            }
        }
        return read_entry.call();

    }

    void extract_read(const uint8_t category_index, const record_type &record) {
#pragma omp critical(extract_read)
        extract_handles_[category_index][0].push_back(record);
    }

    void extract_paired_read(const uint8_t category_index, const record_type &record, const record_type &record2) {
#pragma omp critical(extract_read)
        extract_handles_[category_index][0].push_back(record);
#pragma omp critical(extract_read2)
        extract_handles_[category_index][1].push_back(record2);
    }

    void add_read(ReadEntry &read_entry, const record_type &record, bool dehost = false) {
        if (stats_model_.ready()) {
            auto category_index = classify_read(read_entry, dehost);
            if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()) {
                extract_read(category_index, record);
            }
        } else {
            PLOG_VERBOSE << "Add read " << read_entry.read_id() << " to training ";
#pragma omp critical(add_to_cache)
            {
                bool training_complete = false;
                bool cached = false;
                if (cache_has_room()) {
                    cache_bytes_ += MemoryManager::cached_read_bytes(read_entry.length(), read_entry.read_id().size(),
                                                                     input_summary_.num_categories(), read_entry.num_hashes());
                    cached_reads_.emplace_back(ReadRecord(false, read_entry, record, record_type{}));
                    cached = true;
                    training_complete = stats_model_.add_read_to_training_data(read_entry.unique_proportions());
                } else {
                    stats_model_.force_ready();
                    training_complete = true;
                }

                if (training_complete)
                    take_cache();
                if (not cached) {
                    // the cache was full, so this read was neither cached nor classified above
                    auto category_index = classify_read(read_entry, dehost);
                    if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()) {
                        extract_read(category_index, record);
                    }
                }
                if (stats_model_.ready())
                    model_ready_.store(true, std::memory_order_release);
            }
        }
    }

    void add_paired_read(ReadEntry &read_entry, const record_type &record, const record_type &record2,
                         const bool dehost = false) {
        if (stats_model_.ready()) {
            auto category_index = classify_read(read_entry, dehost);
            if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()) {
                extract_paired_read(category_index, record, record2);
            }
        } else {
            PLOG_VERBOSE << "Add read " << read_entry.read_id() << " to training ";
#pragma omp critical(add_to_cache)
            {
                bool training_complete = false;
                bool cached = false;
                if (cache_has_room()) {
                    cache_bytes_ += MemoryManager::cached_read_bytes(read_entry.length(), read_entry.read_id().size(),
                                                                     input_summary_.num_categories(), read_entry.num_hashes(), true);
                    cached_reads_.emplace_back(ReadRecord(true, read_entry, record, record2));
                    cached = true;
                    training_complete = stats_model_.add_read_to_training_data(read_entry.unique_proportions());
                } else {
                    stats_model_.force_ready();
                    training_complete = true;
                }

                if (training_complete)
                    take_cache();
                if (not cached) {
                    // the cache was full, so this read was neither cached nor classified above
                    auto category_index = classify_read(read_entry, dehost);
                    if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()) {
                        extract_paired_read(category_index, record, record2);
                    }
                }
                if (stats_model_.ready())
                    model_ready_.store(true, std::memory_order_release);
            }
        }
    }

    // Reads are cached while the model is fitted. Stop caching (forcing the model to be used as it
    // is) at either the read-count capacity or the byte budget, whichever comes first.
    [[nodiscard]] bool cache_has_room() const {
        if (cached_reads_.size() >= cached_reads_.capacity())
            return false;
        return cache_byte_budget_ == 0 || cache_bytes_ < cache_byte_budget_;
    }

    // Called with the add_to_cache lock held, once the model is trained: hand the cached reads to
    // classify_pending. The byte count drops with them, as before, since the cache is not refilled.
    void take_cache() {
        pending_reads_.insert(pending_reads_.end(), std::make_move_iterator(cached_reads_.begin()),
                              std::make_move_iterator(cached_reads_.end()));
        cached_reads_.clear();
        cached_reads_.shrink_to_fit();
        cache_bytes_ = 0;
    }

    // Classify the reads taken from the cache, in parallel. Must be called outside any parallel
    // region (after each batch, and from complete). classify_read and the extract functions take
    // their own locks, and the model is only read once trained.
    void classify_pending(const bool dehost = false) {
        if (pending_reads_.empty())
            return;
        PLOG_VERBOSE << "Classify " << pending_reads_.size() << " cached reads";

#pragma omp parallel for num_threads(num_threads_) schedule(dynamic, 64)
        for (size_t i = 0; i < pending_reads_.size(); ++i) {
            auto &read_record = pending_reads_[i];
            auto &read_entry = read_record.read;
            const auto &record = read_record.record;
            auto category_index = classify_read(read_entry, dehost);
            if (run_extract_ and extract_handles_.find(category_index) != extract_handles_.end()) {
                if (read_record.is_paired) {
                    const auto &record2 = read_record.record2;
                    extract_paired_read(category_index, record, record2);
                } else {
                    extract_read(category_index, record);
                }
            }
        }
        pending_reads_.clear();
        pending_reads_.shrink_to_fit();
    }

    void complete(const bool dehost = false) {
        // Input ended before every model collected enough training reads: train what we can
        // (models without enough data fall back to the default KDE) before classifying the cache.
        if (not stats_model_.ready())
            stats_model_.force_ready();
        take_cache();
        classify_pending(dehost);
    }


    void print_summary() const {

        PLOG_INFO << "Results summary: ";
        for (auto i = 0; i < result_summary_.classified_counts.size(); i++) {
            const auto &category = input_summary_.categories.at(i);
            PLOG_INFO << category << " :\t\t" << result_summary_.classified_counts.at(i);
        }
        PLOG_INFO << "unclassified :\t" << result_summary_.unclassified_count;
    }
};


#endif // CHARON_RESULT_H
