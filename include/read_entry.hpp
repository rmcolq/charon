#ifndef CHARON_ENTRY_H
#define CHARON_ENTRY_H

#pragma once

#include <sstream>
#include <string>
#include <algorithm>
#include <bit>
#include <ankerl/unordered_dense.h>

#include <plog/Log.h>

#include <counts.hpp>
#include <input_summary.hpp>
#include <classify_stats.hpp>


class ReadEntry {
private:
    std::string read_id_;
    uint32_t length_;
    float mean_quality_;
    float compression_;

    uint32_t num_hashes_{0};
    // IBF hit bits for every hash, stored flat: hash i occupies words [i * bin_words_, (i + 1) * bin_words_)
    std::vector<uint64_t> bits_;
    size_t bin_words_{0};
    std::vector<uint32_t> counts_;
    std::vector<uint32_t> unique_counts_;
    std::vector<float> proportions_; // this collects over categories the proportion of all hashes which were from the given category
    std::vector<float> unique_proportions_; // this collects over categories the proportion of all hashes which were unique to the given category
    std::vector<double> probabilities_; // this collects over categories the probability of this read given the data come from that category
    uint8_t call_ = std::numeric_limits<uint8_t>::max();
    uint8_t confidence_score_ = 0;
public:
    ReadEntry() = default;

    ReadEntry(ReadEntry const &) = default;

    ReadEntry(ReadEntry &&) = default;

    ReadEntry &operator=(ReadEntry const &) = default;

    ReadEntry &operator=(ReadEntry &&) = default;

    ReadEntry(const std::string &read_id, const uint32_t &length, const float &mean_quality, const float &compression,
              const InputSummary &summary) :
            read_id_(read_id),
            length_(length),
            mean_quality_(mean_quality),
            compression_(compression),
            counts_(summary.num_categories(), 0),
            proportions_(summary.num_categories(), 0),
            unique_proportions_(summary.num_categories(), 0),
            unique_counts_(summary.num_categories(), 0),
            probabilities_(summary.num_categories(), 1) {
        PLOG_DEBUG << "Initialize entry with read_id " << read_id << " and length " << length;
        PLOG_VERBOSE << "Initializing complete for read_id " << read_id;
    }

    // Explicit destructor to ensure proper cleanup (RAII)
    ~ReadEntry() = default;

    [[nodiscard]] const std::string &read_id() const {
        return read_id_;
    }

    [[nodiscard]] const std::vector<float> &proportions() const {
        return proportions_;
    }

    [[nodiscard]] const std::vector<float> &unique_proportions() const {
        return unique_proportions_;
    }

    [[nodiscard]] const uint8_t call() const {
        return call_;
    }

    [[nodiscard]] const uint8_t confidence_score() const {
        return confidence_score_;
    }

    [[nodiscard]] uint32_t num_hashes() const {
        return num_hashes_;
    }

    [[nodiscard]] uint32_t length() const {
        return length_;
    }

    void update_entry(const auto &entry) {
        // this "entry" is a bitvector with a 1 or 0 for each bin in the ibf
        const auto &bitvector = entry.raw_data();
        if (bin_words_ == 0) {
            bin_words_ = (bitvector.size() + 63) / 64;
            bits_.reserve((length_ / 8 + 16) * bin_words_);
        }
        bits_.insert(bits_.end(), bitvector.data(), bitvector.data() + bin_words_);
        num_hashes_ += 1;
    };

    [[nodiscard]] bool bit(const size_t hash, const size_t bin) const {
        return (bits_[hash * bin_words_ + (bin >> 6)] >> (bin & 63)) & 1u;
    }

    void get_counts(const InputSummary &summary) {
        PLOG_DEBUG << "Get max bits per category for read " << read_id_;

        // get totals in each bin
        seqan3::counting_vector<uint64_t> total_bits_per_bin(summary.num_bins, 0);
        for (size_t i = 0; i < num_hashes_; ++i) {
            for (size_t word = 0; word < bin_words_; ++word) {
                uint64_t remaining = bits_[i * bin_words_ + word];
                while (remaining) {
                    const size_t bin = (word << 6) + std::countr_zero(remaining);
                    if (bin < summary.num_bins)
                        total_bits_per_bin[bin] += 1;
                    remaining &= remaining - 1;
                }
            }
        }
        PLOG_DEBUG << "Have the following total hits per bin: " << total_bits_per_bin;

        // identify max bin per category
        std::vector<uint8_t> index_per_category(summary.num_categories(), std::numeric_limits<uint8_t>::max());
        for (auto bin = 0; bin < total_bits_per_bin.size(); ++bin) {
            const auto &bits_in_bin = total_bits_per_bin[bin];
            const auto &category = summary.bin_to_category.at(bin);
            const auto &index = summary.category_index(category);
            if (index_per_category.at(index) == std::numeric_limits<uint8_t>::max() or
                total_bits_per_bin[bin] > total_bits_per_bin[index_per_category[index]]) {
                index_per_category.at(index) = bin;
                counts_.at(index) = total_bits_per_bin[bin];
                PLOG_DEBUG << "Category " << category << " has max index " << +bin;
            }
        }

        // count hashes found in exactly one category's best bin
        for (size_t i = 0; i < num_hashes_; ++i) {
            uint8_t num_found = 0;
            uint8_t found_category = 0;
            for (auto category = 0; category < index_per_category.size(); ++category) {
                if (bit(i, index_per_category[category])) {
                    ++num_found;
                    found_category = category;
                }
            }
            if (num_found == 1) {
                unique_counts_.at(found_category) += 1;
            }
        }
        PLOG_DEBUG << "Found unique counts " << unique_counts_;

        // The per-hash, per-bin bit vectors are only needed to derive the counts above. They are by far
        // the largest part of an entry (one heap-allocated vector per hash), and entries are kept in
        // memory while the model is fitted, so release them now.
        decltype(bits_)().swap(bits_);
    };

    void get_proportions() {
        PLOG_DEBUG << "Collect proportions for read_id " << read_id_;
        const auto num_categories = proportions_.size();
        for (auto i = 0; i < num_categories; ++i) {
            proportions_.at(i) = static_cast< float >(counts_.at(i)) / static_cast< float >(num_hashes_);
            unique_proportions_.at(i) = static_cast< float >(unique_counts_.at(i)) / static_cast< float >(num_hashes_);
        }
        PLOG_DEBUG << "Found proportions " << proportions_;
        PLOG_DEBUG << "Found unique_proportions " << unique_proportions_;
        return;
    }

    // Alternative to get_counts: a hash counts for a category if ANY bin of that category contains it.
    void get_counts_category_or(const InputSummary &summary) {
        const auto num_categories = summary.num_categories();
        std::vector<uint64_t> masks(num_categories * bin_words_, 0);
        for (size_t bin = 0; bin < summary.num_bins; ++bin) {
            const auto category = summary.category_index(summary.bin_to_category.at(bin));
            masks[category * bin_words_ + (bin >> 6)] |= uint64_t{1} << (bin & 63);
        }
        for (size_t i = 0; i < num_hashes_; ++i) {
            uint8_t num_found = 0;
            uint8_t found_category = 0;
            for (uint8_t category = 0; category < num_categories; ++category) {
                bool hit = false;
                for (size_t word = 0; word < bin_words_; ++word)
                    hit |= (bits_[i * bin_words_ + word] & masks[category * bin_words_ + word]) != 0;
                if (hit) {
                    counts_.at(category) += 1;
                    ++num_found;
                    found_category = category;
                }
            }
            if (num_found == 1)
                unique_counts_.at(found_category) += 1;
        }
    };

    void post_process(const InputSummary &summary, const bool category_or = false) {
        if (category_or)
            get_counts_category_or(summary);
        else
            get_counts(summary);
        get_proportions();
    }

    void call_category(const StatsModel &stats_model) {
        // We want to compare based on unique_counts_ because confidence is defined based on the difference of the
        // top 2 unique_counts_. However, we will reject the result if the probabilities do not support this ordering.
        uint8_t first_pos = 0;
        uint8_t second_pos = 1;
        if (unique_counts_.at(second_pos) > unique_counts_.at(first_pos)) {
            std::swap(first_pos, second_pos);
        }

        for (auto i = 2; i < unique_counts_.size(); ++i) {
            if (unique_counts_.at(i) > unique_counts_.at(second_pos)) {
                second_pos = i;
                if (unique_counts_.at(second_pos) > unique_counts_.at(first_pos)) {
                    std::swap(first_pos, second_pos);
                }
            }
        }

        auto raw_confidence = unique_counts_.at(first_pos) - unique_counts_.at(second_pos);
        if (raw_confidence > std::numeric_limits<uint8_t>::max())
            confidence_score_ = std::numeric_limits<uint8_t>::max();
        else
            confidence_score_ = static_cast<uint8_t>(raw_confidence);

        if (mean_quality_ < stats_model.min_quality()) {
            return;
        }

        if (length_ < stats_model.min_length()) {
            return;
        }

        if (compression_ < stats_model.min_compression()) {
            return;
        }

        if (probabilities_.at(second_pos) == 0 and probabilities_.at(first_pos) > 0) {
            call_ = first_pos;
        } else {
            if (confidence_score_ > stats_model.confidence_threshold() and
                probabilities_.at(first_pos) > probabilities_.at(second_pos))
                call_ = first_pos;
        }

        const auto &first_count = counts_.at(first_pos);
        const auto &second_count = counts_.at(second_pos);
        if (second_count > first_count or first_count - second_count < stats_model.min_num_hits()) {
            PLOG_DEBUG << read_id_ << " has first count " << +first_count << " and second count " << +second_count
                       << " which have differences less than " << +stats_model.min_num_hits();
            call_ = std::numeric_limits<uint8_t>::max(); // if we don't see at least this number of hits difference, then no call
        }

        const auto &first_prop = proportions_.at(first_pos);
        const auto &second_prop = proportions_.at(second_pos);
        if (second_prop > first_prop or first_prop - second_prop < stats_model.min_proportion_difference()) {
            PLOG_DEBUG << read_id_ << " has first proportion " << first_prop << " and second proportion " << second_prop
                       << " which have differences less than " << stats_model.min_proportion_difference();
            call_ = std::numeric_limits<uint8_t>::max(); // if we don't see at least this level of discrepancy between the proportion hitting against each database don't call
        }
    }

    void call_host(const StatsModel &stats_model, const uint8_t host_index) {
        assert(probabilities_.size() == 2);
        const uint8_t other_index = 1 - host_index;

        double host_unique_prop = unique_proportions_.at(host_index);
        double other_unique_prop = unique_proportions_.at(other_index);
        double host_prob = probabilities_.at(host_index);
        double other_prob = probabilities_.at(other_index);

        auto first_pos = host_index;
        auto second_pos = other_index;
        if (host_unique_prop < other_unique_prop) {
            first_pos = other_index;
            second_pos = host_index;
        }
        auto raw_confidence = unique_counts_.at(first_pos) - unique_counts_.at(second_pos);
        if (raw_confidence > std::numeric_limits<uint8_t>::max())
            confidence_score_ = std::numeric_limits<uint8_t>::max();
        else
            confidence_score_ = static_cast<uint8_t>(raw_confidence);
        if (confidence_score_ < stats_model.confidence_threshold()) {
            return;
        }

        if (mean_quality_ < stats_model.min_quality()) {
            return;
        }

        if (length_ < stats_model.min_length()) {
            return;
        }

        if (compression_ < stats_model.min_compression()) {
            return;
        }

        if (host_unique_prop > other_unique_prop
            and host_unique_prop - other_unique_prop > stats_model.min_proportion_difference()
            and host_prob > other_prob
            and host_prob - other_prob > stats_model.min_prob_difference()
            and std::max(host_prob * confidence_score_, static_cast<double>(confidence_score_)) >=
                stats_model.confidence_probability_threshold())
            call_ = host_index;
        else if (host_unique_prop < stats_model.host_unique_prop_lo_threshold()
                 and host_unique_prop < other_unique_prop
                 and other_unique_prop - host_unique_prop > stats_model.min_proportion_difference()
                 and host_prob < other_prob
                 and other_prob - host_prob > stats_model.min_prob_difference()
                 and std::max(other_prob * confidence_score_, static_cast<double>(confidence_score_)) >=
                     stats_model.confidence_probability_threshold())
            call_ = other_index;
    }

    // Host-only index (a single host category, no microbial category). There is nothing to compare the host against,
    // so the read is called host when the model supports it being from the host distribution rather than the
    // negative one; every other read is left uncalled, i.e. it is not host. This is the host-only analogue of the
    // first branch of call_host, with the missing category having no unique hits and probability 1 - host_prob.
    void call_host_only(const StatsModel &stats_model, const uint8_t host_index) {
        assert(probabilities_.size() == 1);

        const double host_unique_prop = unique_proportions_.at(host_index);
        const double host_prob = probabilities_.at(host_index);
        const double other_prob = 1 - host_prob;

        const auto raw_confidence = unique_counts_.at(host_index);
        if (raw_confidence > std::numeric_limits<uint8_t>::max())
            confidence_score_ = std::numeric_limits<uint8_t>::max();
        else
            confidence_score_ = static_cast<uint8_t>(raw_confidence);
        if (confidence_score_ < stats_model.confidence_threshold()) {
            return;
        }

        if (mean_quality_ < stats_model.min_quality()) {
            return;
        }

        if (length_ < stats_model.min_length()) {
            return;
        }

        if (compression_ < stats_model.min_compression()) {
            return;
        }

        if (host_unique_prop > stats_model.min_proportion_difference()
            and host_prob > other_prob
            and host_prob - other_prob > stats_model.min_prob_difference()
            and std::max(host_prob * confidence_score_, static_cast<double>(confidence_score_)) >=
                stats_model.confidence_probability_threshold())
            call_ = host_index;
    }

    void apply_model(const StatsModel &stats_model) {
        for (auto i = 0; i < unique_proportions_.size(); ++i) {
            const auto &read_proportion = unique_proportions_.at(i);
            const auto result_pair = stats_model.classify(i, read_proportion);
            PLOG_DEBUG << "Pos " << +i << " has read proportion " << read_proportion << " yielding probs "
                       << result_pair.pos << " and " << result_pair.neg << " for read " << read_id_;
            probabilities_.at(i) *= result_pair.pos;
        }
    }

    void classify(const StatsModel &stats_model) {
        PLOG_DEBUG << "Classify read " << read_id_;
        apply_model(stats_model);
        call_category(stats_model);
    }

    void dehost(const StatsModel &stats_model, const uint8_t host_index) {
        PLOG_DEBUG << "Classify read " << read_id_;
        apply_model(stats_model);
        if (probabilities_.size() == 1)
            call_host_only(stats_model, host_index);
        else
            call_host(stats_model, host_index);
    }

    void print_result(const InputSummary &summary) {
        std::cout << read_id_ << "\t" << num_hashes_ << "\t" << summary.category_name(call_) << "\t"
                  << +confidence_score_ << "\t";
        for (auto i = 0; i < summary.num_categories(); i++) {
            std::cout << summary.categories.at(i) << ":" << counts_.at(i) << ":" << probabilities_.at(i) << "\t";
        }
        /*for (auto i = 0; i < summary.num_categories(); i++) {
            for (auto j = 0; j < i; j++) {
                std::cout << summary.categories.at(i) << "_x_" << summary.categories.at(j) << ":" << counts_(i, j)
                          << "\t";
            }
        }*/
        /*for (auto i=0; i<bits_.size(); i++){
            for (auto j=0; j<bits_.at(i).size(); j++) {
                std::cout << +bits_.at(i)[j];
            }
            std::cout << "\t";
        }*/
        std::cout << '\n';
    };

    // Formats the line printed by print_assignment_result, so it can be built outside any lock.
    [[nodiscard]] std::string assignment_line(const InputSummary &summary) const {
        // mimic the kraken assignment format with tab separated columns classification status, read_id, call, length, num_hashes, details
        std::ostringstream out;
        if (call_ == std::numeric_limits<uint8_t>::max())
            out << "U" << "\t";
        else
            out << "C" << "\t";
        out.precision(6);
        out << read_id_ << "\t" << summary.category_name(call_) << "\t" << length_ << "\t" << num_hashes_ << "\t"
            << mean_quality_ << "\t" << +confidence_score_ << "\t" << compression_ << "\t";
        for (auto i = 0; i < summary.num_categories(); i++) {
            out << summary.categories.at(i) << ":" << counts_.at(i) << ":" << proportions_.at(i)
                << ":" << unique_proportions_.at(i) << ":" << probabilities_.at(i) << " ";
        }
        out << '\n';
        return out.str();
    }

    void print_assignment_result(const InputSummary &summary) const {
        std::cout << assignment_line(summary);
    };
};

#endif // CHARON_ENTRY_H
