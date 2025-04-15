#ifndef CHARON_ENTRY_H
#define CHARON_ENTRY_H

#pragma once

#include <string>
#include <algorithm>

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
    std::vector<seqan3::interleaved_bloom_filter< seqan3::compressed >::membership_agent_type::binning_bitvector> bits_;// this collects over all bins
    std::unordered_map<uint8_t, std::vector<bool>> max_bits_; // this summarizes over categories (which may have multiple bins)
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

    ~ReadEntry() = default;

    ReadEntry(const std::string& read_id, const uint32_t& length, const float& mean_quality, const float& compression, const InputSummary &summary) :
            read_id_(read_id),
            length_(length),
            mean_quality_(mean_quality),
            compression_(compression),
            counts_(summary.num_categories(), 0),
            proportions_(summary.num_categories(),0),
            unique_proportions_(summary.num_categories(),0),
            unique_counts_(summary.num_categories(), 0),
            probabilities_(summary.num_categories(),1)
    {
        PLOG_DEBUG << "Initialize entry with read_id " << read_id << " and length " << length;
        bits_.reserve(length);

        for (auto i = 0; i < summary.num_categories(); ++i) {
            max_bits_[i].reserve(length);
        }
        PLOG_VERBOSE << "Initializing complete for read_id " << read_id;
    }

    const std::string& read_id() const {
        return read_id_;
    }

    const std::vector<float>& proportions() const {
        return proportions_;
    }

    const std::vector<float>& unique_proportions() const {
        return unique_proportions_;
    }

    const uint8_t call() const {
        return call_;
    }

    const uint8_t confidence_score() const {
        return confidence_score_;
    }

    void update_entry(const auto &entry) {
        // this "entry" is a bitvector with a 1 or 0 for each bin in the ibf
        bits_.emplace_back(entry);
        num_hashes_ += 1;
    };

    void get_counts(const InputSummary &summary) {
        PLOG_DEBUG << "Get max bits per category for read " << read_id_;

        // get totals in each bin
        seqan3::counting_vector<uint64_t> total_bits_per_bin(summary.num_bins, 0);
        for (const auto & entry : bits_)
        {
            total_bits_per_bin += entry;
        }
        PLOG_DEBUG << "Have the following total hits per bin: " << total_bits_per_bin;

        // identify max bin per category
        std::vector<uint8_t> index_per_category(summary.num_categories(),std::numeric_limits<uint8_t>::max());
        for ( auto bin=0; bin < total_bits_per_bin.size(); ++bin)
        {
            const auto & bits_in_bin = total_bits_per_bin[bin];
            const auto & category = summary.bin_to_category.at(bin);
            const auto & index = summary.category_index(category);
            if (index_per_category.at(index) == std::numeric_limits<uint8_t>::max() or total_bits_per_bin[bin] > total_bits_per_bin[index_per_category[index]]) {
                index_per_category.at(index) = bin;
                counts_.at(index) = total_bits_per_bin[bin];
                PLOG_DEBUG << "Category " << category << " has max index " << +bin;
            }
        }

        //
        for (auto category = 0; category < summary.num_categories(); ++category) {
            max_bits_.at(category).resize(num_hashes_);
        }

        // collect together the bitvector for max bits and the unique_counts
        std::vector<uint8_t> found;
        for (const auto & entry : bits_){
            found.clear();
            for (auto category=0; category < index_per_category.size(); ++category)
            {
                const auto & category_index = index_per_category[category];
                max_bits_.at(category).emplace_back(entry[category_index]);
                if (entry[category_index] == 1)
                {
                    found.push_back(category);
                }
            }
            if (found.size() == 1) {
                auto i = found.front();
                unique_counts_.at(i) += 1;
            }
        }
        PLOG_DEBUG << "Found unique counts " << unique_counts_;
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

    void post_process(const InputSummary &summary) {
        get_counts(summary);
        get_proportions();
    }

    void call_category(const StatsModel &stats_model) {
        //TODO extend this to work with more than 2 categories
        assert(probabilities_.size() == 2);

        double first = probabilities_.at(0);
        uint8_t first_pos = 0;
        double second = probabilities_.at(1);
        uint8_t second_pos = 1;
        if (unique_counts_.at(second_pos) > unique_counts_.at(first_pos))
        {
            std::swap(first, second);
            std::swap(first_pos, second_pos);
        }

        /*for (auto i = 0; i < probabilities_.size(); ++i) {
            const double &val = probabilities_.at(i);
            if (val >= second) {
                second = val;
                second_pos = i;
                if (second >= first)
                {
                    std::swap(first, second);
                    std::swap(first_pos, second_pos);
                }
            }
        }*/

        auto raw_confidence = unique_counts_.at(first_pos) - unique_counts_.at(second_pos);
        if (raw_confidence > std::numeric_limits<uint8_t>::max())
            confidence_score_ = std::numeric_limits<uint8_t>::max();
        else
            confidence_score_ = static_cast<uint8_t>(raw_confidence);

        if (mean_quality_ < stats_model.min_quality()){
            return;
        }

        if (length_ < stats_model.min_length()){
            return;
        }

        if (compression_ < stats_model.min_compression()){
            return;
        }

        if (second == 0 and first > 0) {
            call_ = first_pos;
        } else {
            /*auto ratio = first/second;
            PLOG_VERBOSE << "confidence score " << ratio << " from " << first << "/" << second;
            if (ratio < std::numeric_limits<uint8_t>::max())
                confidence_score_ = static_cast<uint8_t>(ratio);
            else
                confidence_score_ = std::numeric_limits<uint8_t>::max();*/
            if (confidence_score_ > stats_model.confidence_threshold() and first > second)
                call_ = first_pos;
        }

        const auto & first_count = counts_.at(first_pos);
        const auto & second_count = counts_.at(second_pos);
        if (second_count > first_count or first_count - second_count < stats_model.min_num_hits()){
            PLOG_DEBUG << read_id_ << " has first count " << +first_count << " and second count " << +second_count << " which have differences less than " << +stats_model.min_num_hits();
            call_ = std::numeric_limits<uint8_t>::max(); // if we don't see at least this number of hits difference, then no call
        }

        const auto & first_prop = proportions_.at(first_pos);
        const auto & second_prop = proportions_.at(second_pos);
        if (second_prop > first_prop or first_prop - second_prop < stats_model.min_proportion_difference()){
            PLOG_DEBUG << read_id_ << " has first proportion " << first_prop << " and second proportion " << second_prop << " which have differences less than " << stats_model.min_proportion_difference();
            call_ = std::numeric_limits<uint8_t>::max(); // if we don't see at least this level of discrepancy between the proportion hitting against each database don't call
        }
    }

    void call_host(const StatsModel &stats_model, const uint8_t host_index) {
        assert(probabilities_.size() == 2);
        const uint8_t other_index = 1-host_index;

        double host_unique_prop = unique_proportions_.at(host_index);
        double other_unique_prop = unique_proportions_.at(other_index);

        auto first_pos = host_index;
        auto second_pos = other_index;
        if (host_unique_prop < other_unique_prop)
        {
            first_pos = other_index;
            second_pos = host_index;
        }
        auto raw_confidence = unique_counts_.at(first_pos) - unique_counts_.at(second_pos);
        if (raw_confidence > std::numeric_limits<uint8_t>::max())
            confidence_score_ = std::numeric_limits<uint8_t>::max();
        else
            confidence_score_ = static_cast<uint8_t>(raw_confidence);
        if (confidence_score_ < stats_model.confidence_threshold()){
            return;
        }

        if (mean_quality_ < stats_model.min_quality()){
            return;
        }

        if (length_ < stats_model.min_length()){
            return;
        }

        if (compression_ < stats_model.min_compression()){
            return;
        }

        if (host_unique_prop > other_unique_prop and host_unique_prop - other_unique_prop > stats_model.min_proportion_difference())
            call_ = host_index;
        else if (host_unique_prop < stats_model.host_unique_prop_lo_threshold() and host_unique_prop < other_unique_prop and other_unique_prop - host_unique_prop > stats_model.min_proportion_difference())
            call_ = other_index;
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
        call_host(stats_model, host_index);
    }

    void print_result(const InputSummary &summary) {
        std::cout << read_id_ << "\t" << num_hashes_ << "\t" << summary.category_name(call_) << "\t" << +confidence_score_ << "\t" ;
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
        for (auto i = 0; i < summary.num_categories(); i++) {
            //const auto & category = summary.categories.at(i);
            const auto &bitvector = max_bits_.at(i);
            for (auto j = 0; j < bitvector.size(); j++) {
                std::cout << +bitvector.at(j);
            }
            std::cout << "\t";
        }
        std::cout << std::endl;
    };

    void print_assignment_result(const InputSummary &summary) const {
        // mimic the kraken assignment format with tab separated columns classification status, read_id, call, length, num_hashes, details
        if (call_ == std::numeric_limits<uint8_t>::max())
            std::cout << "U" << "\t";
        else
            std::cout << "C" << "\t";
        std::cout.precision(6);
        std::cout << read_id_ << "\t" << summary.category_name(call_) << "\t" << length_ << "\t" << num_hashes_ << "\t" << mean_quality_ << "\t" << +confidence_score_ << "\t" << compression_ << "\t";
        for (auto i = 0; i < summary.num_categories(); i++) {
            std::cout << summary.categories.at(i) << ":" << counts_.at(i) << ":" << proportions_.at(i)
                      << ":" << unique_proportions_.at(i) << ":" << probabilities_.at(i) << " ";
        }

        std::cout << std::endl;
    };
};

#endif // CHARON_ENTRY_H
