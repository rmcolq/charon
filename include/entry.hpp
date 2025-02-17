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
    float mean_quality_;

    uint32_t num_hashes_{0};
    std::unordered_map<uint8_t, std::vector<bool>> bits_; // this collects over all bins
    std::unordered_map<uint8_t, std::vector<bool>> max_bits_; // this summarizes over categories (which may have multiple bins)
    Counts<uint32_t> counts_;
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

    ReadEntry(const std::string& read_id, const uint32_t& length, const double& mean_quality, const InputSummary &summary) :
            read_id_(read_id),
            mean_quality_(mean_quality),
            proportions_(summary.num_categories(),0),
            unique_proportions_(summary.num_categories(),0),
            unique_counts_(summary.num_categories(), 0),
            probabilities_(summary.num_categories(),1)
    {
        PLOG_DEBUG << "Initialize entry with read_id " << read_id << " and length " << length;
        counts_.set_size(summary.num_categories());
        for (auto i = 0; i < summary.num_bins; ++i) {
            bits_[i].reserve(length);
        }

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
        /*PLOG_DEBUG << " entry [";
        for (const auto i : entry){
            PLOG_DEBUG << i;
        }
        PLOG_DEBUG << "]" ;*/

        assert(entry.size() == bits_.size());
        num_hashes_ += 1;
        for (auto bucket = 0; bucket < entry.size(); ++bucket) {
            bits_.at(bucket).emplace_back(entry[bucket]);
            if (bits_.at(bucket).size() != num_hashes_){
                PLOG_ERROR << read_id_ << " bucket k=" << +bucket << " out of " << entry.size() << " final bin/hash sizes " << bits_.at(bucket).size() << " and " << num_hashes_;
            }
            assert(bits_.at(bucket).size() == num_hashes_);
        }
    };

    void get_max_bits(const InputSummary &summary) {
        PLOG_DEBUG << "Get max bits per category for read " << read_id_;
        const auto num_categories = max_bits_.size();
        for (auto i = 0; i < num_categories; ++i) {
            max_bits_.at(i).resize(num_hashes_, 0);
        }
        for (const auto &[bin, bitmap]: bits_) {
            assert(bitmap.size() == num_hashes_);
            const auto & category = summary.bin_to_category.at(bin);
            const auto & index = summary.category_index(category);
            assert(index < num_categories);
            PLOG_VERBOSE << +bin << " belongs to " << category << " with index " << +index << " for read " << read_id_;

            auto current_bit_count = std::count(bitmap.begin(), bitmap.end(), true);
            const auto &max_bitmap = max_bits_.at(index);
            auto max_bit_count = std::count(max_bitmap.begin(), max_bitmap.end(), true);
            PLOG_VERBOSE << current_bit_count << " " << max_bit_count << " for read " << read_id_;
            if (current_bit_count > max_bit_count) {
                PLOG_VERBOSE << "redefine max for read " << read_id_;
                max_bits_.at(index) = bitmap;
            }
        }
    };

    void get_counts() {
        PLOG_DEBUG << "Collect counts for read_id " << read_id_ << " which has " << max_bits_.size() << " categories";
        const auto num_categories = max_bits_.size();
        for (auto i = 0; i < num_categories; ++i) {
            for (auto j = 0; j <= i; ++j) {
                const auto & row = max_bits_.at(i);
                const auto & col = max_bits_.at(j);
                assert(row.size() == num_hashes_ and col.size() == num_hashes_);
                for (auto k = 0; k < num_hashes_; ++k) {
                    //PLOG_DEBUG << "(" << i << "," << j << ")" << k << row.size() << col.size();
                    if (row.at(k) and col.at(k)) {
                        counts_(i, j) += 1;
                    }
                }
                //PLOG_DEBUG << "done hashes";
            }
            //PLOG_DEBUG << "done category " << i;
        }
        //PLOG_DEBUG << "done";
    };

    void get_proportions() {
        PLOG_DEBUG << "Collect proportions for read_id " << read_id_;
        const auto num_categories = proportions_.size();
        for (auto i = 0; i < num_categories; ++i) {
            proportions_.at(i) = static_cast< float >(counts_(i, i)) / static_cast< float >(num_hashes_);
        }
        return;
    };

    void get_unique_counts() {
        PLOG_DEBUG << "Collect unique counts for read_id " << read_id_;
        const auto num_categories = max_bits_.size();
        std::vector<uint8_t> found;
        for (auto k = 0; k < num_hashes_; ++k) {
            found.clear();
            for (auto i = 0; i < num_categories; ++i) {
                const auto & row = max_bits_.at(i);
                if (k > row.size())
                    PLOG_ERROR << "k > row.size()";
                if (row.at(k)) {
                    found.push_back(i);
                }
            }
            if (found.size() == 1) {
                auto i = found.front();
                unique_counts_.at(i) += 1;
            }
        }

    };

    void get_unique_proportions() {
        PLOG_DEBUG << "Collect unique proportions for read_id " << read_id_;
        const auto num_categories = max_bits_.size();
        for (auto i = 0; i < num_categories; ++i) {
            unique_proportions_.at(i) = static_cast< float >(unique_counts_.at(i)) / static_cast< float >(num_hashes_);
        }
        return;
    };

    void post_process(const InputSummary &summary) {
        get_max_bits(summary);
        get_counts();
        get_unique_counts();
        get_unique_proportions();
        get_proportions();
    }

    void call_category(const float& min_quality, const uint8_t& confidence_threshold, const uint8_t& min_num_hits, const float& min_proportion_difference) {
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

        if (mean_quality_ < min_quality){
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
            if (confidence_score_ > confidence_threshold and first > second)
                call_ = first_pos;
        }

        const auto & first_count = counts_(first_pos,first_pos);
        const auto & second_count = counts_(second_pos,second_pos);
        if (second_count > first_count or first_count - second_count < min_num_hits){
            PLOG_DEBUG << read_id_ << " has first count " << +first_count << " and second count " << +second_count << " which have differences less than " << +min_num_hits;
            call_ = std::numeric_limits<uint8_t>::max(); // if we don't see at least this number of hits difference, then no call
        }

        const auto & first_prop = proportions_.at(first_pos);
        const auto & second_prop = proportions_.at(second_pos);
        if (second_prop > first_prop or first_prop - second_prop < min_proportion_difference){
            PLOG_DEBUG << read_id_ << " has first proportion " << first_prop << " and second proportion " << second_prop << " which have differences less than " << min_proportion_difference;
            call_ = std::numeric_limits<uint8_t>::max(); // if we don't see at least this level of discrepancy between the proportion hitting against each database don't call
        }
    }

    void classify(const StatsModel &stats_model) {
        PLOG_DEBUG << "Classify read " << read_id_;;
        for (auto i = 0; i < unique_proportions_.size(); ++i) {
            const auto &read_proportion = unique_proportions_.at(i);
            const auto result_pair = stats_model.classify(i, read_proportion);
            PLOG_DEBUG << "Pos " << +i << " has read proportion " << read_proportion << " yielding probs "
                       << result_pair.pos << " and " << result_pair.neg << " for read " << read_id_;;
            for (auto j = 0; j < unique_proportions_.size(); ++j) {
                if (i == j)
                    probabilities_.at(j) *= result_pair.pos;
                //else
                //    probabilities_.at(j) *= result_pair.g_neg;
            }
        }
        call_category(stats_model.min_quality(), stats_model.confidence_threshold(), stats_model.min_num_hits(), stats_model.min_proportion_difference());
    }

    void print_result(const InputSummary &summary) {
        std::cout << read_id_ << "\t" << num_hashes_ << "\t" << summary.category_name(call_) << "\t" << +confidence_score_ << "\t" ;
        for (auto i = 0; i < summary.num_categories(); i++) {
            std::cout << summary.categories.at(i) << ":" << counts_(i, i) << ":" << probabilities_.at(i) << "\t";
        }
        for (auto i = 0; i < summary.num_categories(); i++) {
            for (auto j = 0; j < i; j++) {
                std::cout << summary.categories.at(i) << "_x_" << summary.categories.at(j) << ":" << counts_(i, j)
                          << "\t";
            }
        }
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
        // mimic the kraken assignment format with tab separated columns classification status, read_id call, num_hashes, details
        if (call_ == std::numeric_limits<uint8_t>::max())
            std::cout << "U" << "\t";
        else
            std::cout << "C" << "\t";
        std::cout.precision(6);
        std::cout << read_id_ << "\t" << summary.category_name(call_) << "\t" << num_hashes_ << "\t" << mean_quality_ << "\t" << +confidence_score_ << "\t" ;
        for (auto i = 0; i < summary.num_categories(); i++) {
            std::cout << summary.categories.at(i) << ":" << counts_(i, i) << ":" << proportions_.at(i)
                      << ":" << unique_proportions_.at(i) << ":" << probabilities_.at(i) << " ";
        }

        std::cout << std::endl;
    };
};

#endif // CHARON_ENTRY_H
