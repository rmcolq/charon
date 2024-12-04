#ifndef SIFTER_ENTRY_H
#define SIFTER_ENTRY_H

#pragma once

#include <string>

#include <plog/Log.h>

#include <counts.hpp>
#include <input_summary.hpp>


class ReadEntry
{
    private:
        std::string read_id_;
        uint16_t length_;
        uint32_t num_hashes_{0};
        Counts<uint32_t> counts_;
        std::unordered_map<uint8_t, std::vector<bool>> bits_;
        std::unordered_map<std::string, std::vector<bool>> any_bits_;
        std::unordered_map<std::string, std::vector<bool>> max_bits_;

    public:
        ReadEntry() = default;
        ReadEntry(ReadEntry const &) = default;
        ReadEntry(ReadEntry &&) = default;
        ReadEntry & operator=(ReadEntry const &) = default;
        ReadEntry & operator=(ReadEntry &&) = default;
        ~ReadEntry() = default;

        ReadEntry(const std::string read_id, const uint16_t length, const InputSummary & summary):
            read_id_(read_id),
            length_(length)
            {
                PLOG_DEBUG << "Initialize entry with read_id " << read_id << " and length " << length;
                counts_.set_size(summary.num_categories());
                for (auto i=0; i<summary.num_bins; ++i){
                    bits_[i].reserve(length);
                }
                for (auto i=0; i<summary.num_categories(); ++i){
                    any_bits_[summary.categories.at(i)].reserve(length);
                    max_bits_[summary.categories.at(i)].reserve(length);
                }
            }

        void update_entry(const auto & entry) {
            // this "entry" is a bitvector with a 1 or 0 for each bin in the ibf
            /*PLOG_DEBUG << " entry [";
            for (const auto i : entry){
                PLOG_DEBUG << i;
            }
            PLOG_DEBUG << "]" ;*/

            num_hashes_ += 1;
            for (auto i = 0; i < entry.size(); ++i) {
                const auto row = entry[i];
                //PLOG_DEBUG << "row " << i << " " << row;
                bits_[i].push_back(row);
            }
        };

        void categorize(const InputSummary & summary) {
            PLOG_DEBUG << "categorize";
            for (auto i=0; i<summary.num_categories(); ++i){
                any_bits_[summary.categories.at(i)].resize(num_hashes_, 0);
                max_bits_[summary.categories.at(i)].resize(num_hashes_, 0);
            }
            for (const auto &[bin, bitmap]: bits_) {
                const auto category = summary.bin_to_category.at(bin);
                PLOG_DEBUG << +bin << " belongs to " << category;

                PLOG_DEBUG << num_hashes_ << " " << any_bits_[category].size() << " " << bitmap.size();
                for (auto j = 0; j < num_hashes_; ++j) {
                    any_bits_[category][j] = bitmap[j] | any_bits_[category][j];
                }
                PLOG_DEBUG << "done any bits";
                auto current_bit_count = std::count(bitmap.begin(), bitmap.end(), true);
                const auto &max_bitmap = max_bits_[category];
                auto max_bit_count = std::count(max_bitmap.begin(), max_bitmap.end(), true);
                PLOG_DEBUG << current_bit_count << " " << max_bit_count;
                if (current_bit_count > max_bit_count) {
                    PLOG_DEBUG << "redefine max";
                    max_bits_[category] = bitmap;
                }
            }
        };

        void collect_counts(const InputSummary & summary)
        {
            PLOG_DEBUG << "collect_counts for read_id " << read_id_;
            for (auto i=0; i<summary.num_categories(); ++i){
                for (auto j=0; j<=i; ++j){
                    const auto row = max_bits_[summary.categories.at(i)];
                    const auto col = max_bits_[summary.categories.at(j)];
                    for (auto k=0; k<num_hashes_; ++k){
                        PLOG_DEBUG << "(" << i << "," << j << ")" << k << row.size() << col.size();
                        if (row[k] and col[k]){
                            counts_(i,j) += 1;
                        }
                    }
                    PLOG_DEBUG << "done hashes";
                }
                PLOG_DEBUG << "done category " << i;
            }
            PLOG_DEBUG << "done";
        };

        void print_result(const InputSummary & summary){
            std::cout << read_id_ << "\t" << num_hashes_ << "\t";
            for (auto i=0; i<summary.num_categories(); i++){
                std::cout << summary.categories.at(i) << ":" << counts_(i,i) << "\t";
            }
            for (auto i=0; i<summary.num_categories(); i++){
                for (auto j=0; j<i; j++) {
                    std::cout << summary.categories.at(i) << "x" << summary.categories.at(j) << ":" << counts_(i, j) << "\t";
                }
            }
            /*for (auto i=0; i<bits_.size(); i++){
                for (auto j=0; j<bits_[i].size(); j++) {
                    std::cout << +bits_[i][j];
                }
                std::cout << "\t";
            }*/
            for (auto i=0; i<summary.num_categories(); i++){
                const auto & category = summary.categories.at(i);
                const auto & bitvector = max_bits_.at(category);
                for (auto j=0; j<bitvector.size(); j++) {
                    std::cout << +bitvector.at(j);
                }
                std::cout << "\t";
            }
            std::cout << std::endl;
        };
    };

#endif // SIFTER_ENTRY_H
