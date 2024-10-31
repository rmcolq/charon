#ifndef SIFTER_ENTRY_H
#define SIFTER_ENTRY_H

#pragma once

#include <string>

#include <plog/Log.h>

#include <counts.hpp>

class ReadEntry
{
    private:
        std::string read_id_;
        uint16_t length_;
        uint32_t num_hashes_{0};
        Counts<uint32_t> counts_;
        std::vector<std::vector<bool>> bits_;

    public:
        ReadEntry() = default;
        ReadEntry(ReadEntry const &) = default;
        ReadEntry(ReadEntry &&) = default;
        ReadEntry & operator=(ReadEntry const &) = default;
        ReadEntry & operator=(ReadEntry &&) = default;
        ~ReadEntry() = default;

        ReadEntry(const std::string read_id, const uint16_t length, const uint8_t num_bins):
            read_id_(read_id),
            length_(length)
            {counts_.set_size(num_bins);}

        void update_entry(const auto & entry){
            /*PLOG_DEBUG << " entry [";
            for (const auto i : entry){
                PLOG_DEBUG << i;
            }
            PLOG_DEBUG << "]" ;*/

            num_hashes_ += 1;
            for (auto i=0; i<entry.size(); ++i)
            {
                const auto row = entry[i];
                PLOG_DEBUG << "row " << i << " " << row;
                bits_[i].push_back(row);
                if (row == 0){
                    continue;
                }
                for (auto j=0; j<=i; ++j){
                    const auto col = entry[j];
                    PLOG_DEBUG << "col " << j << " " << col;
                    if (col == 0){
                        continue;
                    }
                    PLOG_DEBUG << "(" << i << "," << j << ")";
                    counts_(i,j) += 1;
                }
            }
        };
        void print_result(){
            std::cout << read_id_ << "\t" << num_hashes_ << "\t";
            for (auto i=0; i<counts_.rows(); i++){
                std::cout << +i << ":" << counts_(i,i) << "\t";
            }
            for (auto i=0; i<counts_.rows(); i++){
                for (auto j=0; j<i; j++) {
                    std::cout << +i << "x" << +j << ":" << counts_(i, j) << "\t";
                }
            }
            for (auto i=0; i<bits_.size(); i++){
                for (auto j=0; j<bits_[i].size(); j++) {
                    std::cout << +bits_[i][j];
                }
                std::cout << "\t";
            }
            std::cout << std::endl;
        };
    };

#endif // SIFTER_ENTRY_H
