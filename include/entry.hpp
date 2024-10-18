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
        uint32_t num_hashes_{0};
        Counts<uint32_t> counts_;

    public:
        ReadEntry() = default;
        ReadEntry(ReadEntry const &) = default;
        ReadEntry(ReadEntry &&) = default;
        ReadEntry & operator=(ReadEntry const &) = default;
        ReadEntry & operator=(ReadEntry &&) = default;
        ~ReadEntry() = default;

        ReadEntry(const std::string read_id, const uint8_t num_bins):
            read_id_(read_id)
            {counts_.set_size(num_bins);}

        void update_entry(const auto & entry){
            num_hashes_ += 1;
            for (const auto i: entry)
            {
                if (i == 0){
                    continue;
                }
                for (const auto j: entry){
                    if (j == 0){
                        continue;
                    } else if ( i < j ){
                        continue;
                    }
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
            std::cout << std::endl;
        };
    };

#endif // SIFTER_ENTRY_H
