#ifndef SIFTER_INDEX_MAIN_H
#define SIFTER_INDEX_MAIN_H

#pragma once

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"

class Index;
class InputStats;
class InputSummary;

/// Collection of all options of index subcommand.
struct IndexArguments {
    // IO options
    std::string input_file;
    std::string prefix;

    // kmer/sketching
    uint8_t window_size { 41 };
    uint8_t kmer_size { 19 };

    // IBF options
    uint8_t bins {2};
    mutable uint64_t bits {409600000}; // Allow to change bits for each partition
    uint8_t hash {2};
    double max_fpr {0.05};

    // General options
    std::string log_file {"sifter.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };
};

void setup_index_subcommand(CLI::App& app);

InputSummary parse_input_file(const std::filesystem::path& input_file);

InputStats estimate_index_size(const InputSummary& summary, const IndexArguments& opt);

InputStats summarise_input(const InputSummary& summary, const IndexArguments& opt);

Index build_index(const InputSummary& summary, InputStats& stats, const IndexArguments& opt);

int index_main(IndexArguments & opt);


#endif // SIFTER_INDEX_MAIN_H
