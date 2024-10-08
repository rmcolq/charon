#ifndef SIFTER_INDEX_MAIN_H
#define SIFTER_INDEX_MAIN_H

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"

/// Collection of all options of index subcommand.
struct IndexOptions {
    std::string input_file;
    uint8_t window_size { 14 };
    uint8_t kmer_size { 15 };
    uint8_t threads { 1 };
    std::string prefix;
    uint8_t verbosity { 0 };
};

void setup_index_subcommand(CLI::App& app);
int sifter_index(IndexOptions const& opt);

#endif // SIFTER_INDEX_MAIN_H
