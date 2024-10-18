#ifndef SIFTER_INDEX_MAIN_H
#define SIFTER_INDEX_MAIN_H

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"

class Index;

/// Collection of all options of index subcommand.
struct IndexArguments {
    // IO options
    std::string input_file;
    std::string prefix;

    // kmer/sketching
    uint8_t window_size { 31 };
    uint8_t kmer_size { 15 };

    // IBF options
    uint8_t bins {2};
    mutable uint64_t bits {4096}; // Allow to change bits for each partition
    uint8_t hash {2};
    double fpr {0.05};

    // General options
    std::string log_file {"sifter.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };
};

struct InputFileMap {
    std::unordered_map<std::string, uint8_t> filepath_to_bin;
    std::unordered_map<uint8_t, std::string> bin_to_name;
};

struct InputSummary
{
    uint8_t num_bins                                        = 0;
    uint32_t num_files                                      = 0;
    std::unordered_map<uint8_t, uint64_t> records_per_bin   = {};
    std::unordered_map<uint8_t, uint64_t> hashes_per_bin    = {};
};

void setup_index_subcommand(CLI::App& app);

InputFileMap parse_input_file(const std::filesystem::path& input_file);

InputSummary summarise_input(const InputFileMap& input, const InputSummary& summary, const IndexArguments& opt);

Index build_index(const InputFileMap& input, InputSummary& summary, const IndexArguments& opt);

int index_main(IndexArguments & opt);


#endif // SIFTER_INDEX_MAIN_H
