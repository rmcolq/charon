#ifndef CHARON_INDEX_ARGUMENTS_H
#define CHARON_INDEX_ARGUMENTS_H

#pragma once

#include <cstring>

/// Collection of all options of index subcommand.
struct IndexArguments {
    // IO options
    std::string input_file;
    std::string prefix;
    std::string tmp_dir;

    // kmer/sketching
    uint8_t window_size { 41 };
    uint8_t kmer_size { 19 };

    // IBF options
    mutable size_t bits {std::numeric_limits<uint32_t>::max()}; // Allow to change bits for each partition
    uint8_t num_hash {3};
    double max_fpr {0.02};

    // General options
    std::string log_file {"charon.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };
    bool optimize { false };

    std::string to_string()
    {
        std::string ss;

        ss += "\n\nIndex Arguments:\n\n";
        ss += "\tinput_file:\t\t" + input_file + "\n";
        ss += "\tprefix:\t\t\t" + prefix + "\n";
        ss += "\ttmp_dir:\t\t" + tmp_dir + "\n\n";

        ss += "\twindow_size:\t\t" + std::to_string(window_size) + "\n";
        ss += "\tkmer_size:\t\t" + std::to_string(kmer_size) + "\n\n";

        ss += "\tnum_hash:\t\t" + std::to_string(num_hash) + "\n";
        ss += "\tmax_fpr:\t\t" + std::to_string(max_fpr) + "\n\n";

        ss += "\toptimize:\t\t" + std::to_string(optimize) + "\n\n";

        ss += "\tlog_file:\t\t" + log_file + "\n";
        ss += "\tthreads:\t\t" + std::to_string(threads) + "\n";
        ss += "\tverbosity:\t\t" + std::to_string(verbosity) + "\n\n";

        return ss;
    }
};

#endif // CHARON_INDEX_ARGUMENTS_H
