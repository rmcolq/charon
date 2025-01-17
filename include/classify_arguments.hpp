#ifndef CHARON_CLASSIFY_ARGUMENTS_H
#define CHARON_CLASSIFY_ARGUMENTS_H

#pragma once

#include <cstring>

/// Collection of all options of index subcommand.
struct ClassifyArguments {
    // IO options
    std::string read_file;
    std::string db;
    uint8_t chunk_size { 100 };

    // Stats options
    float lo_hi_threshold {0.2};
    uint16_t num_reads_to_fit {1000};

    // Output options
    std::string category_to_extract;
    std::string extract_file;

    // General options
    std::string log_file {"charon.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };
};

#endif // CHARON_CLASSIFY_ARGUMENTS_H
