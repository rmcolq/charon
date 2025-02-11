#ifndef CHARON_CLASSIFY_ARGUMENTS_H
#define CHARON_CLASSIFY_ARGUMENTS_H

#pragma once

#include <cstring>

/// Collection of all options of index subcommand.
struct ClassifyArguments {
    // IO options
    std::filesystem::path read_file;
    std::string db;
    uint8_t chunk_size { 100 };

    // Stats options
    float lo_hi_threshold {0.2};
    uint16_t num_reads_to_fit {5000};
    uint8_t confidence_threshold{0};
    uint8_t min_hits{3};
    float min_proportion_difference { 0.0 };
    std::string dist{"beta"};

    // Output options
    bool run_extract {false};
    std::string category_to_extract;
    std::filesystem::path extract_file;

    // General options
    std::string log_file {"charon.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };

    std::string to_string()
    {
        std::string ss;

        ss += "\n\nClassify Arguments:\n\n";
        ss += "\tread_file:\t\t" + read_file.string() + "\n";
        ss += "\tdb:\t\t\t" + db + "\n\n";

        ss += "\tchunk_size:\t\t" + std::to_string(chunk_size) + "\n";
        ss += "\tlo_hi_threshold:\t\t" + std::to_string(lo_hi_threshold) + "\n\n";
        ss += "\tnum_reads_to_fit:\t" + std::to_string(num_reads_to_fit) + "\n";
        ss += "\tconfidence_threshold:\t" + std::to_string(confidence_threshold) + "\n";
        ss += "\tmin_hits:\t\t" + std::to_string(min_hits) + "\n\n";
        ss += "\tmin_diff:\t\t" + std::to_string(min_proportion_difference) + "\n\n";
        ss += "\tdist:\t\t" + dist + "\n\n";

        ss += "\tcategory_to_extract:\t" + category_to_extract + "\n";
        ss += "\textract_file:\t\t" + extract_file.string() + "\n\n";

        ss += "\tlog_file:\t\t" + log_file + "\n";
        ss += "\tthreads:\t\t" + std::to_string(threads) + "\n";
        ss += "\tverbosity:\t\t" + std::to_string(verbosity) + "\n\n";

        return ss;
    }
};

#endif // CHARON_CLASSIFY_ARGUMENTS_H
