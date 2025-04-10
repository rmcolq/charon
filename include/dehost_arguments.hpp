#ifndef CHARON_DEHOST_ARGUMENTS_H
#define CHARON_DEHOST_ARGUMENTS_H

#pragma once

#include <cstring>

/// Collection of all options of dehost subcommand.
struct DehostArguments {
    // IO options
    std::filesystem::path read_file;
    std::filesystem::path read_file2;
    bool is_paired { false };
    std::string db;
    uint8_t chunk_size { 100 };


    // Stats options
    float lo_hi_threshold {0.15};
    uint16_t num_reads_to_fit {5000};
    std::string dist{"beta"};

    // thresholds for filtering
    float min_quality { 10.0 };
    uint32_t min_length { 500 };
    uint8_t confidence_threshold{2};
    float host_unique_prop_hi_threshold{ 0.1 };
    float host_unique_prop_lo_threshold{ 0.05 };
    float min_proportion_difference { 0.04 };

    // Output options
    bool run_extract {false};
    std::string category_to_extract;
    std::filesystem::path extract_file;
    std::filesystem::path extract_file2;

    // General options
    std::string log_file {"charon.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };

    std::string to_string()
    {
        std::string ss;

        ss += "\n\nDehost Arguments:\n\n";
        ss += "\tread_file:\t\t" + read_file.string() + "\n";
        ss += "\tdb:\t\t\t" + db + "\n\n";

        ss += "\tchunk_size:\t\t" + std::to_string(chunk_size) + "\n\n";

        ss += "\tlo_hi_threshold:\t\t" + std::to_string(lo_hi_threshold) + "\n";
        ss += "\tnum_reads_to_fit:\t" + std::to_string(num_reads_to_fit) + "\n";
        ss += "\tdist:\t\t" + dist + "\n\n";

        ss += "\tmin_length:\t\t" + std::to_string(min_length) + "\n";
        ss += "\tmin_quality:\t\t" + std::to_string(min_quality) + "\n";
        ss += "\tconfidence_threshold:\t" + std::to_string(confidence_threshold) + "\n";
        ss += "\thost_unique_prop_lo_hreshold:\t" + std::to_string(host_unique_prop_lo_threshold) + "\n";
        ss += "\thost_unique_prop_hi_threshold:\t" + std::to_string(host_unique_prop_hi_threshold) + "\n";
        ss += "\tmin_proportion_difference:\t\t" + std::to_string(min_proportion_difference) + "\n\n";
        
        ss += "\tcategory_to_extract:\t" + category_to_extract + "\n";
        ss += "\textract_file:\t\t" + extract_file.string() + "\n\n";

        ss += "\tlog_file:\t\t" + log_file + "\n";
        ss += "\tthreads:\t\t" + std::to_string(threads) + "\n";
        ss += "\tverbosity:\t\t" + std::to_string(verbosity) + "\n\n";

        return ss;
    }
};

#endif // CHARON_DEHOST_ARGUMENTS_H
