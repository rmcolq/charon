#ifndef CHARON_DEHOST_ARGUMENTS_H
#define CHARON_DEHOST_ARGUMENTS_H

#pragma once

#include <cstring>
#include <ankerl/unordered_dense.h>

/// Collection of all options of dehost subcommand.
struct DehostArguments {
    // IO options
    std::filesystem::path read_file;
    std::filesystem::path read_file2;
    bool is_paired{false};
    std::string db;

    // Output options
    bool run_extract{false};
    std::string category_to_extract;
    std::string prefix;
    ankerl::unordered_dense::map<uint8_t, std::vector<std::filesystem::path>> extract_category_to_file;

    uint16_t chunk_size{0};  // 0 = auto-detect
    uint16_t max_memory_gb{0};  // 0 = auto-detect
    // Set from the memory manager at run time (not command line options); 0 = unlimited.
    uint64_t chunk_base_budget{0};  // maximum bases (summed read length) in one chunk
    uint64_t cache_byte_budget{0};  // maximum bytes of reads cached while the model is fitted

    // Stats options
    float lo_hi_threshold{0.15};
    uint16_t num_reads_to_fit{5000};
    std::string dist{"kde"};

    // thresholds for filtering; defaults are for single-end (nanopore) reads
    float min_quality{15.0};
    uint32_t min_length{140};
    float min_compression{0};
    uint8_t confidence_threshold{7};
    float confidence_probability_threshold{0};
    float host_unique_prop_lo_threshold{0.04};
    float min_proportion_difference{0.07};
    float min_prob_difference{0};

    // Defaults for paired (Illumina) reads, used for any of these options not given explicitly.
    // Chosen so that no human pair left called microbial in a real-read Illumina mixture is human by BLAST.
    // min_length applies to each mate, so it must be below the mate length.
    static constexpr uint32_t paired_min_length{70};
    static constexpr float paired_min_compression{0.3};
    static constexpr float paired_host_unique_prop_lo_threshold{0.05};
    static constexpr float paired_min_proportion_difference{0.28};

    // true when the option was given explicitly, so a paired default does not replace it
    bool min_length_set{false};
    bool min_compression_set{false};
    bool host_unique_prop_lo_threshold_set{false};
    bool min_proportion_difference_set{false};

    void apply_paired_defaults() {
        if (not min_length_set)
            min_length = paired_min_length;
        if (not min_compression_set)
            min_compression = paired_min_compression;
        if (not host_unique_prop_lo_threshold_set)
            host_unique_prop_lo_threshold = paired_host_unique_prop_lo_threshold;
        if (not min_proportion_difference_set)
            min_proportion_difference = paired_min_proportion_difference;
    }


    // General options
    std::string log_file{"charon.log"};
    uint8_t threads{1};
    uint8_t verbosity{0};

    std::string to_string() {
        std::string ss;

        ss += "\n\nDehost Arguments:\n\n";
        ss += "\tread_file:\t\t\t" + read_file.string() + "\n";
        ss += "\tread_file2:\t\t\t" + read_file2.string() + "\n";
        ss += "\tdb:\t\t\t\t" + db + "\n\n";

        ss += "\tcategory_to_extract:\t\t" + category_to_extract + "\n";
        ss += "\tprefix:\t\t\t\t" + prefix + "\n\n";

        ss += "\tchunk_size:\t\t\t" + std::to_string(chunk_size) + "\n";
        ss += "\tmax_memory_gb:\t\t\t" + (max_memory_gb > 0 ? std::to_string(max_memory_gb) : "auto") + "\n\n";
        ss += "\tlo_hi_threshold:\t\t" + std::to_string(lo_hi_threshold) + "\n";
        ss += "\tnum_reads_to_fit:\t\t" + std::to_string(num_reads_to_fit) + "\n";
        ss += "\tdist:\t\t\t\t" + dist + "\n\n";

        ss += "\tmin_length:\t\t\t" + std::to_string(min_length) + "\n";
        ss += "\tmin_quality:\t\t\t" + std::to_string(min_quality) + "\n";
        ss += "\tmin_compression:\t\t" + std::to_string(min_compression) + "\n";
        ss += "\tconfidence_threshold:\t\t" + std::to_string(confidence_threshold) + "\n";
        ss += "\tconfidence_probability_threshold:\t\t" + std::to_string(confidence_probability_threshold) + "\n";
        ss += "\thost_unique_prop_lo_threshold:\t" + std::to_string(host_unique_prop_lo_threshold) + "\n";
        ss += "\tmin_proportion_difference:\t" + std::to_string(min_proportion_difference) + "\n";
        ss += "\tmin_prob_difference:\t\t" + std::to_string(min_prob_difference) + "\n\n";


        ss += "\tlog_file:\t\t\t" + log_file + "\n";
        ss += "\tthreads:\t\t\t" + std::to_string(threads) + "\n";
        ss += "\tverbosity:\t\t\t" + std::to_string(verbosity) + "\n\n";

        return ss;
    }
};

#endif // CHARON_DEHOST_ARGUMENTS_H
