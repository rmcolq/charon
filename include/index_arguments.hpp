#ifndef CHARON_INDEX_ARGUMENTS_H
#define CHARON_INDEX_ARGUMENTS_H

#pragma once

#include <cstring>
#include <map>
#include <string>
#include <vector>

/// Collection of all options of index subcommand.
struct IndexArguments {
    // IO options
    std::string input_file;
    std::string prefix;
    std::string tmp_dir;

    // kmer/sketching
    uint8_t window_size{41};
    uint8_t kmer_size{19};

    // IBF options
    mutable size_t bits{std::numeric_limits<uint32_t>::max() - 2}; // Allow to change bits for each partition
    uint8_t num_hash{3};
    double max_fpr{0.01};

    // General options
    std::string log_file{"charon.log"};
    uint8_t threads{1};
    uint8_t verbosity{0};
    bool optimize{false};
    bool stats_only{false};  // count minimisers per file and report sizes, without building the IBF
    bool balance{false};       // merge small and split large bins so technical bins hold similar numbers of minimisers
    uint64_t bin_capacity{0};  // minimisers per technical bin when balancing (0 = total / technical_bins)
    uint32_t technical_bins{256};
    std::vector<std::string> category_fpr_args{};  // "category=fpr" strings from the command line
    std::map<std::string, double> category_fpr{};  // with --balance: lower target fpr for these categories (e.g. host)
    std::string profile{};          // "nanopore" or "illumina": preset w, k, fpr, balance and host fpr
    double profile_host_fpr{0};     // host fpr from the profile, applied to the "host"/"human" category once known

    std::string to_string() {
        std::string ss;

        ss += "\n\nIndex Arguments:\n\n";
        ss += "\tinput_file:\t\t" + input_file + "\n";
        ss += "\tprofile:\t\t" + (profile.empty() ? std::string{"none"} : profile) + "\n";
        ss += "\tprefix:\t\t\t" + prefix + "\n";
        ss += "\ttmp_dir:\t\t" + tmp_dir + "\n\n";

        ss += "\twindow_size:\t\t" + std::to_string(window_size) + "\n";
        ss += "\tkmer_size:\t\t" + std::to_string(kmer_size) + "\n\n";

        ss += "\tnum_hash:\t\t" + std::to_string(num_hash) + "\n";
        ss += "\tmax_fpr:\t\t" + std::to_string(max_fpr) + "\n\n";

        ss += "\toptimize:\t\t" + std::to_string(optimize) + "\n";
        ss += "\tstats_only:\t\t" + std::to_string(stats_only) + "\n";
        ss += "\tbalance:\t\t" + std::to_string(balance) + "\n";
        ss += "\tbin_capacity:\t\t" + std::to_string(bin_capacity) + "\n";
        ss += "\ttechnical_bins:\t\t" + std::to_string(technical_bins) + "\n";
        ss += "\tcategory_fpr:\t\t";
        for (const auto &[category, fpr]: category_fpr)
            ss += category + "=" + std::to_string(fpr) + " ";
        ss += "\n\n";

        ss += "\tlog_file:\t\t" + log_file + "\n";
        ss += "\tthreads:\t\t" + std::to_string(threads) + "\n";
        ss += "\tverbosity:\t\t" + std::to_string(verbosity) + "\n\n";

        return ss;
    }
};

#endif // CHARON_INDEX_ARGUMENTS_H
