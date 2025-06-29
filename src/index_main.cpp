#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "index_main.hpp"
#include "utils.hpp"
#include "index.hpp"
#include "store_index.hpp"
#include "input_summary.hpp"
#include "version.h"

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>


void setup_index_subcommand(CLI::App &app) {
    auto opt = std::make_shared<IndexArguments>();
    auto *index_subcommand = app.add_subcommand(
            "index", "Build an index (IBF) for a number of references split into a small number of bins.");

    index_subcommand->add_option("<input>", opt->input_file,
                                 "Tab separated file with columns for filename and category")
            ->required()
            ->transform(make_absolute)
            ->check(CLI::ExistingFile.description(""))
            ->type_name("FILE");

    index_subcommand
            ->add_option(
                    "-w", opt->window_size, "Window size for (w,k,s)-minimers (must be <=k).")
            ->type_name("INT")
            ->capture_default_str();

    index_subcommand->add_option("-k", opt->kmer_size, "K-mer size for (w,k,s)-minimers.")
            ->type_name("INT")
            ->capture_default_str();

    index_subcommand
            ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use.")
            ->type_name("INT")
            ->capture_default_str();

    index_subcommand->add_option("-p,--prefix", opt->prefix, "Prefix for the output index.")
            ->type_name("FILE")
            ->check(CLI::NonexistentPath.description(""))
            ->default_str("<prefix>");

    index_subcommand->add_option("--temp", opt->tmp_dir, "Temporary directory for index construction files.")
            ->type_name("DIR")
            ->default_str("<dir>");

    index_subcommand->add_option("--log", opt->log_file, "File for log")
            ->transform(make_absolute)
            ->type_name("FILE");

    index_subcommand->add_flag(
            "--optimize", opt->optimize, "Compress the number of bins for improved classification run time");

    index_subcommand->add_flag(
            "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");


    // Set the function that will be called when this subcommand is issued.
    index_subcommand->callback([opt]() { index_main(*opt); });
}


InputSummary parse_input_file(const std::filesystem::path &input_file) {
    PLOG_INFO << "Parsing input file " << input_file;
    InputSummary summary;
    uint8_t next_bin = 0;
    std::unordered_set<std::string> categories;

    std::ifstream input_ifstream;
    input_ifstream.open(input_file);
    if (!input_ifstream.is_open()) {
        PLOG_ERROR << "Error opening file " << input_file;
        exit(1);
    }

    std::string line;
    while (std::getline(input_ifstream, line)) {
        if (!line.empty()) {
            auto parts = split(line, "\t");
            if (parts.size() >= 2) {
                auto path = make_absolute(parts[0]).string();
                auto name = parts[1];
                summary.bin_to_category[next_bin] = name;
                categories.insert(name);
                summary.filepath_to_bin.emplace_back(std::make_pair(path, next_bin));
                if (next_bin == std::numeric_limits<uint8_t>::max()) {
                    PLOG_WARNING << "User has reached the maximum number of files which is " << +next_bin
                                 << " - ignoring any additional lines!";
                    break;
                }
                next_bin++;
            }
        }
    }
    input_ifstream.close();

    summary.num_bins = next_bin;
    summary.categories.insert(summary.categories.end(), categories.begin(), categories.end());

    PLOG_INFO << "Found " << summary.filepath_to_bin.size() << " files corresponding to " << +summary.num_categories()
              << " categories";

    return summary;
}

InputStats count_and_store_hashes(const IndexArguments &opt, const InputSummary &summary) {
    PLOG_INFO << "Extracting hashes from files";
    const auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{opt.kmer_size}},
                                                            seqan3::window_size{opt.window_size});
    PLOG_DEBUG << "Defined hash adaptor";
    InputStats stats;
    PLOG_DEBUG << "Defined stats";

    const auto max_num_hashes = max_num_hashes_for_fpr(opt);
    PLOG_INFO << "Maximum hashes permitted per bin for fpr rate " << opt.max_fpr << " is " << max_num_hashes;

#pragma omp parallel for num_threads(opt.threads)
    for (const auto pair: summary.filepath_to_bin) {
        const auto &fasta_file = pair.first;
        const auto &bin = pair.second;
        stats.records_per_bin[bin] += 0;

        PLOG_DEBUG << "Adding file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
#pragma omp critical
        stats.num_files += 1;

        auto record_count = 0;
        std::unordered_set<uint64_t> hashes;
        for (const auto &record: fin) {
            stats.records_per_bin[bin] += 1;
            const auto mh = record.sequence() | hash_adaptor | std::views::common;
            hashes.insert(mh.begin(), mh.end());
            record_count++;
        }
        store_hashes(std::to_string(bin), hashes, opt.tmp_dir);
        stats.hashes_per_bin[bin] += hashes.size();
        PLOG_INFO << "Added file " << fasta_file << " with " << record_count << " records and " << hashes.size()
                  << " hashes to bin " << +bin;

        if (stats.hashes_per_bin[bin] > max_num_hashes) {
            PLOG_WARNING
                        << "File " << fasta_file << " with " << hashes.size() << " will exceed max_fpr " << opt.max_fpr;
        }
    }

    return stats;
}

std::unordered_map<uint8_t, std::vector<uint8_t>>
optimize_layout(const IndexArguments &opt, InputSummary &summary, InputStats &stats) {
    std::unordered_map<uint8_t, uint8_t> bin_to_bucket_map;
    std::unordered_map<uint8_t, std::vector<uint8_t>> bucket_to_bins_map;

    if (stats.hashes_per_bin.size() == 0) {
        return bucket_to_bins_map;
    } else if (not opt.optimize) {
        for (auto bin = 0; bin < summary.num_bins; ++bin) {
            bucket_to_bins_map[bin].push_back(bin);
        }
        return bucket_to_bins_map;
    }

    PLOG_INFO << "Optimize index bin layout";
    auto new_summary = InputSummary();
    auto new_stats = InputStats();
    new_stats.num_files = stats.num_files;

    auto sorted_pairs = stats.bins_by_size();
    auto max_num_hashes = sorted_pairs.back().second / 2;
    PLOG_DEBUG << "Max hashes found for bin " << +sorted_pairs.back().first << " : " << max_num_hashes;

    uint8_t next_bin = 0;
    std::unordered_map<std::string, uint8_t> last_bin;

    // update stats
    PLOG_INFO << "Reassign bins";
    for (const auto &pair: sorted_pairs) {
        const auto &bin = pair.first;
        PLOG_DEBUG << "Reassign bin " << +bin;
        assert(summary.bin_to_category.find(bin) != summary.bin_to_category.end());
        const auto &category = summary.bin_to_category.at(bin);
        const auto &num_hashes = pair.second;

        auto assigned_bucket = next_bin;
        if (last_bin.find(category) != last_bin.end()) {
            const auto &last_bin_in_category = last_bin[category];
            const auto &num_hashes_in_last_bin = new_stats.hashes_per_bin[last_bin_in_category];
            if (num_hashes_in_last_bin + num_hashes < max_num_hashes) {
                assigned_bucket = last_bin_in_category;
            } else {
                next_bin++;
            }
        } else {
            next_bin++;
        }
        last_bin[category] = assigned_bucket;
        bin_to_bucket_map[bin] = assigned_bucket;
        bucket_to_bins_map[assigned_bucket].push_back(bin);
        new_stats.hashes_per_bin[assigned_bucket] += num_hashes;
        PLOG_INFO << "Bin " << +bin << " assigned to bucket " << +assigned_bucket << " which now has "
                  << new_stats.hashes_per_bin[assigned_bucket] << " hashes";
        assert(stats.records_per_bin.find(bin) != stats.records_per_bin.end());
        new_stats.records_per_bin[assigned_bucket] += stats.records_per_bin.at(bin);
    }
    stats = new_stats;

    // update summary
    PLOG_INFO << "Update summary";
    new_summary.num_bins = next_bin;
    new_summary.categories = summary.categories;
    for (const auto &pair: summary.filepath_to_bin) {
        const auto &filepath = pair.first;
        const auto &bin = pair.second;
        PLOG_DEBUG << "Update summary for bin " << +bin;
        const auto &bucket = bin_to_bucket_map[bin];
        new_summary.filepath_to_bin.emplace_back(std::make_pair(filepath, bucket));
        assert(summary.bin_to_category.find(bin) != summary.bin_to_category.end());
        new_summary.bin_to_category[bucket] = summary.bin_to_category.at(bin);
    }
    summary = new_summary;

    return bucket_to_bins_map;
}

Index build_index(const IndexArguments &opt, const InputSummary &summary, InputStats &stats,
                  const std::unordered_map<uint8_t, std::vector<uint8_t>> &bucket_to_bins_map) {
    const auto max_num_hashes = stats.max_num_hashes();
    const auto num_bits = bin_size_in_bits(opt, max_num_hashes);
    PLOG_INFO << "Create new IBF with " << +summary.num_bins << " bins and " << +num_bits << " bits";
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{summary.num_bins},
                                         seqan3::bin_size{num_bits},
                                         seqan3::hash_function_count{opt.num_hash}};

#pragma omp parallel for
    for (uint8_t bucket = 0; bucket < summary.num_bins; ++bucket) {
        const auto &bins = bucket_to_bins_map.at(bucket);
        for (auto const &bin: bins) {
            const auto &hashes = load_hashes(std::to_string(bin), opt.tmp_dir);
#pragma omp critical
            for (auto &&value: hashes) {
                ibf.emplace(value, seqan3::bin_index{bucket});
            }
            PLOG_DEBUG << "Added " << hashes.size() << " hashes to bin " << +bucket;
        }
        delete_hashes(bins, opt.tmp_dir);
    }


    return Index(opt, summary, stats, ibf);
}

int index_main(IndexArguments &opt) {
    auto log_level = plog::info;
    if (opt.verbosity == 1) {
        log_level = plog::debug;
    } else if (opt.verbosity > 1) {
        log_level = plog::verbose;
    }
    plog::init(log_level, opt.log_file.c_str(), 1000000, 5);

    if (opt.window_size < opt.kmer_size) {
        throw std::logic_error("W must be greater than K");
    }
    if (opt.window_size <= 0) {
        throw std::logic_error("W must be a positive integer");
    }
    if (opt.kmer_size <= 0) {
        throw std::logic_error("K must be a positive integer");
    }
    if (opt.prefix != "") {
        opt.prefix += ".idx";
    } else {
        opt.prefix = opt.input_file + ".idx";
    }
    if (opt.tmp_dir == "") {
        opt.tmp_dir = opt.input_file + ".tmp_idx";
    }
    std::filesystem::create_directory(opt.tmp_dir);

    auto args = opt.to_string();
    LOG_INFO << "Running charon index\n\nCharon version: " << SOFTWARE_VERSION << "\n" << args;


    auto summary = parse_input_file(opt.input_file);
    auto stats = count_and_store_hashes(opt, summary);
    auto bucket_to_bins_map = optimize_layout(opt, summary, stats);
    auto index = build_index(opt, summary, stats, bucket_to_bins_map);

    store_index(opt.prefix, std::move(index));

    return 0;
}
