#include <unordered_set>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cinttypes>
#include <ankerl/unordered_dense.h>

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

    index_subcommand->add_option("--fpr", opt->max_fpr, "Target false positive rate per bin.")
            ->type_name("FLOAT")
            ->check(CLI::Range(0.0001, 0.5))
            ->capture_default_str();

    index_subcommand->add_option("--num-hash", opt->num_hash, "Number of hash functions in the IBF.")
            ->type_name("INT")
            ->check(CLI::Range(1, 5))
            ->capture_default_str();

    index_subcommand->add_flag(
            "--balance", opt->balance,
            "Balance the IBF: pack small bins of a category together and split large bins into technical bins.")
            ->excludes("--optimize");

    index_subcommand->add_option("--bin-capacity", opt->bin_capacity,
                                 "Minimisers per technical bin with --balance (default: total / --technical-bins).")
            ->type_name("INT");

    index_subcommand->add_option("--technical-bins", opt->technical_bins,
                                 "Approximate number of technical bins with --balance when --bin-capacity is not set.")
            ->type_name("INT")
            ->check(CLI::Range(1, 4096))
            ->capture_default_str();

    index_subcommand->add_option("--merge-categories", opt->merge_categories,
                                 "With --balance: put all bins of these categories into a single scoring group "
                                 "(split into technical bins as needed).")
            ->type_name("STRING ...");

    index_subcommand->add_option("--category-fpr", opt->category_fpr_args,
                                 "With --balance: target fpr for the bins of a category, e.g. human=0.002. Its bins are "
                                 "split into more, less loaded technical bins (scoring takes the best bin per category, "
                                 "so the host category's false positives are amplified across its bins).")
            ->type_name("CATEGORY=FPR ...");

    index_subcommand->add_flag(
            "--stats-only", opt->stats_only,
            "Count distinct minimisers per file and write <prefix>.stats.tsv without building the index.");

    index_subcommand->add_flag(
            "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");


    // Set the function that will be called when this subcommand is issued.
    index_subcommand->callback([opt]() { index_main(*opt); });
}


InputSummary parse_input_file(const std::filesystem::path &input_file) {
    PLOG_INFO << "Parsing input file " << input_file;
    InputSummary summary;
    uint8_t next_bin = 0;
    ankerl::unordered_dense::set<std::string> categories;
    categories.reserve(10);  // Reserve space for typical number of categories

    std::ifstream input_ifstream{input_file};
    if (!input_ifstream.is_open()) {
        throw std::runtime_error{"Failed to open input file: " + input_file.string()};
    }

    std::string line;
    while (std::getline(input_ifstream, line)) {
        if (!line.empty()) {
            auto parts = split(line, "\t");
            if (parts.has_value() && parts->size() >= 2) {
            auto path = make_absolute(parts->at(0)).string();
            auto name = parts->at(1);
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
    summary.categories.reserve(categories.size());  // Pre-allocate
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
        try {
            const auto &fasta_file = pair.first;
            const auto &bin = pair.second;
            
            PLOG_DEBUG << "Adding file " << fasta_file;
            seqan3::sequence_file_input fin{fasta_file};
            
            auto record_count = 0;
            ankerl::unordered_dense::set<uint64_t> hashes;
            hashes.reserve(10000);  // Pre-allocate to reduce reallocations
            
            for (const auto &record: fin) {
                const auto mh = record.sequence() | hash_adaptor | std::views::common;
                hashes.insert(mh.begin(), mh.end());
                record_count++;
            }
            
            // Thread-safe updates to shared stats
#pragma omp critical
            {
                stats.records_per_bin[bin] += record_count;
                stats.num_files += 1;
                if (not opt.stats_only)
                    store_hashes(std::to_string(bin), hashes, opt.tmp_dir);
                stats.hashes_per_bin[bin] += hashes.size();
                PLOG_INFO << "Added file " << fasta_file << " with " << record_count << " records and " << hashes.size()
                          << " hashes to bin " << +bin;

                if (stats.hashes_per_bin[bin] > max_num_hashes) {
                    PLOG_WARNING
                                << "File " << fasta_file << " with " << hashes.size() << " will exceed max_fpr " << opt.max_fpr;
                }
            }
        } catch (const std::exception& e) {
            PLOG_ERROR << "Error processing file " << pair.first << ": " << e.what();
            // Re-throw to fail the build process
            throw;
        }
    }

    return stats;
}

ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>>
optimize_layout(const IndexArguments &opt, InputSummary &summary, InputStats &stats) {
    ankerl::unordered_dense::map<uint8_t, uint8_t> bin_to_bucket_map;
    ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>> bucket_to_bins_map;

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
    ankerl::unordered_dense::map<std::string, uint8_t> last_bin;

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

void write_stats(const IndexArguments &opt, const InputSummary &summary, const InputStats &stats) {
    const auto stats_file = opt.prefix.substr(0, opt.prefix.size() - 4) + ".stats.tsv";  // strip ".idx"
    std::ofstream out{stats_file};
    out << "file\tcategory\tbin\trecords\tminimisers\n";
    uint64_t total = 0;
    for (const auto &[path, bin]: summary.filepath_to_bin) {
        const auto hashes = stats.hashes_per_bin.contains(bin) ? stats.hashes_per_bin.at(bin) : 0;
        const auto records = stats.records_per_bin.contains(bin) ? stats.records_per_bin.at(bin) : 0;
        out << path << "\t" << summary.bin_to_category.at(bin) << "\t" << +bin << "\t" << records << "\t" << hashes << "\n";
        total += hashes;
    }
    const auto max_hashes = stats.max_num_hashes();
    // Bloom filter bits per element for the target fpr (uncapped, unlike bin_size_in_bits)
    const double bits_per_element = -static_cast<double>(opt.num_hash) /
                                    std::log(1 - std::exp(std::log(opt.max_fpr) / opt.num_hash));
    const double current_gb = bits_per_element * static_cast<double>(max_hashes) * summary.num_bins / 8e9;
    const double balanced_gb = bits_per_element * static_cast<double>(total) / 8e9;
    PLOG_INFO << "Wrote " << stats_file << ": " << +summary.num_bins << " bins, " << total
              << " minimisers in total, largest bin " << max_hashes;
    PLOG_INFO << "Uncompressed IBF size at fpr " << opt.max_fpr << " with " << +opt.num_hash
              << " hashes: current layout " << current_gb << " GB, perfectly balanced " << balanced_gb << " GB";
}

Index build_index(const IndexArguments &opt, const InputSummary &summary, InputStats &stats,
                  const ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>> &bucket_to_bins_map) {
    const auto max_num_hashes = stats.max_num_hashes();
    const auto num_bits = bin_size_in_bits(opt, max_num_hashes);
    PLOG_INFO << "Create new IBF with " << +summary.num_bins << " bins and " << +num_bits << " bits";
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{summary.num_bins},
                                         seqan3::bin_size{num_bits},
                                         seqan3::hash_function_count{opt.num_hash}};

#pragma omp parallel for
    for (uint8_t bucket = 0; bucket < summary.num_bins; ++bucket) {
        try {
            const auto &bins = bucket_to_bins_map.at(bucket);
            for (auto const &bin: bins) {
                const auto &hashes = load_hashes(std::to_string(bin), opt.tmp_dir);
#pragma omp critical
                {
                    for (auto &&value: hashes) {
                        ibf.emplace(value, seqan3::bin_index{bucket});
                    }
                }
                PLOG_DEBUG << "Added " << hashes.size() << " hashes to bin " << +bucket;
            }
            delete_hashes(bins, opt.tmp_dir);
        } catch (const std::exception& e) {
            PLOG_ERROR << "Error building bucket " << +bucket << ": " << e.what();
            throw;  // Re-throw to ensure proper cleanup
        }
    }


    return Index(opt, summary, stats, std::move(ibf));
}

std::pair<ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>>, Layout>
balance_layout(const IndexArguments &opt, InputSummary &summary, InputStats &stats) {
    uint64_t total = 0;
    for (const auto &[bin, num_hashes]: stats.hashes_per_bin)
        total += num_hashes;
    const uint64_t capacity = opt.bin_capacity > 0 ? opt.bin_capacity
                                                   : std::max<uint64_t>(1, (total + opt.technical_bins - 1) / opt.technical_bins);
    PLOG_INFO << "Balance index layout with capacity " << capacity << " minimisers per technical bin";

    // A category with a lower target fpr gets a smaller capacity per technical bin: technical bins all have
    // the height needed for `capacity` at opt.max_fpr, so loading them less lowers their fpr.
    const double bin_bits = static_cast<double>(capacity) * -static_cast<double>(opt.num_hash) /
                            std::log(1 - std::exp(std::log(opt.max_fpr) / opt.num_hash));
    auto capacity_for = [&](const std::string &category) -> uint64_t {
        const auto it = opt.category_fpr.find(category);
        if (it == opt.category_fpr.end() or it->second >= opt.max_fpr)
            return capacity;
        const double per_element = -static_cast<double>(opt.num_hash) /
                                   std::log(1 - std::exp(std::log(it->second) / opt.num_hash));
        return std::max<uint64_t>(1, static_cast<uint64_t>(bin_bits / per_element));
    };
    for (const auto &[category, fpr]: opt.category_fpr)
        PLOG_INFO << "Category " << category << " uses fpr " << fpr << ": capacity " << capacity_for(category)
                  << " minimisers per technical bin";

    // Pack bins smaller than the capacity together within each category (first fit decreasing);
    // larger bins keep a group of their own and are split into parts below.
    auto sorted_pairs = stats.bins_by_size();
    std::reverse(sorted_pairs.begin(), sorted_pairs.end());
    std::vector<uint64_t> group_load;
    std::vector<std::string> group_category;
    ankerl::unordered_dense::map<uint8_t, uint8_t> bin_to_group;
    ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>> group_to_bins;
    for (const auto &[bin, num_hashes]: sorted_pairs) {
        const auto &category = summary.bin_to_category.at(bin);
        const auto category_capacity = capacity_for(category);
        size_t group = group_load.size();
        const bool merge_all = std::find(opt.merge_categories.begin(), opt.merge_categories.end(), category) !=
                               opt.merge_categories.end();
        if (merge_all) {
            const auto it = std::find(group_category.begin(), group_category.end(), category);
            if (it != group_category.end())
                group = static_cast<size_t>(it - group_category.begin());
        } else if (num_hashes < category_capacity) {
            for (size_t g = 0; g < group_load.size(); ++g) {
                if (group_category[g] == category and group_load[g] + num_hashes <= category_capacity) {
                    group = g;
                    break;
                }
            }
        }
        if (group == group_load.size()) {
            group_load.push_back(0);
            group_category.push_back(category);
        }
        group_load[group] += num_hashes;
        bin_to_group[bin] = static_cast<uint8_t>(group);
        group_to_bins[static_cast<uint8_t>(group)].push_back(bin);
    }

    Layout layout;
    for (size_t group = 0; group < group_load.size(); ++group) {
        const auto group_capacity = capacity_for(group_category[group]);
        const auto parts = std::max<uint64_t>(1, (group_load[group] + group_capacity - 1) / group_capacity);
        layout.first_bin.push_back(static_cast<uint16_t>(layout.num_technical_bins));
        layout.num_parts.push_back(static_cast<uint16_t>(parts));
        layout.num_technical_bins += static_cast<uint32_t>(parts);
        PLOG_DEBUG << "Group " << group << " (" << group_category[group] << ") has " << group_to_bins[group].size()
                   << " bins, " << group_load[group] << " minimisers and " << parts << " parts";
    }
    layout.finalise();
    PLOG_INFO << "Balanced " << +summary.num_bins << " bins into " << group_load.size() << " scoring groups and "
              << layout.num_technical_bins << " technical bins";

    // Groups become the bins of the summary/stats used for scoring.
    InputSummary new_summary;
    InputStats new_stats;
    new_stats.num_files = stats.num_files;
    new_summary.num_bins = static_cast<uint8_t>(group_load.size());
    new_summary.categories = summary.categories;
    for (const auto &[filepath, bin]: summary.filepath_to_bin) {
        const auto group = bin_to_group.at(bin);
        new_summary.filepath_to_bin.emplace_back(filepath, group);
        new_summary.bin_to_category[group] = summary.bin_to_category.at(bin);
        new_stats.records_per_bin[group] += stats.records_per_bin.at(bin);
    }
    for (size_t group = 0; group < group_load.size(); ++group)
        new_stats.hashes_per_bin[static_cast<uint8_t>(group)] = group_load[group];
    summary = new_summary;
    stats = new_stats;
    return {group_to_bins, layout};
}

Index build_balanced_index(const IndexArguments &opt, const InputSummary &summary, InputStats &stats,
                           const ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>> &group_to_bins_map,
                           const Layout &layout) {
    // Exact load of every technical bin (merged groups may share minimisers, so this is an upper bound for them)
    std::vector<uint64_t> technical_bin_load(layout.num_technical_bins, 0);
#pragma omp parallel for num_threads(opt.threads) schedule(dynamic)
    for (size_t group = 0; group < layout.num_groups(); ++group) {
        for (auto const &bin: group_to_bins_map.at(static_cast<uint8_t>(group))) {
            const auto hashes = load_hashes(std::to_string(bin), opt.tmp_dir);
            for (const auto value: hashes)
                technical_bin_load[layout.technical_bin(group, value)] += 1;  // each group owns its own technical bins
        }
    }
    const auto max_load = *std::max_element(technical_bin_load.begin(), technical_bin_load.end());
    const auto num_bits = bin_size_in_bits(opt, max_load);
    PLOG_INFO << "Create new IBF with " << layout.num_technical_bins << " technical bins of " << num_bits
              << " bits (largest technical bin holds " << max_load << " minimisers)";
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{layout.num_technical_bins},
                                         seqan3::bin_size{num_bits},
                                         seqan3::hash_function_count{opt.num_hash}};

    for (size_t group = 0; group < layout.num_groups(); ++group) {
        const auto &bins = group_to_bins_map.at(static_cast<uint8_t>(group));
        for (auto const &bin: bins) {
            const auto hashes = load_hashes(std::to_string(bin), opt.tmp_dir);
            for (const auto value: hashes)
                ibf.emplace(value, seqan3::bin_index{layout.technical_bin(group, value)});
            PLOG_DEBUG << "Added " << hashes.size() << " hashes from bin " << +bin << " to group " << group;
        }
        delete_hashes(bins, opt.tmp_dir);
    }
    return Index(opt, summary, stats, std::move(ibf), layout);
}

int index_main(IndexArguments &opt) {
    auto log_level = plog::info;
    if (opt.verbosity == 1) {
        log_level = plog::debug;
    } else if (opt.verbosity > 1) {
        log_level = plog::verbose;
    }
    plog::init(log_level, opt.log_file.c_str(), 1000000, 5);

    for (const auto &arg: opt.category_fpr_args) {
        const auto parts = split(arg, "=");
        if (not parts.has_value() or parts->size() != 2)
            throw std::invalid_argument("--category-fpr expects CATEGORY=FPR, got " + arg);
        opt.category_fpr[parts->at(0)] = std::stod(parts->at(1));
    }

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
    if (not opt.stats_only)
        std::filesystem::create_directory(opt.tmp_dir);

    auto args = opt.to_string();
    LOG_INFO << "Running charon index\n\nCharon version: " << SOFTWARE_VERSION << "\n" << args;


    auto summary = parse_input_file(opt.input_file);
    auto stats = count_and_store_hashes(opt, summary);
    if (opt.stats_only) {
        write_stats(opt, summary, stats);
        return 0;
    }
    if (opt.balance) {
        auto [group_to_bins_map, layout] = balance_layout(opt, summary, stats);
        store_index(opt.prefix, build_balanced_index(opt, summary, stats, group_to_bins_map, layout));
        return 0;
    }

    auto bucket_to_bins_map = optimize_layout(opt, summary, stats);
    auto index = build_index(opt, summary, stats, bucket_to_bins_map);

    store_index(opt.prefix, std::move(index));

    return 0;
}
