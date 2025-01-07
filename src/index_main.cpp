#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>

#include "index_main.hpp"
#include "index.hpp"
#include "store_index.hpp"
#include "utils.hpp"
#include "input_summary.hpp"

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>


void setup_index_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<IndexArguments>();
    auto* index_subcommand = app.add_subcommand(
        "index", "Build an index (IBF) for a number of references split into a small number of bins.");

    index_subcommand->add_option("<input>", opt->input_file, "Tab separated file with columns for filename and category")
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
        ->default_str("<input>");

    index_subcommand->add_option("--log", opt->log_file, "File for log")
            ->transform(make_absolute)
            ->type_name("FILE");

    index_subcommand->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    index_subcommand->add_flag(
            "--optimize", opt->optimize, "Compress the number of bins for improved classification run time");

    // Set the function that will be called when this subcommand is issued.
    index_subcommand->callback([opt]() { index_main(*opt); });
}

inline constexpr size_t bin_size_in_bits(const IndexArguments & opt, const uint64_t & num_elements)
{
    assert(opt.num_hash > 0);
    assert(opt.max_fpr > 0.0);
    assert(opt.max_fpr < 1.0);

    double const numerator{-static_cast<double>(num_elements * opt.num_hash)};
    double const denominator{std::log(1 - std::exp(std::log(opt.max_fpr) / opt.num_hash))};
    double const result{std::ceil(numerator / denominator)};
    return result;
}

InputSummary parse_input_file(const std::filesystem::path& input_file){
    std::cout << "Parsing input file " << std::endl;
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
    while (std::getline(input_ifstream, line))
    {
        if (!line.empty()) {
            auto parts = split(line, "\t");
            if (parts.size() >= 2) {
                auto path = make_absolute(parts[0]).string();
                auto name = parts[1];
                summary.bin_to_category[next_bin] = name;
                categories.insert(name);
                summary.filepath_to_bin.emplace_back(std::make_pair(path, next_bin));
                next_bin++;
            }
        }
    }
    input_ifstream.close();

    summary.num_bins = next_bin;
    summary.categories.insert(summary.categories.end(), categories.begin(), categories.end());

    PLOG_INFO << "Found " << summary.filepath_to_bin.size() << " files corresponding to " << +summary.num_categories() << " categories over " << +summary.num_bins << " bins";
    std::cout << "Found " << summary.filepath_to_bin.size() << " files corresponding to " << +summary.num_categories() << " categories over " << +summary.num_bins << " bins" << std::endl;

    return summary;
}

InputStats estimate_index_size(const InputSummary& summary, const IndexArguments& opt){
    PLOG_INFO << "Estimate size from files";
    InputStats stats;

    std::unordered_map<uint8_t, uint64_t> lengths;

    for (const auto& pair : summary.filepath_to_bin) {
        const auto& fasta_file = pair.first;
        const auto& bin = pair.second;

        PLOG_INFO << "Checking file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
        stats.num_files += 1;

        auto length = 0;
        for (auto & record : fin)
        {
            length += record.sequence().size();
            stats.records_per_bin[bin]++;
        }
        lengths[bin] += length;
    }

    // estimate max bin size
    uint64_t max_count = 0;
    for(auto kv : lengths) {
        //stats.num_bins += 1;
        stats.hashes_per_bin[kv.first] = static_cast<uint64_t>(round((kv.second * 2)/(opt.window_size + 1)));
        max_count = std::max(max_count, static_cast<uint64_t>(round(1.1*stats.hashes_per_bin[kv.first])));
    }
    opt.bits = std::ceil( ( max_count * std::log( opt.max_fpr ) ) / std::log( 1.0 / std::pow( 2, std::log( 2 ) ) ) );

    //PLOG_DEBUG << "num_bins: " << +stats.num_bins << ", num_files: " << stats.num_files;
    for (const auto& item: stats.records_per_bin){
        PLOG_DEBUG << "Bin " << +item.first << " has " << item.second << " records";
    }
    for (const auto& item: stats.hashes_per_bin){
        PLOG_DEBUG << "Bin " << +item.first << " has " << item.second << " hashes";
    }
    PLOG_DEBUG << "using: " << opt.bits << " to achieve a max_fpr of " << opt.max_fpr;

    return stats;
}

InputStats summarise_input(const InputSummary& summary, const IndexArguments& opt){
    PLOG_INFO << "Estimate size and stats information from files";
    InputStats stats;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{opt.kmer_size}}, seqan3::window_size{opt.window_size});
    std::unordered_map<uint8_t, std::unordered_set<uint64_t>> hashes;

    for (const auto& pair : summary.filepath_to_bin) {
        const auto& fasta_file = pair.first;
        const auto& bin = pair.second;

        bool bin_seen_before = hashes.contains(bin);
        if (bin_seen_before){
            continue;
        }

        PLOG_INFO << "Checking file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
        //stats.num_files += 1;

        auto record_count = 0;
        for (auto & record : fin)
        {
            //stats.records_per_bin[bin] += 1;
            const auto mh = record.sequence() | hash_adaptor | std::views::common;
            hashes[bin].insert( mh.begin(), mh.end() );
            record_count++;
        }
        //PLOG_INFO << "Added file " << fasta_file << " with " << record_count << " records to bin " << bin << std::endl;
    }
    for (const auto& pair : hashes) {
        const auto& bin = pair.first;
        const auto& num_hashes = pair.second.size();
        stats.hashes_per_bin[bin] = num_hashes;
        //stats.num_bins += 1;
    }

    // estimate max bin size
    uint64_t max_count = 0;
    for(auto kv : stats.hashes_per_bin) {
        max_count = std::max(max_count, static_cast<uint64_t>(round(kv.second*1.1)));
    }
    opt.bits = max_count;

    //std::cout << "num_bins: " << +stats.num_bins << ", num_files: " << stats.num_files << std::endl;
    for (const auto& item: stats.records_per_bin){
        std::cout << "Bin " << +item.first << " has " << item.second << " records" << std::endl;
    }
    for (const auto& item: stats.hashes_per_bin){
        std::cout << "Bin " << +item.first << " has " << item.second << " hashes" << std::endl;
    }

    return stats;
}

Index build_index(const InputSummary& summary, InputStats& stats, const IndexArguments& opt)
{
    PLOG_INFO << "Build index from files";
    const auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{opt.kmer_size}}, seqan3::window_size{opt.window_size});
    PLOG_DEBUG << "Defined hash adaptor";
    PLOG_DEBUG << "Trying to initialize ibf with " << +summary.num_bins << " bins and " << opt.bits << " bits";
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{summary.num_bins},
                                         seqan3::bin_size{opt.bits},
                                         seqan3::hash_function_count{opt.num_hash}};
    PLOG_DEBUG << "Initialized ibf";
    InputStats input_stats;
    PLOG_DEBUG << "Defined input_stats";
#pragma omp parallel for
    for (const auto pair : summary.filepath_to_bin) {
        const auto& fasta_file = pair.first;
        const auto& bin = pair.second;

        PLOG_INFO << "Adding file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
#pragma omp critical
        input_stats.num_files += 1;

        auto record_count = 0;
        std::unordered_set<uint64_t> hashes;
        for (const auto & record : fin){
#pragma omp critical
            input_stats.records_per_bin[bin] += 1;
            const auto mh = record.sequence() | hash_adaptor | std::views::common;
            hashes.insert( mh.begin(), mh.end() );
            record_count++;
        }
#pragma omp critical
        for (auto && value : hashes){
            ibf.emplace(value, seqan3::bin_index{bin});
            input_stats.hashes_per_bin[bin] += 1;
        }

        PLOG_INFO << "Added file " << fasta_file << " with " << record_count << " records and " << hashes.size() << " hashes to bin " << +bin << std::endl;
    }

    if (opt.optimize){
        return optimize_ibf(opt, summary, input_stats, ibf);
    } else {
        return Index(opt, summary, stats, ibf);
    }
}

Index optimize_ibf(const IndexArguments& opt, const InputSummary& summary, const InputStats& stats, const seqan3::interleaved_bloom_filter<seqan3::uncompressed>& ibf)
{
    PLOG_INFO << "Optimize index bin layout";
    if (stats.hashes_per_bin.size() == 0) {
        return Index(opt, summary, stats, ibf);
    }

    auto new_summary = InputSummary();
    auto new_stats = InputStats();
    new_stats.num_files = stats.num_files;

    auto sorted_pairs = stats.bins_by_size();
    auto max_num_hashes = sorted_pairs.back().second;

    uint8_t next_bin = 0;
    std::unordered_map<std::string, uint8_t> last_bin;
    std::unordered_map<uint8_t, uint8_t> bin_to_bucket_map;

    // update stats
    PLOG_INFO << "Reassign bins";
    for (const auto & pair : sorted_pairs){
        const auto& bin = pair.first;
        const auto& category = summary.bin_to_category.at(bin);
        const auto& num_hashes = pair.second;

        auto assigned_bucket = next_bin;
        if (last_bin.find(category) != last_bin.end()){
            const auto& last_bin_in_category = last_bin[category];
            const auto& num_hashes_in_last_bin = new_stats.hashes_per_bin[last_bin_in_category];
            if (num_hashes_in_last_bin + num_hashes < max_num_hashes){
                assigned_bucket = last_bin_in_category;
            } else {
                next_bin++;
            }
        } else {
            next_bin++;
        }
        last_bin[category] = assigned_bucket;
        bin_to_bucket_map[bin] = assigned_bucket;
        PLOG_INFO << "Bin " << +bin << " assigned to bucket " << +assigned_bucket;
        new_stats.hashes_per_bin[assigned_bucket] += num_hashes;
        new_stats.records_per_bin[assigned_bucket] += stats.records_per_bin.at(bin);
    }

    // update summary
    PLOG_INFO << "Update summary";
    new_summary.num_bins = next_bin;
    new_summary.categories = summary.categories;
    for (const auto & pair : summary.filepath_to_bin){
        const auto & filepath = pair.first;
        const auto & bin = pair.second;
        const auto & bucket = bin_to_bucket_map[bin];
        new_summary.filepath_to_bin.emplace_back(std::make_pair(filepath, bucket));
        new_summary.bin_to_category[bucket] = summary.bin_to_category.at(bin);
    }

    // update ibf
    PLOG_INFO << "Create new IBF with " << +new_summary.num_bins << " bins and up to " << opt.bits << " bits";
    const auto num_bits = std::min(bin_size_in_bits(opt, max_num_hashes),opt.bits);
    PLOG_DEBUG << "Trying to initialize ibf with " << +new_summary.num_bins << " bins and " << num_bits << " bits";
    seqan3::interleaved_bloom_filter new_ibf{seqan3::bin_count{new_summary.num_bins},
                                         seqan3::bin_size{num_bits},
                                         seqan3::hash_function_count{opt.num_hash}};
    auto agent = ibf.membership_agent();
    for (auto bit=0; bit<opt.bits; bit++){
        auto & result = agent.bulk_contains(bit);
#pragma omp parallel for
        for (auto bin=0; bin<summary.num_bins; bin++){
            if (result[bin] > 0){
                const auto & bucket = bin_to_bucket_map[bin];
#pragma omp critical
                new_ibf.emplace(bit, seqan3::bin_index{bucket});
            }
        }
    }
    PLOG_INFO << "Filled in IBF";
    return Index(opt, summary, stats, ibf);

}

int index_main(IndexArguments & opt)
{
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

    LOG_INFO << "Running sifter index!";

    auto summary = parse_input_file(opt.input_file);
    auto stats = InputStats(); //estimate_index_size(input, opt);
    auto index = build_index(summary, stats, opt);
    store_index(opt.prefix, std::move(index));

    return 0;
}