#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <string>

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
        "index", "Build an index (HIBF) for a number of references split into a small number of bins.");

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

    // Set the function that will be called when this subcommand is issued.
    index_subcommand->callback([opt]() { index_main(*opt); });
}

InputFileMap parse_input_file(const std::filesystem::path& input_file){
    std::cout << "Parsing input file " << std::endl;
    PLOG_INFO << "Parsing input file " << input_file;
    std::vector<std::pair<std::string, uint8_t>> filepath_to_bin;
    std::unordered_map<uint8_t, std::string> bin_to_name;
    std::unordered_set<std::string> categories;
    uint8_t next_bin = 0;

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
                auto path = parts[0];
                auto name = make_absolute(parts[1]).string();
                bin_to_name[next_bin] = name;
                categories.insert(name);
                filepath_to_bin.emplace_back(std::make_pair(path, next_bin));
                next_bin++;
            }
        }
    }
    input_ifstream.close();

    PLOG_INFO << "Found " << filepath_to_bin.size() << " files corresponding to " << categories.size() << " categories";
    std::cout << "Found " << filepath_to_bin.size() << " files corresponding to " << categories.size() << " categories" << std::endl;

    return InputFileMap {filepath_to_bin, bin_to_name};
}

InputSummary estimate_index_size(const InputFileMap& input, const IndexArguments& opt){
    PLOG_INFO << "Estimate size from files";
    InputSummary summary;

    std::unordered_map<uint8_t, uint64_t> lengths;

    for (const auto& pair : input.filepath_to_bin) {
        const auto& fasta_file = pair.first;
        const auto& bin = pair.second;

        PLOG_INFO << "Checking file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
        summary.num_files += 1;

        auto length = 0;
        for (auto & record : fin)
        {
            length += record.sequence().size();
            summary.records_per_bin[bin]++;
        }
        lengths[bin] += length;
    }

    // estimate max bin size
    uint64_t max_count = 0;
    for(auto kv : lengths) {
        summary.num_bins += 1;
        summary.hashes_per_bin[kv.first] = static_cast<uint64_t>(round((kv.second * 2)/(opt.window_size + 1)));
        max_count = std::max(max_count, static_cast<uint64_t>(round(1.1*summary.hashes_per_bin[kv.first])));
    }
    opt.bits = std::ceil( ( max_count * std::log( opt.max_fpr ) ) / std::log( 1.0 / std::pow( 2, std::log( 2 ) ) ) );

    PLOG_DEBUG << "num_bins: " << +summary.num_bins << ", num_files: " << summary.num_files;
    for (const auto& item: summary.records_per_bin){
        PLOG_DEBUG << "Bin " << +item.first << " has " << item.second << " records";
    }
    for (const auto& item: summary.hashes_per_bin){
        PLOG_DEBUG << "Bin " << +item.first << " has " << item.second << " hashes";
    }
    PLOG_DEBUG << "using: " << opt.bits << " to achieve a max_fpr of " << opt.max_fpr;

    return summary;
}

InputSummary summarise_input(const InputFileMap& input, const IndexArguments& opt){
    PLOG_INFO << "Estimate size and summary information from files";
    InputSummary summary;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{opt.kmer_size}}, seqan3::window_size{opt.window_size});
    std::unordered_map<uint8_t, std::unordered_set<uint64_t>> hashes;

    for (const auto& pair : input.filepath_to_bin) {
        const auto& fasta_file = pair.first;
        const auto& bin = pair.second;

        bool bin_seen_before = hashes.contains(bin);
        if (bin_seen_before){
            continue;
        }

        PLOG_INFO << "Checking file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
        //summary.num_files += 1;

        auto record_count = 0;
        for (auto & record : fin)
        {
            //summary.records_per_bin[bin] += 1;
            const auto mh = record.sequence() | hash_adaptor | std::views::common;
            hashes[bin].insert( mh.begin(), mh.end() );
            record_count++;
        }
        //PLOG_INFO << "Added file " << fasta_file << " with " << record_count << " records to bin " << bin << std::endl;
    }
    for (const auto& pair : hashes) {
        const auto& bin = pair.first;
        const auto& num_hashes = pair.second.size();
        summary.hashes_per_bin[bin] = num_hashes;
        summary.num_bins += 1;
    }

    // estimate max bin size
    uint64_t max_count = 0;
    for(auto kv : summary.hashes_per_bin) {
        max_count = std::max(max_count, static_cast<uint64_t>(round(kv.second*1.1)));
    }
    opt.bits = max_count;

    std::cout << "num_bins: " << +summary.num_bins << ", num_files: " << summary.num_files << std::endl;
    for (const auto& item: summary.records_per_bin){
        std::cout << "Bin " << +item.first << " has " << item.second << " records" << std::endl;
    }
    for (const auto& item: summary.hashes_per_bin){
        std::cout << "Bin " << +item.first << " has " << item.second << " hashes" << std::endl;
    }

    return summary;
}

Index build_index(const InputFileMap& input, InputSummary& summary, const IndexArguments& opt)
{
    PLOG_INFO << "Build index from files";
    const auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{opt.kmer_size}}, seqan3::window_size{opt.window_size});
    PLOG_DEBUG << "Defined hash adaptor";
    PLOG_DEBUG << "Trying to initialize ibf with " << +summary.num_bins << " bins and " << opt.bits << " bits";
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{summary.num_bins},
                                         seqan3::bin_size{opt.bits},
                                         seqan3::hash_function_count{2u}};
    PLOG_DEBUG << "Initialized ibf";
    InputSummary index_summary;
    PLOG_DEBUG << "Defined index_summary";
#pragma omp parallel for
    for (const auto pair : input.filepath_to_bin) {
        const auto& fasta_file = pair.first;
        const auto& bin = pair.second;

        PLOG_INFO << "Adding file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
#pragma omp critical
        index_summary.num_files += 1;

        auto record_count = 0;
        std::unordered_set<uint64_t> hashes;
        for (const auto & record : fin){
#pragma omp critical
            index_summary.records_per_bin[bin] += 1;
            const auto mh = record.sequence() | hash_adaptor | std::views::common;
            hashes.insert( mh.begin(), mh.end() );
            record_count++;
        }
#pragma omp critical
        for (auto && value : hashes){
            ibf.emplace(value, seqan3::bin_index{bin});
            index_summary.hashes_per_bin[bin] += 1;
        }

        PLOG_INFO << "Added file " << fasta_file << " with " << record_count << " records and " << hashes.size() << " hashes to bin " << +bin << std::endl;
    }

    return Index(opt, input, summary, ibf);
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

    auto input = parse_input_file(opt.input_file);
    auto summary = InputSummary(); //estimate_index_size(input, opt);
    summary.num_bins = input.bin_to_name.size();
    auto index = build_index(input, summary, opt);
    store_index(opt.prefix, std::move(index));

    return 0;
}