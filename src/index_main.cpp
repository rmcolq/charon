#include <unordered_map>
#include <iostream>
#include <fstream>
#include <string>

#include "index_main.hpp"
#include "index.hpp"
#include "utils.hpp"

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

using namespace seqan3::literals;

void setup_index_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<IndexOptions>();
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

    index_subcommand->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    index_subcommand->callback([opt]() { sifter_index(*opt); });
}

int build_index(IndexOptions const& opt)
{
    auto log_level = plog::info;
    if (opt.verbosity == 1) {
        log_level = plog::debug;
    } else if (opt.verbosity > 1) {
        log_level = plog::verbose;
    }
    plog::init(log_level, "sifter.log", 1000000, 5);

    if (opt.window_size < opt.kmer_size) {
        throw std::logic_error("W must be greater than K");
    }
    if (opt.window_size <= 0) {
        throw std::logic_error("W must be a positive integer");
    }
    if (opt.kmer_size <= 0) {
        throw std::logic_error("K must be a positive integer");
    }

    LOG_INFO << "Hello log!";
    PLOG_VERBOSE << "verbose";
    PLOG_DEBUG << "debug";
    PLOG_INFO << "info";
    PLOG_WARNING << "warning";
    PLOG_ERROR << "error";
    PLOG_FATAL << "fatal";
    PLOG_NONE << "none";

    auto input = parse_input_file(opt.input_file);
    auto ibf = construct_ibf(input, opt.window_size, opt.kmer_size);

    return 0;
}

Input parse_input_file(const std::filesystem::path& input_file){
    PLOG_INFO << "Parsing input file " << input_file;
    std::unordered_map<std::string, uint8_t> filepath_to_bin;
    std::unordered_map<uint8_t, std::string> bin_to_name;
    std::unordered_map<std::string, uint8_t> name_to_bin;
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
                auto found_bin = name_to_bin.find(name);
                if (found_bin == name_to_bin.end()){
                    name_to_bin[name] = next_bin;
                    bin_to_name[next_bin] = name;
                    next_bin++;
                }
                auto bin = name_to_bin[name];
                filepath_to_bin[path] = bin;
            }
        }
    }
    input_ifstream.close();

    PLOG_INFO << "Found " << filepath_to_bin.size() << " files corresponding to " << bin_to_name.size() << " bins";

    return Input {filepath_to_bin, bin_to_name};
}

seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> construct_ibf(const Input& input, InputOptions const & opt){
    PLOG_INFO << "Building IBF from files";

    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{8u},
                                         seqan3::bin_size{8192u},
                                         seqan3::hash_function_count{2u}};

    IndexSummary summary;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{opt.kmer_size}}, seqan3::window_size{opt.window_size});

    for (const auto& pair : input.filepath_to_bin) {
        const auto& fasta_file = pair.first;
        const auto& bin = pair.second;

        PLOG_INFO << "Adding file " << fasta_file;
        seqan3::sequence_file_input fin{fasta_file};
        summary.files += 1;

        auto record_count = 0;
        for (auto & record : fin)
        {
            summary.records += 1;
            for (auto && value : record.sequence() | hash_adaptor)
                ibf.emplace(value, seqan3::bin_index{bin});
                summary.kmers_per_bin[bin] += 1;
                record_count++;
        }
        PLOG_INFO << "Added file " << fasta_file << " with " << record_count << " records to bin " << bin << std::endl;
    }
    summary.bins = summary.kmers_per_bin.size();

    // Construct an immutable, compressed Interleaved Bloom Filter.
    auto index = sifter_index(input, opt, summary, ibf);

    return index;
}
