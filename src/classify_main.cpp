#include <unordered_map>
#include <unordered_set>
#include <iostream>

#include "classify_main.hpp"
#include "index.hpp"
#include "load_index.hpp"
#include "utils.hpp"

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>


void setup_classify_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<ClassifyArguments>();
    auto* classify_subcommand = app.add_subcommand(
        "classify", "Classify read file using index.");

    classify_subcommand->add_option("<fastaq>", opt->read_file, "Fasta/q file")
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    classify_subcommand
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use.")
        ->type_name("INT")
        ->capture_default_str();

    classify_subcommand->add_option("--db", opt->db, "Prefix for the index.")
        ->type_name("FILE")
        ->check(CLI::ExistingPath.description(""));

    classify_subcommand->add_option("--log", opt->log_file, "File for log")
            ->transform(make_absolute)
            ->type_name("FILE");

    classify_subcommand->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    classify_subcommand->callback([opt]() { classify_main(*opt); });
}

Result classify_reads(const Index& index, const ClassifyArguments& opt){
    PLOG_INFO << "Classifying file " << opt.read_file;
    auto result = Result(index.num_bins());
    PLOG_INFO << "Defined Result";

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{index.kmer_size()}}, seqan3::window_size{index.window_size()});
    PLOG_INFO << "Defined hash_adaptor";
    auto agent = index.agent();
    PLOG_INFO << "Defined agent";

    PLOG_INFO << "Classifying file " << opt.read_file;
    seqan3::sequence_file_input fin{opt.read_file};

    for (auto & record : fin)
    {
        PLOG_INFO << "Processing read " << record.id();
        auto read_id = split(record.id(), " ")[0];
        for (auto && value : record.sequence() | hash_adaptor)
        {
            auto & entry = agent.bulk_contains(value);
            result.update_entry(read_id, entry);
        }
        result.print_result(read_id);
    }
    PLOG_INFO << "Read all reads";
    return result;
}


int classify_main(ClassifyArguments & opt)
{
    auto log_level = plog::info;
    if (opt.verbosity == 1) {
        log_level = plog::debug;
    } else if (opt.verbosity > 1) {
        log_level = plog::verbose;
    }
    plog::init(log_level, opt.log_file.c_str(), 1000000, 5);

    if (!ends_with(opt.db, ".idx")) {
        opt.db += ".idx";
    }

    LOG_INFO << "Hello log!";
    PLOG_VERBOSE << "verbose";
    PLOG_DEBUG << "debug";
    PLOG_INFO << "info";
    PLOG_WARNING << "warning";
    PLOG_ERROR << "error";
    PLOG_FATAL << "fatal";
    PLOG_NONE << "none";

    auto index = Index();
    load_index(index, opt.db);

    auto result = classify_reads(index, opt);

    return 0;
}