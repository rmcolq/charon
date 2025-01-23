#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <algorithm>

#include "classify_main.hpp"
#include "index.hpp"
#include "load_index.hpp"
#include "utils.hpp"

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/sequence_file/record.hpp>
#include <seqan3/utility/views/enforce_random_access.hpp>
#include <seqan3/alphabet/quality/phred_base.hpp>
#include <seqan3/utility/range/concept.hpp>
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
            ->add_option("--chunk_size", opt->chunk_size, "Read file is read in chunks of this size, to be processed in parallel within a chunk.")
            ->type_name("INT")
            ->capture_default_str();

    classify_subcommand
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use.")
        ->type_name("INT")
        ->capture_default_str();

    classify_subcommand->add_option("--db", opt->db, "Prefix for the index.")
        ->type_name("FILE")
        ->required()
        ->check(CLI::ExistingPath.description(""));

    classify_subcommand->add_option("-e,--extract", opt->category_to_extract, "Reads from this category in the index will be extracted to file.")
            ->type_name("STRING");

    classify_subcommand->add_option("--extract_file", opt->extract_file, "Fasta/q file for output")
            ->transform(make_absolute)
            ->check(CLI::NonexistentPath.description(""))
            ->type_name("FILE");

    classify_subcommand->add_option("--log", opt->log_file, "File for log")
            ->transform(make_absolute)
            ->type_name("FILE");

    classify_subcommand->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    classify_subcommand->callback([opt]() { classify_main(*opt); });
}

void classify_reads(const ClassifyArguments& opt, const Index& index, Result& result){
    PLOG_INFO << "Classifying file " << opt.read_file;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{index.kmer_size()}}, seqan3::window_size{index.window_size()});
    PLOG_VERBOSE << "Defined hash_adaptor";

    auto agent = index.agent();
    PLOG_VERBOSE << "Defined agent";

    seqan3::sequence_file_input<my_traits> fin{opt.read_file};
    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};

    for (auto && chunk : fin | seqan3::views::chunk(opt.chunk_size))
    {
        // You can use a for loop:
        for (auto & record : chunk)
        {
            records.push_back(std::move(record));
        }

#pragma omp master
        result.check_entries_size(opt.chunk_size);
#pragma omp barrier

#pragma omp parallel for firstprivate(agent, hash_adaptor) num_threads(opt.threads) shared(result)
        for (auto i=0; i<records.size(); ++i){

            const auto & record = records[i];
            //PLOG_INFO << "Processing read " << record.id();
            const auto read_id = split(record.id(), " ")[0];
            const auto read_length = record.sequence().size();
            if (read_length > std::numeric_limits<uint32_t>::max()){
                PLOG_INFO << "Ignoring read " << record.id() << " as too long!";
                continue;
            }
            //auto read_quality = std::max(record.base_qualities());
            result.add_read(read_id, read_length);
            for (auto && value : record.sequence() | hash_adaptor) {
                const auto & entry = agent.bulk_contains(value);
                result.update_read(read_id, entry);
            }
            PLOG_VERBOSE << "Finished adding raw hash counts for read " << read_id;
            result.post_process_read(read_id);
        }
        records.clear();
    }
    result.complete();
    result.print_summary();
}

void extract(const ClassifyArguments& opt, Result& result)
{
    LOG_INFO << "Extracting reads matching category " << opt.category_to_extract;
    const auto category = result.category_index(opt.category_to_extract);
    LOG_DEBUG << "Corresponds to index " << +category;

    seqan3::sequence_file_input<my_traits> fin{opt.read_file};
    seqan3::sequence_file_output fout{opt.extract_file};

    // std::views::filter takes a function object (a lambda in this case) as input that returns a boolean
    auto call_matches_category_filter = std::views::filter(
        [&category, &result](auto const & record)
        {
            auto read_id = split(record.id(), " ")[0];;
            return result.call(read_id) == category;
        });

    auto total = 0;
    for (auto & record : fin | call_matches_category_filter)
    {
        fout.push_back(record);
        total += 1;
    }

    LOG_INFO << "Wrote " << total << " reads to file " << opt.extract_file.c_str() ;
}


int classify_main(ClassifyArguments & opt)
{
    auto log_level = plog::info;
    if (opt.verbosity == 1) {
        log_level = plog::debug;
    } else if (opt.verbosity > 1) {
        log_level = plog::verbose;
    }
    plog::init(log_level, opt.log_file.c_str(), 10000000, 5);

    if (!ends_with(opt.db, ".idx")) {
        opt.db += ".idx";
    }

    LOG_INFO << "Running charon classify";

    auto index = Index();
    load_index(index, opt.db);

    bool run_extract = (opt.category_to_extract != "");
    const auto categories = index.categories();
    if (run_extract and std::find(categories.begin(), categories.end(), opt.category_to_extract) == categories.end())
    {
        std::string options = "";
        for (auto i: categories)
            options += i + " ";
        PLOG_ERROR << "Cannot extract " << opt.category_to_extract << ", please chose one of [ " << options << "]";
        return 1;
    } else if (run_extract and opt.extract_file == ""){
        opt.extract_file = opt.read_file + ".extracted.fq";
    }

    auto result = Result(opt, index.summary());
    PLOG_DEBUG << "Defined Result with " << +index.num_bins() << " bins";

    classify_reads(opt, index, result);

    if (run_extract)
        extract(opt, result);

    return 0;
}