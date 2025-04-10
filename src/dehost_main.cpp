#include <unordered_set>
#include <iostream>
#include <algorithm>

#include "dehost_main.hpp"
#include "classify_stats.hpp"
#include "index.hpp"
#include "load_index.hpp"
#include "utils.hpp"
#include "version.h"

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

void setup_dehost_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<DehostArguments>();
    auto* dehost_subcommand = app.add_subcommand(
        "dehost", "Dehost read file into host and other using index.");

    dehost_subcommand->add_option("<fastaq>", opt->read_file, "Fasta/q file")
        ->required()
        ->transform(make_absolute)
        ->check(CLI::ExistingFile.description(""))
        ->type_name("FILE");

    dehost_subcommand->add_option("<fastaq>", opt->read_file2, "Paired Fasta/q file")
            ->transform(make_absolute)
            ->check(CLI::ExistingFile.description(""))
            ->type_name("FILE");

    dehost_subcommand
            ->add_option("--chunk_size", opt->chunk_size, "Read file is read in chunks of this size, to be processed in parallel within a chunk.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
        ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use.")
        ->type_name("INT")
        ->capture_default_str();

    dehost_subcommand->add_option("--db", opt->db, "Prefix for the index.")
        ->type_name("FILE")
        ->required()
        ->check(CLI::ExistingPath.description(""));

    dehost_subcommand->add_option("-e,--extract", opt->category_to_extract, "Reads from this category in the index will be extracted to file.")
            ->type_name("STRING");

    dehost_subcommand->add_option("--extract_file", opt->extract_file, "Fasta/q file for output")
            ->transform(make_absolute)
            ->check(CLI::NonexistentPath.description(""))
            ->type_name("FILE");

    dehost_subcommand->add_option("--extract_file2", opt->extract_file2, "Fasta/q file for output")
            ->transform(make_absolute)
            ->check(CLI::NonexistentPath.description(""))
            ->type_name("FILE");

    dehost_subcommand->add_option("-d,--dist", opt->dist, "Probability distribution to use for modelling.")
            ->type_name("STRING");

    dehost_subcommand
            ->add_option("--confidence", opt->confidence_threshold, "Minimum difference between the top 2 unique hit counts.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--min_hits", opt->confidence_threshold, "Minimum difference between the top 2 (non-unique) hit counts.")
            ->type_name("INT")
            ->capture_default_str();

    dehost_subcommand
            ->add_option("--min_diff", opt->min_proportion_difference, "Minimum difference between the proportion of (non-unique) kmers found in each category.")
            ->type_name("FLOAT")
            ->capture_default_str();

    dehost_subcommand->add_option("--log", opt->log_file, "File for log")
            ->transform(make_absolute)
            ->type_name("FILE");

    dehost_subcommand->add_flag(
        "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");

    // Set the function that will be called when this subcommand is issued.
    dehost_subcommand->callback([opt]() { dehost_main(*opt); });
}

void dehost_reads(const DehostArguments& opt, const Index& index){
    PLOG_INFO << "Dehosting file " << opt.read_file;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{index.kmer_size()}}, seqan3::window_size{index.window_size()});
    PLOG_VERBOSE << "Defined hash_adaptor";

    auto agent = index.agent();
    PLOG_VERBOSE << "Defined agent";

    seqan3::sequence_file_input<my_traits> fin{opt.read_file};
    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};

    using outfile_field_ids = decltype(fin)::field_ids;
    using outfile_format = decltype(fin)::valid_formats;

    auto result = Result<record_type, outfile_field_ids, outfile_format>(opt, index.summary());

    PLOG_DEBUG << "Defined Result with " << +index.num_bins() << " bins";

    for (auto && chunk : fin | seqan3::views::chunk(opt.chunk_size))
    {
        // You can use a for loop:
        for (auto & record : chunk)
        {
            records.push_back(std::move(record));
        }

#pragma omp parallel for firstprivate(agent, hash_adaptor) num_threads(opt.threads) shared(result)
        for (auto i=0; i<records.size(); ++i){

            const record_type & record = records[i];
            const auto read_id = split(record.id(), " ")[0];
            const auto read_length = record.sequence().size();
            if (read_length > std::numeric_limits<uint32_t>::max()){
                PLOG_WARNING << "Ignoring read " << record.id() << " as too long!";
                continue;
            }
            if (read_length == 0){
                PLOG_WARNING << "Ignoring read " << record.id() << " as has zero length!";
                continue;
            }
            auto qualities = record.base_qualities() | std::views::transform( [](auto quality) { return seqan3::to_phred(quality); });
            auto sum = std::accumulate(qualities.begin(), qualities.end(), 0);
            float mean_quality = 0;
            if (std::ranges::size(qualities) > 0)
                mean_quality = static_cast< float >( sum )/ static_cast< float >(std::ranges::size(qualities));
            PLOG_VERBOSE << "Mean quality of read  " << record.id() << " is " << mean_quality;

            auto read = ReadEntry(read_id, read_length, mean_quality, result.input_summary());
            for (auto && value : record.sequence() | hash_adaptor) {
                const auto & entry = agent.bulk_contains(value);
                read.update_entry(entry);
            }
            PLOG_VERBOSE << "Finished adding raw hash counts for read " << read_id;

            read.post_process(result.input_summary());
#pragma omp critical(add_read_to_results)
            result.add_read(read, record, true);
        }
        records.clear();
    }
    result.complete(true);
    result.print_summary();
}


void dehost_paired_reads(const DehostArguments& opt, const Index& index){
    PLOG_INFO << "Dehosting files " << opt.read_file << " and " << opt.read_file2;

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{index.kmer_size()}}, seqan3::window_size{index.window_size()});
    PLOG_VERBOSE << "Defined hash_adaptor";

    auto agent = index.agent();
    PLOG_VERBOSE << "Defined agent";

    seqan3::sequence_file_input<my_traits> fin1{opt.read_file};
    seqan3::sequence_file_input<my_traits> fin2{opt.read_file2};
    using record_type = decltype(fin1)::record_type;
    std::vector<record_type> records1{};
    std::vector<record_type> records2{};

    using outfile_field_ids = decltype(fin1)::field_ids;
    using outfile_format = decltype(fin1)::valid_formats;

    auto result = Result<record_type, outfile_field_ids, outfile_format>(opt, index.summary());

    PLOG_DEBUG << "Defined Result with " << +index.num_bins() << " bins";

    for (auto && chunk : fin1 | seqan3::views::chunk(opt.chunk_size))
    {
        for (auto & record : chunk)
        {
            records1.push_back(std::move(record));
        }

        // loop in the second file and get same amount of reads
        for ( auto& record2 : fin2 | std::views::take( opt.chunk_size ) )
        {
            records2.push_back( std::move( record2 ));
        }

#pragma omp parallel for firstprivate(agent, hash_adaptor) num_threads(opt.threads) shared(result)
        for (auto i=0; i<records1.size(); ++i){

            const auto & record1 = records1[i];
            const auto & record2 = records2[i];

            auto id1 = record1.id();
            id1.erase(id1.size() - 1);
            auto id2 = record2.id();
            id2.erase(id2.size() - 1);
            if (id1 != id2){
                std::cout << id1 << " " << id2;
                throw std::runtime_error("Your pairs don't match for read ids.");
            }
            const auto read_id = split(record1.id(), " ")[0];
            const auto read_length = record1.sequence().size() + record2.sequence().size();
            if (read_length > std::numeric_limits<uint32_t>::max()){
                PLOG_WARNING << "Ignoring read " << record1.id() << " as too long!";
                continue;
            }
            if (read_length == 0){
                PLOG_WARNING << "Ignoring read " << record1.id() << " as has zero length!";
                continue;
            }
            auto qualities1 = record1.base_qualities() | std::views::transform( [](auto quality) { return seqan3::to_phred(quality); });
            auto qualities2 = record2.base_qualities() | std::views::transform( [](auto quality) { return seqan3::to_phred(quality); });
            auto sum = std::accumulate(qualities1.begin(), qualities1.end(), 0);
            sum = std::accumulate(qualities2.begin(), qualities2.end(), sum);
            float mean_quality = 0;
            if (std::ranges::size(qualities1) + std::ranges::size(qualities2) > 0)
                mean_quality = static_cast< float >( sum )/ static_cast< float >(std::ranges::size(qualities1) + std::ranges::size(qualities2));
            PLOG_VERBOSE << "Mean quality of read  " << record1.id() << " is " << mean_quality;

            auto read = ReadEntry(read_id, read_length, mean_quality, result.input_summary());
            for (auto && value : record1.sequence() | hash_adaptor) {
                const auto & entry = agent.bulk_contains(value);
                read.update_entry(entry);
            }
            for (auto && value : record2.sequence() | hash_adaptor) {
                const auto & entry = agent.bulk_contains(value);
                read.update_entry(entry);
            }
            PLOG_VERBOSE << "Finished adding raw hash counts for read " << read_id;

            read.post_process(result.input_summary());
#pragma omp critical(add_read_to_results)
            result.add_paired_read(read, record1, record2);
        }
        records1.clear();
        records2.clear();
    }
    result.complete();
    result.print_summary();
}


int dehost_main(DehostArguments & opt)
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

    if (opt.read_file2 != "") {
        opt.is_paired = true;
        opt.min_length = 150;
    }

    auto args = opt.to_string();
    LOG_INFO << "Running charon dehost\n\nCharon version: " << SOFTWARE_VERSION << "\n" << args;

    auto index = Index();
    load_index(index, opt.db);
    auto host_index = index.get_host_index();
    LOG_INFO << "Found host at index " << +host_index << " in the index categories";

    opt.run_extract = (opt.category_to_extract != "");
    const auto categories = index.categories();
    if (opt.run_extract and std::find(categories.begin(), categories.end(), opt.category_to_extract) == categories.end())
    {
        std::string options = "";
        for (auto i: categories)
            options += i + " ";
        PLOG_ERROR << "Cannot extract " << opt.category_to_extract << ", please chose one of [ " << options << "]";
        return 1;
    } else if (opt.run_extract and opt.extract_file == ""){
        opt.extract_file = opt.read_file;
        if (opt.is_paired){
            opt.extract_file.replace_extension(opt.category_to_extract + opt.read_file.extension().string());
            opt.extract_file2 = opt.read_file2;
            opt.extract_file2.replace_extension(opt.category_to_extract + opt.read_file.extension().string());
        } else
            opt.extract_file.replace_extension(opt.category_to_extract + opt.read_file.extension().string());
    }

    if (opt.dist != "gamma" and opt.dist != "beta")
    {
        PLOG_ERROR << "Supported distributions are [gamma , beta]";
        return 1;
    }


    if (opt.is_paired)
        dehost_paired_reads(opt, index);
    else
        dehost_reads(opt, index);

    return 0;
}