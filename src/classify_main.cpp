#include <unordered_set>
#include <iostream>
#include <algorithm>
#include <cinttypes>

#include "classify_main.hpp"
#include "read_processor.hpp"
#include "classify_stats.hpp"
#include "index.hpp"
#include "load_index.hpp"
#include "utils.hpp"
#include "version.h"
#include "memory_manager.hpp"

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


void setup_classify_subcommand(CLI::App &app) {
    auto opt = std::make_shared<ClassifyArguments>();
    auto *classify_subcommand = app.add_subcommand(
            "classify", "Classify read file using index.");

    classify_subcommand->add_option("<fastaq>", opt->read_file, "Fasta/q file")
            ->required()
            ->transform(make_absolute)
            ->check(CLI::ExistingFile.description(""))
            ->type_name("FILE");

    classify_subcommand->add_option("<fastaq>", opt->read_file2, "Paired Fasta/q file")
            ->transform(make_absolute)
            ->check(CLI::ExistingFile.description(""))
            ->type_name("FILE");

    classify_subcommand->add_option("--db", opt->db, "Prefix for the index.")
            ->type_name("FILE")
            ->required()
            ->check(CLI::ExistingPath.description(""));

    classify_subcommand->add_option("-e,--extract", opt->category_to_extract,
                                    "Reads from this category in the index will be extracted to file.")
            ->type_name("STRING");

    classify_subcommand->add_option("-p,--prefix", opt->prefix, "Prefix for the output files.")
            ->type_name("FILE")
            ->check(CLI::NonexistentPath.description(""))
            ->default_str("<prefix>");

    classify_subcommand
            ->add_option("--chunk_size", opt->chunk_size,
                         "Read file is read in chunks of this size, to be processed in parallel within a chunk.")
            ->type_name("INT")
            ->capture_default_str();

    classify_subcommand
            ->add_option("--max-memory", opt->max_memory_gb,
                         "Maximum memory to use in GB (0 = auto-detect based on system memory).")
            ->type_name("GB")
            ->capture_default_str();

    classify_subcommand
            ->add_option("--lo_hi_threshold", opt->lo_hi_threshold,
                         "Threshold used during model fitting stage to decide if read should be used to train lo or hi distribution.")
            ->type_name("FLOAT")
            ->capture_default_str();

    classify_subcommand
            ->add_option("--num_reads_to_fit", opt->num_reads_to_fit,
                         "Number of reads to use to train each distribution in the model.")
            ->type_name("INT")
            ->capture_default_str();

    classify_subcommand->add_option("-d,--dist", opt->dist, "Probability distribution to use for modelling.")
            ->type_name("STRING");

    classify_subcommand
            ->add_option("--min_quality", opt->min_quality, "Minimum read quality to classify.")
            ->type_name("INT")
            ->capture_default_str();

    classify_subcommand
            ->add_option("--min_length", opt->min_length, "Minimum read length to classify.")
            ->type_name("INT")
            ->capture_default_str();

    classify_subcommand
            ->add_option("--min_compression", opt->min_compression,
                         "Minimum read gzip compression ratio to classify (a measure of how much information is in the read.")
            ->type_name("FLOAT")
            ->capture_default_str();

    classify_subcommand
            ->add_option("--confidence", opt->confidence_threshold,
                         "Minimum difference between the top 2 unique hit counts.")
            ->type_name("INT")
            ->capture_default_str();

    classify_subcommand
            ->add_option("--min_proportion_diff", opt->min_proportion_difference,
                         "Minimum difference between the proportion of (non-unique) kmers found in each category.")
            ->type_name("FLOAT")
            ->capture_default_str();

    classify_subcommand->add_option("--log", opt->log_file, "File for log")
            ->transform(make_absolute)
            ->type_name("FILE");

    classify_subcommand
            ->add_option("-t,--threads", opt->threads, "Maximum number of threads to use.")
            ->type_name("INT")
            ->capture_default_str();

    classify_subcommand->add_flag(
            "-v", opt->verbosity, "Verbosity of logging. Repeat for increased verbosity");


    // Set the function that will be called when this subcommand is issued.
    classify_subcommand->callback([opt]() { classify_main(*opt); });
}

void classify_reads(const ClassifyArguments &opt, const Index &index) {
    PLOG_INFO << "Classifying file " << opt.read_file;

    // Initialize memory manager and validate configuration
    MemoryManager mem_mgr(opt.max_memory_gb > 0 ? opt.max_memory_gb * 1024 : 0);
    
    // Calculate optimal chunk size if not explicitly set by user
    uint16_t effective_chunk_size = opt.chunk_size;
    if (opt.chunk_size == 100) {
        // User didn't override chunk_size (100 is the default), calculate optimal
        effective_chunk_size = static_cast<uint16_t>(
            mem_mgr.calculate_optimal_chunk_size(opt.threads)
        );
        PLOG_INFO << "Auto-calculated chunk size: " << effective_chunk_size 
                 << " (based on " << (opt.max_memory_gb > 0 ? opt.max_memory_gb : 75) 
                 << "% system memory limit)";
    }
    
    // Validate configuration safety
    if (!mem_mgr.is_configuration_safe(effective_chunk_size, opt.threads)) {
        PLOG_WARNING << "Configuration may exceed memory limits!";
        PLOG_WARNING << "Consider reducing --chunk-size or --threads";
    }
    
    // Display memory usage report
    PLOG_INFO << mem_mgr.get_memory_report(effective_chunk_size, opt.threads);

    // Pre-compute hash parameters for better cache locality
    const auto hash_adaptor = seqan3::views::minimiser_hash(
        seqan3::shape{seqan3::ungapped{index.kmer_size()}},
        seqan3::window_size{index.window_size()}
    );
    PLOG_VERBOSE << "Defined hash_adaptor";

    auto agent = index.agent();
    PLOG_VERBOSE << "Defined agent";

    seqan3::sequence_file_input<MyTraits> fin{opt.read_file};
    using record_type = decltype(fin)::record_type;
    std::vector<record_type> records{};
    records.reserve(effective_chunk_size);

    using outfile_field_ids = decltype(fin)::field_ids;
    using outfile_format = decltype(fin)::valid_formats;

    auto result = Result<record_type, outfile_field_ids, outfile_format>(opt, index.summary());

    PLOG_DEBUG << "Defined Result with " << +index.num_bins() << " bins";

    uint64_t total_reads_processed = 0;
    
    try {
        for (auto &&chunk: fin | seqan3::views::chunk(effective_chunk_size)) {
            records.clear();
            for (auto &record: chunk) {
                records.push_back(std::move(record));
            }

            // Use optimized batch processing with two-pass algorithm
            process_read_batch<record_type, decltype(result), decltype(hash_adaptor), decltype(agent)>(
                records,
                agent,
                hash_adaptor,
                result,
                opt.min_length,
                opt.min_quality,
                opt.threads,
                false,  // is_dehost = false for classify
                total_reads_processed
            );
        }
    } catch (const seqan3::unexpected_end_of_input& e) {
        PLOG_WARNING << "File appears truncated or corrupted: " << e.what();
        PLOG_WARNING << "Successfully processed " << total_reads_processed << " reads before error.";
        PLOG_WARNING << "Consider checking input file integrity with tools like 'fastq_validator'.";
    } catch (const std::exception& e) {
        PLOG_ERROR << "Unexpected error while reading file: " << e.what();
        PLOG_WARNING << "Successfully processed " << total_reads_processed << " reads before error.";
    }
    
    result.complete();
    result.print_summary();
}


void classify_paired_reads(const ClassifyArguments &opt, const Index &index) {
    PLOG_INFO << "Classifying files " << opt.read_file << " and " << opt.read_file2;

    // Pre-compute hash parameters
    const auto hash_adaptor = seqan3::views::minimiser_hash(
        seqan3::shape{seqan3::ungapped{index.kmer_size()}},
        seqan3::window_size{index.window_size()}
    );
    PLOG_VERBOSE << "Defined hash_adaptor";

    auto agent = index.agent();
    PLOG_VERBOSE << "Defined agent";

    seqan3::sequence_file_input<MyTraits> fin1{opt.read_file};
    seqan3::sequence_file_input<MyTraits> fin2{opt.read_file2};
    using record_type = decltype(fin1)::record_type;
    std::vector<record_type> records1{};
    std::vector<record_type> records2{};
    records1.reserve(opt.chunk_size);
    records2.reserve(opt.chunk_size);

    using outfile_field_ids = decltype(fin1)::field_ids;
    using outfile_format = decltype(fin1)::valid_formats;

    auto result = Result<record_type, outfile_field_ids, outfile_format>(opt, index.summary());

    PLOG_DEBUG << "Defined Result with " << +index.num_bins() << " bins";

    uint64_t total_reads_processed = 0;
    
    try {
        for (auto &&chunk: fin1 | seqan3::views::chunk(opt.chunk_size)) {
            records1.clear();
            records2.clear();
            for (auto &record: chunk) {
                records1.push_back(std::move(record));
            }

            // Get corresponding reads from second file
            for (auto &record2: fin2 | std::views::take(opt.chunk_size)) {
                records2.push_back(std::move(record2));
            }

            // Use optimized batch processing for paired reads
            process_paired_read_batch<record_type, decltype(result), decltype(hash_adaptor), decltype(agent)>(
                records1,
                records2,
                agent,
                hash_adaptor,
                result,
                opt.min_length,
                opt.threads,
                total_reads_processed
            );
        }
    } catch (const seqan3::unexpected_end_of_input& e) {
        PLOG_WARNING << "Paired file appears truncated or corrupted: " << e.what();
        PLOG_WARNING << "Successfully processed " << total_reads_processed << " read pairs before error.";
        PLOG_WARNING << "Consider checking input file integrity with tools like 'fastq_validator'.";
    } catch (const std::exception& e) {
        PLOG_ERROR << "Unexpected error while reading paired files: " << e.what();
        PLOG_WARNING << "Successfully processed " << total_reads_processed << " read pairs before error.";
    }
    
    result.complete();
    result.print_summary();
}


int classify_main(ClassifyArguments &opt) {
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
        opt.min_length = 80;
    }

    auto args = opt.to_string();
    LOG_INFO << "Running charon classify\n\nCharon version: " << SOFTWARE_VERSION << "\n" << args;

    auto index = Index();
    load_index(index, opt.db);

    opt.run_extract = (opt.category_to_extract != "");
    const auto categories = index.categories();
    if (opt.run_extract and opt.category_to_extract != "all" and
        std::find(categories.begin(), categories.end(), opt.category_to_extract) == categories.end()) {
        std::string options = "";
        for (auto i: categories)
            options += i + " ";
        PLOG_ERROR << "Cannot extract " << opt.category_to_extract << ", please chose one of [ all " << options << "]";
        return 1;
    } else if (opt.run_extract) {
        if (opt.prefix == "")
            opt.prefix = "charon";
        std::vector<std::string> to_extract;
        if (opt.category_to_extract == "all")
            to_extract = categories;
        else
            to_extract.push_back(opt.category_to_extract);
        const auto extension = get_extension(opt.read_file);
        for (const auto &category: to_extract) {
            const auto category_index = index.get_category_index(category);
            if (opt.is_paired) {
                opt.extract_category_to_file[category_index].push_back(
                        opt.prefix + "_" + category + "_1" + extension + ".gz");
                opt.extract_category_to_file[category_index].push_back(
                        opt.prefix + "_" + category + "_2" + extension + ".gz");
            } else {
                opt.extract_category_to_file[category_index].push_back(opt.prefix + "_" + category + extension + ".gz");
            }
        }
    }

    if (opt.dist != "gamma" and opt.dist != "beta") {
        PLOG_ERROR << "Supported distributions are [gamma , beta]";
        return 1;
    }


    if (opt.is_paired)
        classify_paired_reads(opt, index);
    else
        classify_reads(opt, index);

    return 0;
}