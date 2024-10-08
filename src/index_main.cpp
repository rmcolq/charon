#include "index_main.hpp"
#include <plog/Log.h>
#include "plog/Initializers/RollingFileInitializer.h"

void setup_index_subcommand(CLI::App& app)
{
    auto opt = std::make_shared<IndexOptions>();
    auto* index_subcommand = app.add_subcommand(
        "index", "Build an index (HIBF) for a number of references split into a small number of bins.");

    index_subcommand->add_option("<input>", opt->input_file, "Tab separated file with columns for filename and category")
        ->required()
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

int sifter_index(IndexOptions const& opt)
{
    auto log_level = plog::info;
    if (opt.verbosity == 1) {
        log_level = plog::debug;
    } else if (opt.verbosity > 1) {
        log_level = plog::verbose;
    }
    plog::init(log_level, "sifter.log", 1000000, 5);

    if (opt.window_size > opt.kmer_size) {
        throw std::logic_error("W must NOT be greater than K");
    }
    if (opt.window_size <= 0) {
        throw std::logic_error("W must be a positive integer");
    }
    if (opt.kmer_size <= 0) {
        throw std::logic_error("K must be a positive integer");
    }

    LOG_INFO << "Hello log!";
    return 0;
}
