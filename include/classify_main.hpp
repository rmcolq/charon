#ifndef SIFTER_CLASSIFY_MAIN_H
#define SIFTER_CLASSIFY_MAIN_H

#pragma once

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>

#include "result.hpp"

class Index;


/// Collection of all options of index subcommand.
struct ClassifyArguments {
    // IO options
    std::string read_file;
    std::string db;
    uint8_t chunk_size { 100 };

    // General options
    std::string log_file {"sifter.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };
};

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using quality_alphabet = seqan3::phred94;
};

void setup_classify_subcommand(CLI::App& app);

Result classify_reads(const Index& index, const ClassifyArguments& opt);

int classify_main(ClassifyArguments & opt);


#endif // SIFTER_CLASSIFY_MAIN_H
