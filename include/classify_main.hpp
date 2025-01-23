#ifndef CHARON_CLASSIFY_MAIN_H
#define CHARON_CLASSIFY_MAIN_H

#pragma once

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>

#include "result.hpp"

class Index;
struct ClassifyArguments;

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using quality_alphabet = seqan3::phred94;
};

void setup_classify_subcommand(CLI::App& app);

void classify_reads(const ClassifyArguments& opt, const Index& index, Result& result);

void extract(const ClassifyArguments& opt, Result& result);

int classify_main(ClassifyArguments & opt);


#endif // CHARON_CLASSIFY_MAIN_H
