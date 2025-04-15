#ifndef CHARON_DEHOST_MAIN_H
#define CHARON_DEHOST_MAIN_H

#pragma once

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>

#include "result.hpp"

class Index;
struct DehostArguments;

void setup_dehost_subcommand(CLI::App& app);

void dehost_reads(const DehostArguments& opt, const Index& index);

int dehost_main(DehostArguments & opt);


#endif // CHARON_DEHOST_MAIN_H
