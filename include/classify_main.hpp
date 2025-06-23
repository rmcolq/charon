#ifndef CHARON_CLASSIFY_MAIN_H
#define CHARON_CLASSIFY_MAIN_H

#pragma once

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"

#include "result.hpp"

class Index;

struct ClassifyArguments;

void setup_classify_subcommand(CLI::App &app);

void classify_reads(const ClassifyArguments &opt, const Index &index);

int classify_main(ClassifyArguments &opt);


#endif // CHARON_CLASSIFY_MAIN_H
