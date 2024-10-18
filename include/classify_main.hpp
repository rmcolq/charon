#ifndef SIFTER_CLASSIFY_MAIN_H
#define SIFTER_CLASSIFY_MAIN_H

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"


/// Collection of all options of index subcommand.
struct ClassifyArguments {
    // IO options
    std::string read_file;
    std::string db;

    // General options
    std::string log_file {"sifter.log"};
    uint8_t threads { 1 };
    uint8_t verbosity { 0 };
};

void setup_classify_subcommand(CLI::App& app);

int classify_main(ClassifyArguments & opt);


#endif // SIFTER_CLASSIFY_MAIN_H
