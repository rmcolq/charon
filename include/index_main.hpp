#ifndef CHARON_INDEX_MAIN_H
#define CHARON_INDEX_MAIN_H

#pragma once

#include <omp.h>
#include <cstring>

#include "CLI11.hpp"
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include "index_arguments.hpp"

class Index;
class InputStats;
class InputSummary;


void setup_index_subcommand(CLI::App& app);

InputSummary parse_input_file(const std::filesystem::path& input_file);

InputStats count_and_store_hashes(const IndexArguments& opt, const InputSummary& summary);

std::unordered_map<uint8_t, std::vector<uint8_t>> optimize_layout(const IndexArguments& arguments, const InputSummary& summary, const InputStats& stats);

Index build_index(const IndexArguments& opt, const InputSummary& summary, InputStats& stats);


int index_main(IndexArguments & opt);


#endif // CHARON_INDEX_MAIN_H
