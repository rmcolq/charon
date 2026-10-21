#ifndef CHARON_INDEX_MAIN_H
#define CHARON_INDEX_MAIN_H

#pragma once

#include <layout.hpp>

#include <omp.h>
#include <cstring>
#include <ankerl/unordered_dense.h>

#include "CLI11.hpp"
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include "index_arguments.hpp"

class Index;

class InputStats;

class InputSummary;


void setup_index_subcommand(CLI::App &app);

void apply_profile(IndexArguments &opt, bool window_set, bool kmer_set, bool fpr_set);

InputSummary parse_input_file(const std::filesystem::path &input_file);

InputStats count_and_store_hashes(const IndexArguments &opt, const InputSummary &summary);

Index build_index(const IndexArguments &opt, const InputSummary &summary, InputStats &stats,
                  const ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>> &bucket_to_bins_map);

void write_stats(const IndexArguments &opt, const InputSummary &summary, const InputStats &stats);

std::pair<ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>>, Layout>
balance_layout(const IndexArguments &opt, InputSummary &summary, InputStats &stats);

Index build_balanced_index(const IndexArguments &opt, const InputSummary &summary, InputStats &stats,
                           const ankerl::unordered_dense::map<uint8_t, std::vector<uint8_t>> &group_to_bins_map,
                           const Layout &layout);

int index_main(IndexArguments &opt);


#endif // CHARON_INDEX_MAIN_H
