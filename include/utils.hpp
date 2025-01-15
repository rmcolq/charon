#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#pragma once

#include <filesystem>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>

class IndexArguments;

// used to transform paths to absolute paths - designed to be used with CLI11 transform
std::filesystem::path make_absolute(std::filesystem::path);

std::vector<std::string> split(const std::string&, const std::string&);

bool ends_with(std::string str, std::string suffix);
bool starts_with(std::string str, std::string prefix);

void store_hashes( const std::string target,
                   const std::unordered_set< uint64_t >& hashes,
                   const std::string tmp_output_folder );
std::vector< uint64_t > load_hashes( const std::string target,
                                     const std::string tmp_output_folder );
void delete_hashes( const std::vector<uint8_t>& targets, const std::string tmp_output_folder );

size_t bin_size_in_bits(const IndexArguments & opt, const uint64_t & num_elements);

double max_num_hashes_for_fpr(const IndexArguments & opt);

#endif
