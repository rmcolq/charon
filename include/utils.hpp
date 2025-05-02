#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#pragma once

#include <filesystem>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>

class IndexArguments;

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using quality_alphabet = seqan3::phred94;
};

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

size_t max_num_hashes_for_fpr(const IndexArguments & opt);

std::string sequence_to_string(const __type_pack_element<0, std::vector<seqan3::dna5>, std::string, std::vector<seqan3::phred94>>& input);

float get_compression_ratio(const std::string& sequence);

std::string get_extension(const std::filesystem::path);
#endif
