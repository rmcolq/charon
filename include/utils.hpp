#pragma once

#include <filesystem>
#include <optional>
#include <vector>
#include <string>
#include <fstream>
#include <unordered_set>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/quality/phred94.hpp>

class IndexArguments;

struct MyTraits : seqan3::sequence_file_input_default_traits_dna {
    using quality_alphabet = seqan3::phred94;
};

// Used to transform paths to absolute paths — designed for use with CLI11 transform
[[nodiscard]] std::filesystem::path make_absolute(const std::filesystem::path& path);

// Split s by delimiter. Returns nullopt if delimiter is not found.
[[nodiscard]] std::optional<std::vector<std::string>> split(std::string_view s, std::string_view delimiter);

// Return the substring before the first occurrence of delimiter,
// or the whole string if delimiter is not found.
[[nodiscard]] std::string first_field(std::string_view s, std::string_view delimiter);

[[nodiscard]] bool ends_with(std::string_view str, std::string_view suffix);

[[nodiscard]] bool starts_with(std::string_view str, std::string_view prefix);

void store_hashes(std::string_view target,
                  const std::unordered_set<uint64_t>& hashes,
                  const std::filesystem::path& tmp_output_folder);

[[nodiscard]] std::vector<uint64_t> load_hashes(std::string_view target,
                                  const std::filesystem::path& tmp_output_folder);

void delete_hashes(const std::vector<uint8_t>& targets, const std::string& tmp_output_folder);

[[nodiscard]] size_t bin_size_in_bits(const IndexArguments& opt, const uint64_t num_elements);

size_t max_num_hashes_for_fpr(const IndexArguments& opt);

// Fixed: replaced compiler-internal __type_pack_element with the explicit resolved type
std::string sequence_to_string(const std::vector<seqan3::dna5>& input);

float get_compression_ratio(const std::string& sequence);

std::string get_extension(const std::filesystem::path& read_file);