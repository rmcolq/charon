#include "utils.hpp"
#include "index_main.hpp"

#include <plog/Log.h>
#include <gzip/compress.hpp>

std::filesystem::path make_absolute(const std::filesystem::path& path) {
    return std::filesystem::absolute(path);
}

std::optional<std::vector<std::string>> split(const std::string& s, const std::string& delimiter) {
    if (s.find(delimiter) == std::string::npos) {
        return std::nullopt;
    }
    std::vector<std::string> substrings;
    size_t start = 0;
    size_t end;
    while ((end = s.find(delimiter, start)) != std::string::npos) {
        substrings.push_back(s.substr(start, end - start));
        start = end + delimiter.size();
    }
    substrings.push_back(s.substr(start));
    return substrings;
}

std::string first_field(const std::string& s, const std::string& delimiter) {
    const size_t pos = s.find(delimiter);
    if (pos == std::string::npos) {
        return s;
    }
    return s.substr(0, pos);
}

bool ends_with(const std::string& str, const std::string& suffix) {
    return str.size() >= suffix.size()
        && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool starts_with(const std::string& str, const std::string& prefix) {
    return str.size() >= prefix.size()
        && str.compare(0, prefix.size(), prefix) == 0;
}

void store_hashes(std::string_view target,
                  const std::unordered_set<uint64_t>& hashes,
                  const std::filesystem::path& tmp_output_folder) {
    const std::filesystem::path outf = tmp_output_folder / (std::string{target} + ".min");
    std::ofstream outfile{outf, std::ios::binary | std::ios::app};
    if (!outfile.is_open()) {
        throw std::runtime_error{"Failed to open file for writing: " + outf.string()};
    }
    for (const auto& h : hashes) {
        outfile.write(reinterpret_cast<const char*>(&h), sizeof(h));
    }
    // outfile automatically closed by RAII
}

std::vector<uint64_t> load_hashes(std::string_view target,
                                  const std::filesystem::path& tmp_output_folder) {
    const std::filesystem::path file = tmp_output_folder / (std::string{target} + ".min");
    std::ifstream infile{file, std::ios::binary};
    if (!infile.is_open()) {
        throw std::runtime_error{"Failed to open file for reading: " + file.string()};
    }
    
    // Get file size to reserve space
    infile.seekg(0, std::ios::end);
    const auto file_size = infile.tellg();
    infile.seekg(0, std::ios::beg);
    
    std::vector<uint64_t> hashes;
    const size_t expected_count = file_size / sizeof(uint64_t);
    hashes.reserve(expected_count);  // Prevent reallocations
    
    uint64_t hash;
    while (infile.read(reinterpret_cast<char*>(&hash), sizeof(hash))) {
        hashes.push_back(hash);
    }
    return hashes;
}

void delete_hashes(const std::vector<uint8_t>& targets, const std::string& tmp_output_folder) {
    for (const auto& target : targets) {
        std::filesystem::path outf{tmp_output_folder};
        outf += "/" + std::to_string(target) + ".min";
        if (std::filesystem::exists(outf)) {
            std::filesystem::remove(outf);
        }
    }
    if (std::filesystem::is_empty(tmp_output_folder)) {
        std::filesystem::remove(tmp_output_folder);
    }
}

// TODO: consider passing num_hash, max_fpr, and bits as plain scalars for testability
size_t bin_size_in_bits(const IndexArguments& opt, const uint64_t num_elements) {
    assert(opt.num_hash > 0);
    assert(opt.max_fpr > 0.0);
    assert(opt.max_fpr < 1.0);

    double const numerator{-static_cast<double>(num_elements * opt.num_hash)};
    double const denominator{std::log(1 - std::exp(std::log(opt.max_fpr) / opt.num_hash))};
    double const result{std::ceil(numerator / denominator)};

    if (result > opt.bits) {
        PLOG_WARNING << "Require " << +result << " bits for max_fpr " << opt.max_fpr
                     << " but only have " << +opt.bits << " bits available";
        return opt.bits;
    }
    return static_cast<size_t>(result);
}

// TODO: same as bin_size_in_bits — consider plain scalar parameters for testability
size_t max_num_hashes_for_fpr(const IndexArguments& opt) {
    assert(opt.bits > 0);
    assert(opt.max_fpr > 0.0);
    assert(opt.max_fpr < 1.0);

    double const numerator{-static_cast<double>(opt.bits / opt.num_hash)};
    double const denominator{std::log(1 - std::exp(std::log(opt.max_fpr) / static_cast<double>(opt.num_hash)))};
    double const result{std::floor(numerator * denominator)};

    return static_cast<size_t>(result);
}

std::string sequence_to_string(const std::vector<seqan3::dna5>& input) {
    std::string str;
    str.reserve(input.size());
    for (const auto& c : input) {
        str.push_back(seqan3::to_char(c));
    }
    return str;
}

float get_compression_ratio(const std::string& sequence) {
    const char* pointer = sequence.data();
    std::size_t initial_size = sequence.size();

    std::string compressed_data = gzip::compress(pointer, initial_size);
    std::size_t compressed_size = compressed_data.size();

    assert(compressed_size != 0);
    return static_cast<float>(compressed_size) / static_cast<float>(initial_size);
}

std::string get_extension(const std::filesystem::path& read_file) {
    auto ext = read_file.extension().string();
    if (ext == ".gz") {
        std::filesystem::path short_read_file = read_file.stem();
        ext = short_read_file.extension().string();
    }
    return ext;
}