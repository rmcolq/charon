#include "utils.hpp"
#include "index_main.hpp"

#include <plog/Log.h>
#include <gzip/compress.hpp>

std::filesystem::path make_absolute(std::filesystem::path path) { return std::filesystem::absolute(path); }

std::vector<std::string> split(const std::string &s, const std::string &delimiter) {
    std::vector<std::string> substrings;

    int start, end = -1 * delimiter.size();
    do {
        start = end + delimiter.size();
        end = s.find(delimiter, start);
        substrings.push_back(s.substr(start, end - start));
    } while (end != -1);

    return substrings;
}

bool ends_with(std::string str, std::string suffix) {
    return str.size() >= suffix.size() && str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool starts_with(std::string str, std::string prefix) {
    return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
}

void store_hashes(const std::string target,
                  const std::unordered_set<uint64_t> &hashes,
                  const std::string tmp_output_folder) {
    /*
     * store hashes from set to disk in the specified folder (or current folder ".")
     */
    std::filesystem::path outf{tmp_output_folder};
    outf += "/" + target + ".min";
    std::ofstream outfile{outf, std::ios::binary | std::ios::app};
    for (auto &&h: hashes) {
        outfile.write(reinterpret_cast< const char * >( &h ), sizeof(h));
    }
    outfile.close();
}


std::vector<uint64_t> load_hashes(const std::string target,
                                  const std::string tmp_output_folder) {
    /*
     * load hashes file from disk and returns them in a vector
     */
    uint64_t hash;
    std::vector<uint64_t> hashes;
    std::filesystem::path file{tmp_output_folder};
    file += "/" + target + ".min";
    std::ifstream infile{file, std::ios::binary};
    while (infile.read(reinterpret_cast< char * >( &hash ), sizeof(hash)))
        hashes.push_back(hash);
    return hashes;
}

void delete_hashes(const std::vector<uint8_t> &targets, const std::string tmp_output_folder) {
    /*
     * delete hashes from disk
     */
    for (const auto &target: targets) {
        std::filesystem::path outf{tmp_output_folder};
        outf += "/" + std::to_string(target) + ".min";
        if (std::filesystem::exists(outf))
            std::filesystem::remove(outf);
    }
    if (std::filesystem::is_empty(tmp_output_folder))
        std::filesystem::remove(tmp_output_folder);
}

size_t bin_size_in_bits(const IndexArguments &opt, const uint64_t &num_elements) {
    assert(opt.num_hash > 0);
    assert(opt.max_fpr > 0.0);
    assert(opt.max_fpr < 1.0);

    double const numerator{-static_cast<double>(num_elements * opt.num_hash)};
    double const denominator{std::log(1 - std::exp(std::log(opt.max_fpr) / opt.num_hash))};
    double const result{std::ceil(numerator / denominator)};

    if (result > opt.bits) {
        PLOG_WARNING << "Require " << +result << " bits for max_fpr " << opt.max_fpr << " but only have " << +opt.bits
                     << " bits available";
        return opt.bits;
    }
    return result;
}

size_t max_num_hashes_for_fpr(const IndexArguments &opt) {
    assert(opt.bits > 0);
    assert(opt.max_fpr > 0.0);
    assert(opt.max_fpr < 1.0);

    double const numerator{-static_cast<double>(opt.bits / opt.num_hash)};
    double const denominator{std::log(1 - std::exp(std::log(opt.max_fpr) / static_cast<double>(opt.num_hash)))};
    double const result{std::floor(numerator * denominator)};

    return result;
}

//std::__tuple_element_t<std::vector<seqan3::dna5>,std::vector<seqan3::phred94>>
std::string sequence_to_string(
        const __type_pack_element<0, std::vector<seqan3::dna5>, std::string, std::vector<seqan3::phred94>> &input) {
    std::vector<char> sequence{};
    for (auto c: input)
        sequence.push_back(seqan3::to_char(c));
    std::string str(sequence.begin(), sequence.end());
    return str;
}

float get_compression_ratio(const std::string &sequence) {
    const char *pointer = sequence.data();
    std::size_t initial_size = sequence.size();

    std::string compressed_data = gzip::compress(pointer, initial_size);
    const char *compressed_pointer = compressed_data.data();
    std::size_t compressed_size = compressed_data.size();

    assert(compressed_size != 0);
    return static_cast<double>(compressed_size) / static_cast<double>(initial_size);
}

std::string get_extension(const std::filesystem::path read_file) {
    auto ext = read_file.extension().string();
    if (ext == ".gz") {
        std::filesystem::path short_read_file = read_file.stem();
        ext = short_read_file.extension().string();
    }
    return ext;
}
