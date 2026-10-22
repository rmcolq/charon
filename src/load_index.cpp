#include <array>
#include <filesystem>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <plog/Log.h>

#include <load_index.hpp>

void load_index(Index &index, std::filesystem::path const &path) {
    PLOG_INFO << "Loading index from file " << path;
    std::ifstream is{path, std::ios::binary};
    std::array<char, Index::format_magic.size()> magic{};
    is.read(magic.data(), magic.size());
    const bool has_magic = is.gcount() == static_cast<std::streamsize>(magic.size()) and magic == Index::format_magic;
    if (not has_magic) {
        // index written by Charon <= v3: the archive starts at the beginning of the file
        is.clear();
        is.seekg(0);
    }
    cereal::BinaryInputArchive iarchive{is};
    index.load(iarchive, has_magic);
    PLOG_INFO << "Index loaded (" << (has_magic ? "v4" : "v3") << " format, "
              << (index.is_compressed() ? "compressed" : "uncompressed") << " IBF with "
              << index.technical_bins() << " technical bins of " << index.bin_size() << " bits, "
              << +index.num_bins() << " scoring bins)";
}
