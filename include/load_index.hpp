#include <filesystem>
#include <fstream>
#include <cereal/archives/binary.hpp>

#include <index.hpp>

void load_index(Index & index, std::filesystem::path const & path)
{
    PLOG_INFO << "Loading index from file " << path;
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(index);
    PLOG_INFO << "Index loaded";
}