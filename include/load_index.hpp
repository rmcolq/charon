#include <filesystem>
#include <fstream>

#include <index.hpp>

void load_index(Index & index, std::filesystem::path const & path)
{
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(index);
}