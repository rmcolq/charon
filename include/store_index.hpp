#include <filesystem>
#include <fstream>

#include <index.hpp>

static inline void store_index(std::filesystem::path const & path, Index && index)
{
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}