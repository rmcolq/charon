#include <filesystem>
#include <fstream>
#include <cereal/archives/binary.hpp>

#include <index.hpp>

static inline void store_index(std::filesystem::path const & path, Index && index)
{
    PLOG_INFO << "Saving index to file " << path;
    std::ofstream os{path, std::ios::binary};
    cereal::BinaryOutputArchive oarchive{os};
    oarchive(index);
}