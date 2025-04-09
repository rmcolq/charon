#include <filesystem>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <plog/Log.h>

#include <load_index.hpp>

void load_index(Index & index, std::filesystem::path const & path)
{
    PLOG_INFO << "Loading index from file " << path;
    std::ifstream is{path, std::ios::binary};
    cereal::BinaryInputArchive iarchive{is};
    iarchive(index);
    PLOG_INFO << "Index loaded";
    //PLOG_DEBUG << "Index has " << index.ibf().bin_count() << " bins and " << index.ibf().bit_size() << " bits";
}
