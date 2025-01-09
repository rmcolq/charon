#include "utils.hpp"

std::filesystem::path make_absolute(std::filesystem::path path) { return std::filesystem::absolute(path); }

std::vector<std::string> split(const std::string& s, const std::string& delimiter){
    std::vector<std::string> substrings;

    int start, end = -1*delimiter.size();
    do {
        start = end + delimiter.size();
        end = s.find(delimiter, start);
        substrings.push_back(s.substr(start, end - start));
    } while (end != -1);

    return substrings;
}

bool ends_with(std::string str, std::string suffix)
{
    return str.size() >= suffix.size() && str.compare(str.size()-suffix.size(), suffix.size(), suffix) == 0;
}

bool starts_with(std::string str, std::string prefix)
{
    return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
}

void store_hashes( const std::string target,
                   const std::unordered_set< uint64_t >& hashes,
                   const std::string tmp_output_folder )
    {
        /*
         * store hashes from set to disk in the specified folder (or current folder ".")
         */
        std::filesystem::path outf{ tmp_output_folder };
        outf += target + ".min";
        std::ofstream outfile{ outf, std::ios::binary | std::ios::app };
        for ( auto&& h : hashes )
        {
            outfile.write( reinterpret_cast< const char* >( &h ), sizeof( h ) );
        }
        outfile.close();
    }


    std::vector< uint64_t > load_hashes( const std::string target,
                                         const std::string tmp_output_folder )
    {
        /*
         * load hashes file from disk and returns them in a vector
         */
        uint64_t                hash;
        std::vector< uint64_t > hashes;
        std::filesystem::path file{ tmp_output_folder };
        file += target + ".min";
        std::ifstream           infile{ file, std::ios::binary };
        while ( infile.read( reinterpret_cast< char* >( &hash ), sizeof( hash ) ) )
            hashes.push_back( hash );
        return hashes;
    }

void delete_hashes( const auto& targets, const std::string tmp_output_folder )
{
    /*
     * delete hashes from disk
     */
    for ( auto const& target: targets )
    {
        std::filesystem::path outf{ tmp_output_folder };
        outf += target + ".min";
        if ( std::filesystem::exists( outf ) )
            std::filesystem::remove( outf );
    }
}