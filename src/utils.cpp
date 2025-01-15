#include "utils.hpp"
#include "index_main.hpp"

#include <plog/Log.h>

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
        outf += "/" + target + ".min";
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
        file += "/" + target + ".min";
        std::ifstream           infile{ file, std::ios::binary };
        while ( infile.read( reinterpret_cast< char* >( &hash ), sizeof( hash ) ) )
            hashes.push_back( hash );
        return hashes;
    }

void delete_hashes( const std::vector<uint8_t>& targets, const std::string tmp_output_folder )
{
    /*
     * delete hashes from disk
     */
    for ( const auto & target: targets )
    {
        std::filesystem::path outf{ tmp_output_folder };
        outf += "/" + std::to_string(target) + ".min";
        if ( std::filesystem::exists( outf ) )
            std::filesystem::remove( outf );
    }
    if ( std::filesystem::is_empty( tmp_output_folder ) )
        std::filesystem::remove(tmp_output_folder);
}

size_t bin_size_in_bits(const IndexArguments & opt, const uint64_t & num_elements)
{
    assert(opt.num_hash > 0);
    assert(opt.max_fpr > 0.0);
    assert(opt.max_fpr < 1.0);

    double const numerator{-static_cast<double>(num_elements * opt.num_hash)};
    double const denominator{std::log(1 - std::exp(std::log(opt.max_fpr) / opt.num_hash))};
    double const result{std::ceil(numerator / denominator)};

    if (result > opt.bits){
        PLOG_WARNING << "Require " << +result << " bits for max_fpr " << opt.max_fpr << " but only have " << +opt.bits << " bits available";
        return opt.bits;
    }
    return result;
}

double max_num_hashes_for_fpr(const IndexArguments & opt)
{
    assert(opt.bits > 0);
    assert(opt.max_fpr > 0.0);
    assert(opt.max_fpr < 1.0);

    double const numerator{-static_cast<double>(opt.num_hash / opt.bits)};
    double const denominator{std::log(1 - std::exp(std::log(opt.max_fpr) / opt.num_hash))};
    double const result{std::ceil(numerator / denominator)};

    return result;
}