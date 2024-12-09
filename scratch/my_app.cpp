// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <filesystem>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

using namespace seqan3::literals;

int main()
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{24u},
                                         seqan3::bin_size{409600000},
                                         seqan3::hash_function_count{2u}};

    const auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{19}}, seqan3::window_size{41});

    std::unordered_map<uint8_t, std::unordered_set<uint64_t>> hashes;
    auto fasta_file = "../../GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz";
    seqan3::sequence_file_input fin{fasta_file};

    uint8_t bin = 0;
    for (auto & record : fin)
    {
        const auto mh = record.sequence() | hash_adaptor | std::views::common;
        hashes[bin].insert( mh.begin(), mh.end() );

        for (auto && value : hashes[bin]){
            ibf.emplace(value, seqan3::bin_index{bin});
        }
        seqan3::debug_stream << "Add ID:  " << record.id() << " to bin " << +bin << " with " << +hashes[bin].size() << " hashes\n";
        bin++;
    }

    // Construct an immutable, compressed Interleaved Bloom Filter.
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf2{ibf};

    auto agent = ibf2.counting_agent<uint32_t>();

    for (const auto & [bin,values] : hashes) {
        seqan3::debug_stream << "bin:  " << bin << '\t';

        // Count all 5-mers of sequence1 for all bins
        seqan3::debug_stream << agent.bulk_count(values) << '\n'; // [11,0,0,0,9,0,0,0]
    }
}
