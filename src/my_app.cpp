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

int app_main()
{
    seqan3::interleaved_bloom_filter ibf{seqan3::bin_count{8u},
                                         seqan3::bin_size{8192u},
                                         seqan3::hash_function_count{2u}};

    auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}}, seqan3::window_size{5});

    auto fasta_file = "../../my.fasta";
    seqan3::sequence_file_input fin{fasta_file};

    for (auto & record : fin)
    {
        //seqan3::debug_stream << "ID:  " << record.id() << '\n';
        //seqan3::debug_stream << "SEQ: " << record.sequence() << '\n';
        // a quality field also exists, but is not printed, because we know it's empty for FASTA files.
        for (auto && value : record.sequence() | hash_adaptor)
            ibf.emplace(value, seqan3::bin_index{0u});
    }

    auto const sequence1 = "ACTGACTGACTGATC"_dna4;

    // Construct an immutable, compressed Interleaved Bloom Filter.
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf2{ibf};

    auto agent = ibf2.counting_agent();

    // Count all 5-mers of sequence1 for all bins
    seqan3::debug_stream << agent.bulk_count(sequence1 | hash_adaptor) << '\n'; // [11,0,0,0,9,0,0,0]

    // Search for specific values
    seqan3::debug_stream << agent.bulk_count(std::views::iota(0u, 1024u)) << '\n'; // [0,0,0,0,7,0,0,10]

}
