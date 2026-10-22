#ifndef CHARON_FASTA_CHUNKS_H
#define CHARON_FASTA_CHUNKS_H

#pragma once

#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

// Whether for_each_fasta_chunk can read this file: FASTA, uncompressed or gzip-compressed (as seqan3 names them).
inline bool is_chunkable_fasta(const std::filesystem::path &file) {
    auto path = file;
    if (path.extension() == ".gz" or path.extension() == ".bgzf")
        path = path.stem();
    const auto extension = path.extension().string();
    for (const auto *fasta: {".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn", ".fas"}) {
        if (extension == fasta)
            return true;
    }
    return false;
}

/**
 * @brief Streams the sequences of a FASTA file in overlapping chunks, so whole records (chromosomes) are never held.
 *
 * Calls on_chunk(chunk) with consecutive pieces of each record's sequence as dna5, each of
 * at most chunk_size + overlap letters (or overlap + 1, if more). Consecutive chunks of a record share `overlap`
 * letters, so with overlap = window size - 1 every window of a record lies entirely within one chunk, and the minimisers of the
 * chunks are exactly those of the record. Parsing matches seqan3's FASTA reader: '>' or ';' starts an id line,
 * whitespace and digits in sequences are skipped, letters must be valid dna15 and are converted to dna5.
 * Returns the number of records.
 */
template<typename on_chunk_t>
uint64_t for_each_fasta_chunk(const std::filesystem::path &file, const size_t chunk_size, const size_t overlap,
                              on_chunk_t &&on_chunk) {
    gzFile input = gzopen(file.c_str(), "rb");
    if (input == nullptr)
        throw std::runtime_error{"Failed to open file for reading: " + file.string()};
    gzbuffer(input, 1 << 20);

    std::vector<char> buffer(1 << 20);
    std::vector<seqan3::dna5> chunk;
    chunk.reserve(chunk_size + overlap);
    size_t new_letters = 0;  // letters added since the last chunk was passed on
    uint64_t record_count = 0;
    bool in_id = false;

    auto end_record = [&]() {
        if (new_letters > 0)
            on_chunk(chunk);
        chunk.clear();
        new_letters = 0;
    };

    int bytes_read;
    while ((bytes_read = gzread(input, buffer.data(), static_cast<unsigned>(buffer.size()))) > 0) {
        for (int i = 0; i < bytes_read; ++i) {
            const char c = buffer[i];
            if (in_id) {
                in_id = c != '\n';
            } else if (c == '>' or c == ';') {
                end_record();
                in_id = true;
                ++record_count;
            } else if (not(seqan3::is_space || seqan3::is_digit)(c)) {
                if (not seqan3::char_is_valid_for<seqan3::dna15>(c)) {
                    gzclose(input);
                    throw std::runtime_error{"Unexpected letter '" + std::string{c} + "' in " + file.string()};
                }
                chunk.push_back(seqan3::assign_char_to(c, seqan3::dna5{}));
                // Pass a chunk on only once it holds a whole window: minimiser_hash also yields a minimiser for a
                // sequence shorter than the window, which the whole record wouldn't have.
                if (++new_letters >= chunk_size and chunk.size() > overlap) {
                    on_chunk(chunk);
                    const auto keep = std::min(overlap, chunk.size());
                    std::copy(chunk.end() - static_cast<std::ptrdiff_t>(keep), chunk.end(), chunk.begin());
                    chunk.resize(keep);
                    new_letters = 0;
                }
            }
        }
    }
    const bool failed = bytes_read < 0;
    gzclose(input);
    if (failed)
        throw std::runtime_error{"Failed to read " + file.string()};
    end_record();
    return record_count;
}

#endif // CHARON_FASTA_CHUNKS_H
