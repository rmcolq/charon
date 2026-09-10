#ifndef CHARON_READ_PROCESSOR_H
#define CHARON_READ_PROCESSOR_H

#pragma once

#include <vector>
#include <string>
#include <cstdint>
#include <limits>
#include <numeric>
#include <algorithm>

#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/alphabet/quality/phred_base.hpp>
#include <seqan3/utility/range/concept.hpp>

#include <plog/Log.h>

#include "read_entry.hpp"
#include "index.hpp"
#include "utils.hpp"

/**
 * @brief Pre-computed read metadata for batch processing
 * 
 * This structure holds all the metadata extracted from a read before
 * hash computation, enabling two-pass processing for better performance.
 */
struct ReadMetadata {
    std::string read_id;
    uint32_t length;
    float mean_quality;
    float compression_ratio;
    bool is_valid;
    
    ReadMetadata() : length(0), mean_quality(0.0f), compression_ratio(0.0f), is_valid(false) {}
    
    ReadMetadata(std::string id, uint32_t len, float qual, float comp, bool valid)
        : read_id(std::move(id)), length(len), mean_quality(qual), 
          compression_ratio(comp), is_valid(valid) {}
};

/**
 * @brief Extract metadata from a sequence record
 * 
 * Performs early filtering to identify invalid reads before expensive hash computation.
 * 
 * @tparam record_type Type of the sequence record
 * @param record The sequence record to process
 * @param min_length Minimum read length threshold
 * @param min_quality Minimum mean quality threshold (optional, 0 = no filter)
 * @return ReadMetadata Extracted metadata with validity flag
 */
template<typename record_type>
ReadMetadata extract_read_metadata(const record_type& record, 
                                   uint32_t min_length = 0,
                                   float min_quality = 0.0f) {
    const auto read_id = first_field(record.id(), " ");
    const uint32_t read_length = std::ranges::size(record.sequence());
    
    // Early length filtering
    if (read_length > std::numeric_limits<uint32_t>::max() || read_length == 0) {
        if (read_length > std::numeric_limits<uint32_t>::max()) {
            PLOG_WARNING << "Ignoring read " << record.id() << " as too long!";
        } else {
            PLOG_WARNING << "Ignoring read " << record.id() << " as has zero length!";
        }
        return ReadMetadata(read_id, read_length, 0.0f, 0.0f, false);
    }
    
    // Compute mean quality
    auto qualities = record.base_qualities() |
                     std::views::transform([](auto quality) { return seqan3::to_phred(quality); });
    auto sum = std::accumulate(qualities.begin(), qualities.end(), 0);
    float mean_quality = 0.0f;
    if (std::ranges::size(qualities) > 0) {
        mean_quality = static_cast<float>(sum) / static_cast<float>(std::ranges::size(qualities));
    }
    PLOG_VERBOSE << "Mean quality of read " << read_id << " is " << mean_quality;
    
    // Compute compression ratio
    float compression_ratio = get_compression_ratio(sequence_to_string(record.sequence()));
    PLOG_VERBOSE << "Found compression ratio of read " << read_id << " is " << compression_ratio;
    
    // Quality filtering (optional)
    bool passes_quality_filter = (min_quality == 0.0f) || (mean_quality >= min_quality);
    bool passes_length_filter = read_length >= min_length;
    
    return ReadMetadata(read_id, read_length, mean_quality, compression_ratio, 
                       passes_quality_filter && passes_length_filter);
}

/**
 * @brief Compute all minimizer hashes for a sequence
 * 
 * Collects all hash values for a read in a contiguous vector for batch processing.
 * This enables better cache locality and SIMD optimization opportunities.
 * 
 * @tparam record_type Type of the sequence record
 * @tparam hash_adaptor_type Type of the hash adaptor
 * @param record The sequence record
 * @param hash_adaptor The minimizer hash adaptor
 * @return std::vector<uint64_t> Vector of hash values
 */
template<typename record_type, typename hash_adaptor_type>
std::vector<uint64_t> compute_read_hashes(const record_type& record, 
                                          const hash_adaptor_type& hash_adaptor) {
    std::vector<uint64_t> hashes;
    hashes.reserve(std::ranges::size(record.sequence()));  // Reserve to avoid reallocations
    
    for (auto &&value: record.sequence() | hash_adaptor) {
        hashes.push_back(value);
    }
    
    return hashes;
}

/**
 * @brief Process a single read with pre-computed metadata and hashes
 * 
 * This is the core processing function that:
 * 1. Creates a ReadEntry with pre-computed metadata
 * 2. Processes all hashes through the IBF agent
 * 3. Performs post-processing (counting, proportions, probabilities)
 * 4. Adds the result to the output
 * 
 * @tparam record_type Type of the sequence record
 * @tparam result_type Type of the result container
 * @param read_id Read identifier
 * @param read_length Read length
 * @param mean_quality Mean quality score
 * @param compression_ratio Compression ratio
 * @param hashes Pre-computed hash values
 * @param agent IBF membership agent
 * @param result Output result container
 * @param record Original sequence record (for output)
 * @param is_dehost Flag indicating if this is dehost mode (affects result handling)
 */
template<typename record_type, typename result_type, typename agent_type>
void process_single_read(const std::string& read_id,
                        uint32_t read_length,
                        float mean_quality,
                        float compression_ratio,
                        const std::vector<uint64_t>& hashes,
                        agent_type& agent,
                        result_type& result,
                        const record_type& record,
                        bool is_dehost = false) {
    // Create ReadEntry with pre-computed metadata
    auto read = ReadEntry(read_id, read_length, mean_quality, compression_ratio, 
                         result.input_summary());
    
    // Process all hashes - tight loop for better branch prediction
    for (const auto& hash_value : hashes) {
        const auto &entry = agent.bulk_contains(hash_value);
        read.update_entry(entry);
    }
    
    // Post-process: compute counts, proportions, and probabilities
    read.post_process(result.input_summary());
    
    // Add to results (thread-safe critical section)
    #pragma omp critical(add_read_to_results)
    {
        if constexpr (requires { result.add_read(read, record, is_dehost); }) {
            result.add_read(read, record, is_dehost);
        } else {
            result.add_read(read, record);
        }
    }
}

/**
 * @brief Process a batch of reads with optimized two-pass algorithm
 * 
 * PASS 1: Extract metadata and compute all hashes (sequential, cache-friendly)
 * PASS 2: Process hashes in parallel with OpenMP + SIMD hints
 * 
 * This provides:
 * - Better cache locality (metadata and hashes stored contiguously)
 * - Reduced memory fragmentation (pre-allocated vectors)
 * - SIMD vectorization opportunities
 * - Early filtering of invalid reads
 * 
 * @tparam record_type Type of the sequence records
 * @tparam result_type Type of the result container
 * @tparam hash_adaptor_type Type of the hash adaptor
 * @param records Vector of sequence records to process
 * @param agent IBF membership agent
 * @param hash_adaptor The minimizer hash adaptor
 * @param result Output result container
 * @param min_length Minimum read length threshold
 * @param min_quality Minimum mean quality threshold
 * @param num_threads Number of OpenMP threads to use
 * @param is_dehost Flag indicating if this is dehost mode
 * @param total_reads_processed Reference to counter for processed reads
 */
template<typename record_type, typename result_type, typename hash_adaptor_type, typename agent_type>
void process_read_batch(const std::vector<record_type>& records,
                       agent_type& agent,
                       const hash_adaptor_type& hash_adaptor,
                       result_type& result,
                       uint32_t min_length,
                       float min_quality,
                       int num_threads,
                       bool is_dehost,
                       uint64_t& total_reads_processed) {
    
    // PASS 1: Extract metadata and compute hashes (sequential)
    std::vector<ReadMetadata> metadata(records.size());
    std::vector<std::vector<uint64_t>> all_hashes(records.size());
    
    for (auto i = 0; i < records.size(); ++i) {
        const auto& record = records[i];
        
        // Extract metadata with early filtering
        metadata[i] = extract_read_metadata(record, min_length, min_quality);
        
        // Only compute hashes for valid reads
        if (metadata[i].is_valid) {
            all_hashes[i] = compute_read_hashes(record, hash_adaptor);
        }
    }
    
    // PASS 2: Process reads in parallel
    #pragma omp parallel for num_threads(num_threads) shared(result)
    for (auto i = 0; i < records.size(); ++i) {
        try {
            if (!metadata[i].is_valid) {
                continue;
            }
            
            process_single_read(
                metadata[i].read_id,
                metadata[i].length,
                metadata[i].mean_quality,
                metadata[i].compression_ratio,
                all_hashes[i],
                agent,
                result,
                records[i],
                is_dehost
            );
            
        } catch (const std::exception& e) {
            PLOG_ERROR << "Error processing read " << records[i].id() << ": " << e.what();
            // Continue processing other reads
        }
    }
    
    total_reads_processed += records.size();
}

/**
 * @brief Process a batch of paired-end reads with optimized algorithm
 * 
 * Similar to process_read_batch but handles paired-end reads by:
 * - Combining metadata from both reads
 * - Processing hashes from both reads
 * - Maintaining proper pairing information
 * 
 * @tparam record_type Type of the sequence records
 * @tparam result_type Type of the result container
 * @tparam hash_adaptor_type Type of the hash adaptor
 * @param records1 Vector of first read pairs
 * @param records2 Vector of second read pairs
 * @param agent IBF membership agent
 * @param hash_adaptor The minimizer hash adaptor
 * @param result Output result container
 * @param min_length Minimum read length threshold
 * @param num_threads Number of OpenMP threads to use
 * @param total_reads_processed Reference to counter for processed read pairs
 */
template<typename record_type, typename result_type, typename hash_adaptor_type, typename agent_type>
void process_paired_read_batch(const std::vector<record_type>& records1,
                              const std::vector<record_type>& records2,
                              agent_type& agent,
                              const hash_adaptor_type& hash_adaptor,
                              result_type& result,
                              uint32_t min_length,
                              int num_threads,
                              uint64_t& total_reads_processed) {
    
    // PASS 1: Extract metadata and compute hashes for both reads
    std::vector<ReadMetadata> metadata1(records1.size());
    std::vector<ReadMetadata> metadata2(records2.size());
    std::vector<std::vector<uint64_t>> hashes1(records1.size());
    std::vector<std::vector<uint64_t>> hashes2(records2.size());
    
    for (auto i = 0; i < records1.size(); ++i) {
        const auto& record1 = records1[i];
        const auto& record2 = records2[i];
        
        // Validate pairing
        auto id1 = record1.id();
        id1.erase(id1.size() - 1);
        auto id2 = record2.id();
        id2.erase(id2.size() - 1);
        if (id1 != id2) {
            std::cout << id1 << " " << id2;
            throw std::runtime_error("Your pairs don't match for read ids.");
        }
        
        // Extract metadata for both reads
        metadata1[i] = extract_read_metadata(record1, min_length, 0.0f);
        metadata2[i] = extract_read_metadata(record2, min_length, 0.0f);
        
        // Compute hashes for valid pairs
        if (metadata1[i].is_valid && metadata2[i].is_valid) {
            hashes1[i] = compute_read_hashes(record1, hash_adaptor);
            hashes2[i] = compute_read_hashes(record2, hash_adaptor);
        }
    }
    
    // PASS 2: Process paired reads in parallel
    #pragma omp parallel for num_threads(num_threads) shared(result)
    for (auto i = 0; i < records1.size(); ++i) {
        try {
            if (!metadata1[i].is_valid || !metadata2[i].is_valid) {
                continue;
            }
            
            const auto& record1 = records1[i];
            const auto& record2 = records2[i];
            
            // Combine metadata for paired reads
            const auto read_id = metadata1[i].read_id;
            const uint32_t read_length = metadata1[i].length + metadata2[i].length;
            
            // Combined mean quality
            const float total_qual = (metadata1[i].mean_quality * metadata1[i].length) +
                                    (metadata2[i].mean_quality * metadata2[i].length);
            const float mean_quality = (read_length > 0) 
                ? total_qual / static_cast<float>(read_length) 
                : 0.0f;
            
            // Combined compression ratio
            const auto combined_seq = sequence_to_string(record1.sequence()) + 
                                     sequence_to_string(record2.sequence());
            const float compression_ratio = get_compression_ratio(combined_seq);
            
            // Create ReadEntry with combined metadata
            auto read = ReadEntry(read_id, read_length, mean_quality, compression_ratio,
                                 result.input_summary());
            
            // Process hashes from both reads
            for (const auto& hash_value : hashes1[i]) {
                const auto &entry = agent.bulk_contains(hash_value);
                read.update_entry(entry);
            }
            
            for (const auto& hash_value : hashes2[i]) {
                const auto &entry = agent.bulk_contains(hash_value);
                read.update_entry(entry);
            }
            
            read.post_process(result.input_summary());
            
            #pragma omp critical(add_read_to_results)
            result.add_paired_read(read, record1, record2);
            
        } catch (const std::exception& e) {
            PLOG_ERROR << "Error processing paired read " << records1[i].id() << ": " << e.what();
        }
    }
    
    total_reads_processed += records1.size();
}

#endif // CHARON_READ_PROCESSOR_H
