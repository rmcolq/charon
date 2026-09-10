#ifndef CHARON_MEMORY_MANAGER_H
#define CHARON_MEMORY_MANAGER_H

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <algorithm>

#include <plog/Log.h>

/**
 * @brief Memory usage estimator and manager for Charon
 * 
 * Provides utilities to:
 * 1. Estimate memory usage for different chunk sizes
 * 2. Calculate optimal chunk size based on available memory
 * 3. Monitor and enforce memory limits
 */
class MemoryManager {
private:
    size_t total_system_memory_;
    size_t max_allowed_memory_;
    size_t index_memory_estimate_;
    size_t overhead_estimate_;
    
    // Constants for memory estimation
    static constexpr size_t MB = 1024 * 1024;
    static constexpr size_t GB = 1024 * 1024 * 1024;
    
    // Average read length in bases (typical for Nanopore)
    static constexpr size_t AVG_READ_LENGTH = 5000;
    
    // Bytes per base for sequence + quality
    static constexpr size_t BYTES_PER_BASE = 2;  // 1 for seq, 1 for qual
    
    // Overhead factor for C++ objects, vectors, etc.
    static constexpr float OVERHEAD_FACTOR = 1.5f;
    
    // Safety margin (use only 80% of available memory)
    static constexpr float SAFETY_MARGIN = 0.8f;

public:
    /**
     * @brief Construct a new Memory Manager
     * 
     * @param max_memory_mb Maximum memory to use in MB (0 = auto-detect)
     * @param index_size_mb Estimated index size in MB (optional, for better estimates)
     */
    MemoryManager(size_t max_memory_mb = 0, size_t index_size_mb = 0)
        : total_system_memory_(get_system_memory())
        , max_allowed_memory_(calculate_max_memory(max_memory_mb))
        , index_memory_estimate_(index_size_mb * MB)
        , overhead_estimate_(500 * MB)  // Base overhead for application, logs, etc.
    {
        PLOG_INFO << "Memory Manager initialized";
        PLOG_INFO << "Total system memory: " << (total_system_memory_ / GB) << " GB";
        PLOG_INFO << "Maximum allowed memory: " << (max_allowed_memory_ / MB) << " MB";
        
        if (index_size_mb > 0) {
            PLOG_INFO << "Index size estimate: " << (index_size_mb / MB) << " MB";
        }
    }
    
    /**
     * @brief Estimate memory usage for a given chunk size
     * 
     * Memory breakdown:
     * 1. Index (loaded, shared)
     * 2. Read data (chunk_size * avg_read_length * bytes_per_base)
     * 3. Hash storage (chunk_size * avg_hashes_per_read * 8 bytes)
     * 4. Result buffers (chunk_size * result_size)
     * 5. Overhead (threads, buffers, etc.)
     * 
     * @param chunk_size Number of reads to process in parallel
     * @param num_threads Number of OpenMP threads
     * @return Estimated memory usage in bytes
     */
    size_t estimate_memory_usage(size_t chunk_size, size_t num_threads = 1) const {
        // Average number of hashes per read (depends on k-mer size and window)
        constexpr size_t AVG_HASHES_PER_READ = 5000;  // ~1 hash per base
        constexpr size_t HASH_SIZE = 8;  // uint64_t
        
        // Estimate per-read memory
        size_t read_data_size = chunk_size * AVG_READ_LENGTH * BYTES_PER_BASE;
        size_t hash_storage_size = chunk_size * AVG_HASHES_PER_READ * HASH_SIZE;
        size_t result_buffer_size = chunk_size * 1024;  // ~1KB per read for results
        
        // Thread-local storage overhead
        size_t thread_overhead = num_threads * 50 * MB;  // ~50MB per thread
        
        // Total estimate with overhead factor
        size_t variable_memory = static_cast<size_t>(
            (read_data_size + hash_storage_size + result_buffer_size) * OVERHEAD_FACTOR
        );
        
        size_t total = index_memory_estimate_ + overhead_estimate_ + 
                      variable_memory + thread_overhead;
        
        PLOG_VERBOSE << "Memory estimate for chunk_size=" << chunk_size 
                    << ", threads=" << num_threads << ": " 
                    << (total / MB) << " MB";
        PLOG_VERBOSE << "  - Index: " << (index_memory_estimate_ / MB) << " MB";
        PLOG_VERBOSE << "  - Read data: " << (read_data_size / MB) << " MB";
        PLOG_VERBOSE << "  - Hash storage: " << (hash_storage_size / MB) << " MB";
        PLOG_VERBOSE << "  - Result buffers: " << (result_buffer_size / MB) << " MB";
        PLOG_VERBOSE << "  - Thread overhead: " << (thread_overhead / MB) << " MB";
        PLOG_VERBOSE << "  - Base overhead: " << (overhead_estimate_ / MB) << " MB";
        
        return total;
    }
    
    /**
     * @brief Calculate optimal chunk size for given memory constraint
     * 
     * @param num_threads Number of OpenMP threads to use
     * @return Optimal chunk size (number of reads)
     */
    size_t calculate_optimal_chunk_size(size_t num_threads = 1) const {
        // Available memory for variable data
        size_t available_memory = max_allowed_memory_ - 
                                 index_memory_estimate_ - 
                                 overhead_estimate_ - 
                                 (num_threads * 50 * MB);
        
        if (available_memory > max_allowed_memory_) {
            // Overflow happened, not enough memory
            PLOG_WARNING << "Insufficient memory for requested configuration";
            PLOG_WARNING << "Consider reducing threads or index size";
            return 50;  // Minimum safe chunk size
        }
        
        // Calculate chunk size that fits in available memory
        // Using the inverse of estimate_memory_usage formula
        constexpr size_t AVG_HASHES_PER_READ = 5000;
        constexpr size_t HASH_SIZE = 8;
        constexpr size_t AVG_READ_LENGTH = 5000;
        constexpr size_t BYTES_PER_BASE = 2;
        constexpr size_t RESULT_SIZE_PER_READ = 1024;
        
        size_t bytes_per_read = (AVG_READ_LENGTH * BYTES_PER_BASE) +
                               (AVG_HASHES_PER_READ * HASH_SIZE) +
                               RESULT_SIZE_PER_READ;
        bytes_per_read = static_cast<size_t>(bytes_per_read * OVERHEAD_FACTOR);
        
        size_t optimal_chunk = available_memory / bytes_per_read;
        
        // Apply safety margin
        optimal_chunk = static_cast<size_t>(optimal_chunk * SAFETY_MARGIN);
        
        // Enforce reasonable bounds
        optimal_chunk = std::max(optimal_chunk, static_cast<size_t>(50));   // Minimum 50
        optimal_chunk = std::min(optimal_chunk, static_cast<size_t>(10000)); // Maximum 10000
        
        PLOG_INFO << "Optimal chunk size calculated: " << optimal_chunk 
                 << " reads (using " << (max_allowed_memory_ / MB) << " MB limit)";
        
        return optimal_chunk;
    }
    
    /**
     * @brief Check if a configuration is safe for the memory limit
     * 
     * @param chunk_size Proposed chunk size
     * @param num_threads Number of threads
     * @return true if configuration is safe, false otherwise
     */
    bool is_configuration_safe(size_t chunk_size, size_t num_threads = 1) const {
        size_t estimated = estimate_memory_usage(chunk_size, num_threads);
        bool safe = (estimated <= max_allowed_memory_);
        
        if (!safe) {
            PLOG_WARNING << "Configuration may exceed memory limit!";
            PLOG_WARNING << "Estimated: " << (estimated / MB) << " MB, "
                        << "Limit: " << (max_allowed_memory_ / MB) << " MB";
        }
        
        return safe;
    }
    
    /**
     * @brief Get recommended configuration for 16GB systems
     * 
     * @return std::pair<chunk_size, threads> Recommended configuration
     */
    static std::pair<size_t, size_t> get_recommended_16gb_config() {
        // For 16GB systems:
        // - Reserve 4GB for system and other applications
        // - Use up to 12GB for Charon
        // - Index typically 2-4GB uncompressed
        // - Leaves 8-10GB for read processing
        
        unsigned int hw_threads = std::thread::hardware_concurrency();
        size_t threads = std::min(static_cast<size_t>(8), 
                                 static_cast<size_t>(hw_threads > 0 ? hw_threads : 8));
        
        // Conservative chunk size for 16GB
        size_t chunk_size = 500;  // Good balance of performance and memory
        
        PLOG_INFO << "Recommended configuration for 16GB system:";
        PLOG_INFO << "  - Threads: " << threads;
        PLOG_INFO << "  - Chunk size: " << chunk_size;
        PLOG_INFO << "  - Estimated memory: ~8-10 GB";
        
        return {chunk_size, threads};
    }
    
    /**
     * @brief Get human-readable memory report
     * 
     * @param chunk_size Current chunk size
     * @param num_threads Current thread count
     * @return Formatted memory report string
     */
    std::string get_memory_report(size_t chunk_size, size_t num_threads) const {
        size_t estimated = estimate_memory_usage(chunk_size, num_threads);
        size_t available = max_allowed_memory_ - estimated;
        
        std::string report;
        report += "\n=== Memory Usage Report ===\n";
        report += "System Memory:      " + format_size(total_system_memory_) + "\n";
        report += "Memory Limit:       " + format_size(max_allowed_memory_) + "\n";
        report += "Estimated Usage:    " + format_size(estimated) + "\n";
        report += "Available Headroom: " + format_size(available) + "\n";
        report += "Utilization:        " + 
                 std::to_string(static_cast<int>(100.0 * estimated / max_allowed_memory_)) + "%\n";
        
        if (estimated > max_allowed_memory_) {
            report += "\n⚠️  WARNING: Estimated usage exceeds limit!\n";
            report += "Recommendation: Reduce chunk size or thread count\n";
        } else if (available < max_allowed_memory_ * 0.2) {
            report += "\n⚠️  CAUTION: Low memory headroom (<20%)\n";
            report += "Recommendation: Monitor system performance\n";
        } else {
            report += "\n✓ Configuration looks safe\n";
        }
        
        return report;
    }

private:
    /**
     * @brief Get total system memory
     */
    static size_t get_system_memory() {
#ifdef __linux__
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        if (pages > 0 && page_size > 0) {
            return static_cast<size_t>(pages) * static_cast<size_t>(page_size);
        }
#elif defined(__APPLE__)
        int mib[2] = {CTL_HW, HW_MEMSIZE};
        uint64_t size;
        size_t len = sizeof(size);
        if (sysctl(mib, 2, &size, &len, NULL, 0) == 0) {
            return static_cast<size_t>(size);
        }
#endif
        // Fallback: assume 16GB if detection fails
        return 16 * GB;
    }
    
    /**
     * @brief Calculate maximum allowed memory
     */
    size_t calculate_max_memory(size_t user_limit_mb) const {
        if (user_limit_mb > 0) {
            // User specified limit
            return user_limit_mb * MB;
        }
        
        // Auto-detect: use 75% of system memory, but cap at reasonable values
        size_t auto_limit = static_cast<size_t>(total_system_memory_ * 0.75);
        
        // Minimum 2GB, maximum based on system
        auto_limit = std::max(auto_limit, 2 * GB);
        
        return auto_limit;
    }
    
    /**
     * @brief Format size in human-readable format
     */
    static std::string format_size(size_t bytes) {
        const char* units[] = {"B", "KB", "MB", "GB", "TB"};
        int unit_index = 0;
        double size = static_cast<double>(bytes);
        
        while (size >= 1024.0 && unit_index < 4) {
            size /= 1024.0;
            unit_index++;
        }
        
        char buffer[64];
        snprintf(buffer, sizeof(buffer), "%7.2f %s", size, units[unit_index]);
        return std::string(buffer);
    }
};

#endif // CHARON_MEMORY_MANAGER_H
