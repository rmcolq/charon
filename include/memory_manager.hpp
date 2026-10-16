#ifndef CHARON_MEMORY_MANAGER_H
#define CHARON_MEMORY_MANAGER_H

#pragma once

#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <algorithm>

#include <plog/Log.h>

#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach/mach.h>
#endif
#ifdef __linux__
#include <unistd.h>
#endif

/**
 * @brief Memory budgeting for Charon.
 *
 * Memory use is dominated by three things, and each is accounted for separately:
 *  1. The index, which is resident for the whole run. Its cost is MEASURED (resident memory
 *     before/after loading, or the file size as a lower bound) rather than assumed.
 *  2. The chunk of reads being processed. Chunks are sized in BASES, not reads, so that memory is
 *     flat whatever the read length (a chunk of 1000 reads can be 1 MB or 100 MB).
 *  3. The cache of reads held while the statistical model is fitted. It is capped in BYTES.
 *
 * Whatever is estimated, actual resident memory can be sampled with current_rss() so that the
 * run can report its peak and a caller can react when it approaches the limit.
 */
class MemoryManager {
private:
    size_t total_system_memory_;
    size_t max_allowed_memory_;
    size_t index_bytes_{0};
    size_t peak_rss_{0};

    static constexpr size_t MB = 1024 * 1024;
    static constexpr size_t GB = 1024 * MB;

    // Fixed overhead of the application, buffers, logs, libraries.
    static constexpr size_t BASE_OVERHEAD = 500 * MB;
    // Per-thread working memory (thread-local agent, per-read hash state, compression buffers).
    static constexpr size_t PER_THREAD_OVERHEAD = 64 * MB;
    // A read held in memory costs sequence + quality, one byte each per base.
    static constexpr size_t BYTES_PER_BASE = 2;
    // Container/allocator overhead applied to buffered read data.
    static constexpr double OVERHEAD_FACTOR = 1.5;
    // Fraction of the memory left after index and overhead that is handed out as working memory.
    static constexpr double SAFETY_MARGIN = 0.8;
    // Share of working memory that the model-fitting read cache may use.
    static constexpr double CACHE_SHARE = 0.4;

    // Smallest useful budgets.
    static constexpr size_t MIN_CHUNK_BASES = 1 * MB;
    static constexpr size_t MIN_CACHE_BYTES = 16 * MB;
    // Upper bound on a chunk so that very large limits do not mean very latent output.
    static constexpr size_t MAX_CHUNK_BASES = 256 * MB;

public:
    /**
     * @param max_memory_mb Maximum memory to use in MB (0 = 75% of system memory)
     */
    explicit MemoryManager(size_t max_memory_mb = 0)
        : total_system_memory_(get_system_memory())
        , max_allowed_memory_(calculate_max_memory(max_memory_mb)) {
        PLOG_VERBOSE << "Memory Manager initialized";
        PLOG_VERBOSE << "Total system memory: " << format_size(total_system_memory_);
        PLOG_VERBOSE << "Maximum allowed memory: " << format_size(max_allowed_memory_);
    }

    [[nodiscard]] size_t max_allowed() const { return max_allowed_memory_; }

    [[nodiscard]] size_t total_system_memory() const { return total_system_memory_; }

    [[nodiscard]] size_t index_bytes() const { return index_bytes_; }

    void set_index_bytes(size_t bytes) {
        index_bytes_ = bytes;
        PLOG_INFO << "Index memory: " << format_size(index_bytes_);
    }

    /**
     * @brief Memory needed before any reads are processed: index, fixed overhead and thread state.
     */
    [[nodiscard]] size_t fixed_memory(size_t num_threads) const {
        return index_bytes_ + BASE_OVERHEAD + num_threads * PER_THREAD_OVERHEAD;
    }

    /**
     * @brief Memory left for reads and the model-fitting cache, with the safety margin applied.
     * Returns 0 if the index and fixed overhead already exceed the limit.
     */
    [[nodiscard]] size_t working_memory(size_t num_threads) const {
        const size_t fixed = fixed_memory(num_threads);
        if (fixed >= max_allowed_memory_)
            return 0;
        return static_cast<size_t>((max_allowed_memory_ - fixed) * SAFETY_MARGIN);
    }

    /**
     * @brief Can a run with this index fit at all?
     *
     * Used before loading the index (with the file size as a lower bound) and after loading it
     * (with the measured size), so a run that cannot fit fails with a clear message instead of
     * being killed by the operating system part-way through.
     */
    [[nodiscard]] bool is_feasible(size_t num_threads) const {
        return working_memory(num_threads) >= MIN_CHUNK_BASES * BYTES_PER_BASE + MIN_CACHE_BYTES;
    }

    /**
     * @brief Budget for the reads of one chunk, in bases.
     */
    [[nodiscard]] size_t chunk_base_budget(size_t num_threads) const {
        const size_t chunk_memory = static_cast<size_t>(working_memory(num_threads) * (1.0 - CACHE_SHARE));
        const size_t bases = static_cast<size_t>(chunk_memory / (BYTES_PER_BASE * OVERHEAD_FACTOR));
        return std::clamp(bases, MIN_CHUNK_BASES, MAX_CHUNK_BASES);
    }

    /**
     * @brief Budget in bytes for reads cached while the model is fitted.
     */
    [[nodiscard]] size_t cache_byte_budget(size_t num_threads) const {
        const size_t cache = static_cast<size_t>(working_memory(num_threads) * CACHE_SHARE);
        return std::max(cache, MIN_CACHE_BYTES);
    }

    /**
     * @brief Estimated bytes held for one cached read (sequence and quality of each mate, id, and
     * the per-category hit vectors kept by the read entry).
     */
    [[nodiscard]] static size_t cached_read_bytes(size_t read_length, size_t id_length, size_t num_categories,
                                                  size_t num_hashes, bool paired = false) {
        const size_t mates = paired ? 2 : 1;
        const size_t record = mates * (read_length * BYTES_PER_BASE + id_length);
        const size_t hit_bits = num_categories * ((num_hashes + 7) / 8);
        const size_t entry = 256 + id_length + hit_bits;
        return static_cast<size_t>((record + entry) * OVERHEAD_FACTOR);
    }

    /**
     * @brief Current resident memory of this process in bytes (0 if it cannot be determined).
     */
    [[nodiscard]] static size_t current_rss() {
#ifdef __linux__
        std::ifstream statm("/proc/self/statm");
        size_t size_pages = 0, resident_pages = 0;
        if (statm >> size_pages >> resident_pages) {
            const long page_size = sysconf(_SC_PAGE_SIZE);
            if (page_size > 0)
                return resident_pages * static_cast<size_t>(page_size);
        }
#elif defined(__APPLE__)
        // phys_footprint is what the OS accounts to the process, including pages held in the memory
        // compressor; resident_size alone under-reports on macOS under memory pressure.
        task_vm_info_data_t info;
        mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
        if (task_info(mach_task_self(), TASK_VM_INFO, reinterpret_cast<task_info_t>(&info), &count) ==
            KERN_SUCCESS) {
            return static_cast<size_t>(info.phys_footprint);
        }
#endif
        return 0;
    }

    /**
     * @brief Sample resident memory, remember the peak, and say whether it is within `fraction` of the limit.
     */
    bool within_limit(double fraction = 0.9) {
        const size_t rss = current_rss();
        peak_rss_ = std::max(peak_rss_, rss);
        // If memory cannot be measured, do not claim a problem.
        return rss == 0 || rss <= static_cast<size_t>(max_allowed_memory_ * fraction);
    }

    [[nodiscard]] size_t peak_rss() const { return peak_rss_; }

    [[nodiscard]] std::string get_memory_report(size_t num_threads) const {
        std::string report;
        report += "\n=== Memory Plan ===\n";
        report += "System memory:      " + format_size(total_system_memory_) + "\n";
        report += "Memory limit:       " + format_size(max_allowed_memory_) + "\n";
        report += "Index:              " + format_size(index_bytes_) + "\n";
        report += "Fixed overhead:     " + format_size(fixed_memory(num_threads) - index_bytes_) + "\n";
        report += "Working memory:     " + format_size(working_memory(num_threads)) + "\n";
        report += "Chunk budget:       " + format_size(chunk_base_budget(num_threads) * BYTES_PER_BASE) +
                  " (" + std::to_string(chunk_base_budget(num_threads) / 1000000) + " Mbases)\n";
        report += "Model-fit cache:    " + format_size(cache_byte_budget(num_threads)) + "\n";
        return report;
    }

    static std::string format_size(size_t bytes) {
        const char *units[] = {"B", "KB", "MB", "GB", "TB"};
        int unit_index = 0;
        double size = static_cast<double>(bytes);
        while (size >= 1024.0 && unit_index < 4) {
            size /= 1024.0;
            unit_index++;
        }
        char buffer[64];
        snprintf(buffer, sizeof(buffer), "%.2f %s", size, units[unit_index]);
        return std::string(buffer);
    }

private:
    static size_t get_system_memory() {
#ifdef __linux__
        long pages = sysconf(_SC_PHYS_PAGES);
        long page_size = sysconf(_SC_PAGE_SIZE);
        if (pages > 0 && page_size > 0) {
            size_t total = static_cast<size_t>(pages) * static_cast<size_t>(page_size);
            // Inside a container the limit that matters is the cgroup's, not the host's.
            for (const char *path: {"/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"}) {
                std::ifstream cg(path);
                size_t limit = 0;
                if (cg >> limit && limit > 0)
                    total = std::min(total, limit);
            }
            return total;
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

    size_t calculate_max_memory(size_t user_limit_mb) const {
        if (user_limit_mb > 0) {
            return user_limit_mb * MB;
        }
        // Auto-detect: 75% of system memory (minimum 2GB)
        return std::max(static_cast<size_t>(total_system_memory_ * 0.75), 2 * GB);
    }
};

/**
 * @brief Before loading the index: fail fast if the index file alone cannot fit in the memory limit.
 *
 * The file size is a lower bound on the index's memory. Loading can take minutes, so it is better
 * to refuse a run that cannot fit than to be killed by the operating system after the load.
 * @return false if the run cannot fit
 */
template<typename arguments_type>
bool plan_memory_before_index(MemoryManager &mem_mgr, const arguments_type &opt) {
    std::error_code ec;
    const auto file_bytes = std::filesystem::file_size(opt.db, ec);
    if (ec)
        return true;  // load_index will report the problem
    mem_mgr.set_index_bytes(static_cast<size_t>(file_bytes));
    if (!mem_mgr.is_feasible(opt.threads)) {
        PLOG_ERROR << "Index " << opt.db << " (" << MemoryManager::format_size(file_bytes) << ") does not fit in the memory limit of "
                   << MemoryManager::format_size(mem_mgr.max_allowed()) << " with " << +opt.threads << " thread(s).";
        PLOG_ERROR << "Increase --max-memory (or the memory available to the container), or reduce --threads.";
        return false;
    }
    return true;
}

/**
 * @brief After loading the index: record its measured memory and derive the chunk and cache budgets.
 * @param measured_index_bytes growth in resident memory caused by loading (0 = could not be measured)
 * @return false if the run cannot fit
 */
template<typename arguments_type>
bool plan_memory_after_index(MemoryManager &mem_mgr, arguments_type &opt, size_t measured_index_bytes) {
    // The index file is deserialised into memory, so its size (recorded before loading) is a lower
    // bound; only take the measurement if it is larger, since resident-memory accounting can
    // under-report (e.g. compressed or paged-out memory).
    if (measured_index_bytes > mem_mgr.index_bytes())
        mem_mgr.set_index_bytes(measured_index_bytes);

    if (!mem_mgr.is_feasible(opt.threads)) {
        PLOG_ERROR << "The index uses " << MemoryManager::format_size(mem_mgr.index_bytes())
                   << ", which leaves too little of the " << MemoryManager::format_size(mem_mgr.max_allowed())
                   << " limit for reads with " << +opt.threads << " thread(s).";
        PLOG_ERROR << "Increase --max-memory (or the memory available to the container), or reduce --threads.";
        return false;
    }

    opt.chunk_base_budget = mem_mgr.chunk_base_budget(opt.threads);
    opt.cache_byte_budget = mem_mgr.cache_byte_budget(opt.threads);

    // Chunks are limited by bases; chunk_size only caps the number of reads in a chunk.
    if (opt.chunk_size == 0) {
        opt.chunk_size = 10000;
        PLOG_INFO << "Using up to " << opt.chunk_size << " reads and " << opt.chunk_base_budget / 1000000
                  << " Mbases per chunk";
    }
    PLOG_INFO << mem_mgr.get_memory_report(opt.threads);
    return true;
}

#endif // CHARON_MEMORY_MANAGER_H
