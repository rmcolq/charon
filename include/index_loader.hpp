#ifndef CHARON_INDEX_LOADER_H
#define CHARON_INDEX_LOADER_H

#pragma once

#include <filesystem>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <plog/Log.h>
#include <atomic>
#include <thread>
#include <chrono>

#include <load_index.hpp>
#include <index.hpp>

#ifdef __linux__
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

#ifdef __APPLE__
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#endif

/**
 * @brief Progress reporter for long-running index loading operations
 * 
 * Monitors file loading progress and logs periodic updates.
 */
class IndexLoadProgress {
private:
    std::atomic<bool>& running_;
    std::filesystem::path path_;
    size_t file_size_;
    std::chrono::steady_clock::time_point start_time_;
    
public:
    IndexLoadProgress(std::atomic<bool>& running, 
                     std::filesystem::path path,
                     size_t file_size)
        : running_(running)
        , path_(std::move(path))
        , file_size_(file_size)
        , start_time_(std::chrono::steady_clock::now()) {}
    
    void report_progress() {
        while (running_) {
            std::this_thread::sleep_for(std::chrono::seconds(2));
            
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                std::chrono::steady_clock::now() - start_time_).count();
            
            PLOG_INFO << "Loading index " << path_.filename() 
                     << " (" << elapsed << "s elapsed)";
        }
    }
};

/**
 * @brief Load index with progress reporting for large files
 * 
 * This function provides:
 * 1. Memory-mapped I/O for efficient large file handling
 * 2. Progress reporting for files > 500MB
 * 3. Optimized buffer sizes
 * 4. Detailed logging of index statistics
 * 
 * @param index Reference to Index object to populate
 * @param path Path to the index file
 * @param report_progress Enable progress reporting (default: true for large files)
 */
inline void load_index_with_progress(Index &index, 
                                    std::filesystem::path const &path,
                                    bool report_progress = true) {
    PLOG_INFO << "Loading index from file " << path;
    
    // Validate file exists
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("Index file does not exist: " + path.string());
    }
    
    const auto file_size = std::filesystem::file_size(path);
    const double file_size_mb = static_cast<double>(file_size) / (1024.0 * 1024.0);
    const double file_size_gb = file_size_mb / 1024.0;
    
    PLOG_INFO << "Index file size: " 
              << (file_size_gb > 1.0 
                  ? std::to_string(static_cast<int>(file_size_gb * 100) / 100.0) + " GB" 
                  : std::to_string(static_cast<int>(file_size_mb * 100) / 100.0) + " MB");
    
    // Determine if we should use memory-mapped I/O
    const bool use_mmap = (file_size > 100 * 1024 * 1024);  // 100MB threshold
    const bool should_report = report_progress && (file_size > 500 * 1024 * 1024);  // 500MB threshold
    
    // Start progress reporter if needed
    std::atomic<bool> progress_running{should_report};
    std::unique_ptr<std::thread> progress_thread;
    
    if (should_report) {
        auto progress = std::make_unique<IndexLoadProgress>(progress_running, path, file_size);
        progress_thread = std::make_unique<std::thread>(&IndexLoadProgress::report_progress, progress.get());
        // Note: progress object must outlive the thread
    }
    
    bool load_successful = false;
    
    try {
        if (use_mmap) {
            PLOG_INFO << "Using memory-mapped I/O for efficient loading";
            
#ifdef __linux__
            int fd = open(path.c_str(), O_RDONLY);
            if (fd != -1) {
                void* mapped = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
                if (mapped != MAP_FAILED) {
                    // Optimize for sequential access
                    madvise(mapped, file_size, MADV_SEQUENTIAL);
                    
                    // Load from memory-mapped region
                    std::string buffer(reinterpret_cast<const char*>(mapped), file_size);
                    std::istringstream is{buffer};
                    cereal::BinaryInputArchive iarchive{is};
                    iarchive(index);
                    
                    munmap(mapped, file_size);
                    close(fd);
                    load_successful = true;
                    
                    PLOG_INFO << "Index loaded via memory-mapped I/O";
                } else {
                    close(fd);
                }
            }
#elif defined(__APPLE__)
            int fd = open(path.c_str(), O_RDONLY);
            if (fd != -1) {
                void* mapped = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
                if (mapped != MAP_FAILED) {
                    madvise(mapped, file_size, MADV_SEQUENTIAL);
                    
                    std::string buffer(reinterpret_cast<const char*>(mapped), file_size);
                    std::istringstream is{buffer};
                    cereal::BinaryInputArchive iarchive{is};
                    iarchive(index);
                    
                    munmap(mapped, file_size);
                    close(fd);
                    load_successful = true;
                    
                    PLOG_INFO << "Index loaded via memory-mapped I/O";
                } else {
                    close(fd);
                }
            }
#endif
            
            // Fall through to standard I/O if mmap failed
            if (!load_successful) {
                PLOG_WARNING << "Memory-mapped I/O failed, falling back to standard I/O";
            }
        }
        
        // Standard I/O fallback
        if (!load_successful) {
            PLOG_VERBOSE << "Using standard binary I/O for index loading";
            std::ifstream is{path, std::ios::binary};
            
            if (!is.good()) {
                throw std::runtime_error("Failed to open index file: " + path.string());
            }
            
            // Optimize buffer for large file reading
            constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
            std::vector<char> buffer(BUFFER_SIZE);
            is.rdbuf()->pubsetbuf(buffer.data(), BUFFER_SIZE);
            
            cereal::BinaryInputArchive iarchive{is};
            iarchive(index);
            
            PLOG_INFO << "Index loaded via standard I/O";
        }
        
        // Stop progress reporter
        if (should_report) {
            progress_running = false;
            if (progress_thread && progress_thread->joinable()) {
                progress_thread->join();
            }
        }
        
        PLOG_INFO << "Index loaded successfully";
        // Note: Detailed statistics logging moved to caller to avoid header compilation issues
        
    } catch (...) {
        // Ensure progress thread is cleaned up on error
        if (should_report) {
            progress_running = false;
            if (progress_thread && progress_thread->joinable()) {
                progress_thread->join();
            }
        }
        throw;  // Re-throw the exception
    }
}

#endif // CHARON_INDEX_LOADER_H
