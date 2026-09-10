#include <filesystem>
#include <fstream>
#include <cereal/archives/binary.hpp>
#include <plog/Log.h>

#include <load_index.hpp>

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

void load_index(Index &index, std::filesystem::path const &path) {
    PLOG_INFO << "Loading index from file " << path;
    
    // Check if file exists and get size
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("Index file does not exist: " + path.string());
    }
    
    const auto file_size = std::filesystem::file_size(path);
    const double file_size_mb = static_cast<double>(file_size) / (1024.0 * 1024.0);
    const double file_size_gb = file_size_mb / 1024.0;
    
    PLOG_INFO << "Index file size: " << (file_size_gb > 1.0 
        ? std::to_string(static_cast<int>(file_size_gb * 100) / 100.0) + " GB" 
        : std::to_string(static_cast<int>(file_size_mb * 100) / 100.0) + " MB");
    
    // Use memory-mapped I/O for large files (>100MB) on supported platforms
    // This allows the OS to handle paging efficiently and avoids copying data
    bool use_mmap = false;
    
#ifdef __linux__
    use_mmap = (file_size > 100 * 1024 * 1024);  // 100MB threshold
#elif defined(__APPLE__)
    use_mmap = (file_size > 100 * 1024 * 1024);  // 100MB threshold
#endif
    
    if (use_mmap) {
        PLOG_INFO << "Using memory-mapped I/O for efficient loading";
        
#ifdef __linux__
        int fd = open(path.c_str(), O_RDONLY);
        if (fd == -1) {
            PLOG_WARNING << "Failed to open file for mmap, falling back to standard I/O";
            use_mmap = false;
        } else {
            void* mapped = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
            if (mapped == MAP_FAILED) {
                PLOG_WARNING << "Failed to mmap file, falling back to standard I/O";
                use_mmap = false;
                close(fd);
            } else {
                // Advise the kernel about sequential access pattern
                madvise(mapped, file_size, MADV_SEQUENTIAL);
                
                // Create a stream from the memory-mapped region
                std::string buffer(reinterpret_cast<const char*>(mapped), file_size);
                std::istringstream is{buffer};
                cereal::BinaryInputArchive iarchive{is};
                iarchive(index);
                
                munmap(mapped, file_size);
                close(fd);
                
                PLOG_INFO << "Index loaded via memory-mapped I/O";
            }
        }
#elif defined(__APPLE__)
        int fd = open(path.c_str(), O_RDONLY);
        if (fd == -1) {
            PLOG_WARNING << "Failed to open file for mmap, falling back to standard I/O";
            use_mmap = false;
        } else {
            void* mapped = mmap(nullptr, file_size, PROT_READ, MAP_PRIVATE, fd, 0);
            if (mapped == MAP_FAILED) {
                PLOG_WARNING << "Failed to mmap file, falling back to standard I/O";
                use_mmap = false;
                close(fd);
            } else {
                // Advise the kernel about sequential access pattern
                madvise(mapped, file_size, MADV_SEQUENTIAL);
                
                // Create a stream from the memory-mapped region
                std::string buffer(reinterpret_cast<const char*>(mapped), file_size);
                std::istringstream is{buffer};
                cereal::BinaryInputArchive iarchive{is};
                iarchive(index);
                
                munmap(mapped, file_size);
                close(fd);
                
                PLOG_INFO << "Index loaded via memory-mapped I/O";
            }
        }
#endif
    }
    
    // Fallback to standard I/O
    if (!use_mmap) {
        PLOG_VERBOSE << "Using standard binary I/O for index loading";
        std::ifstream is{path, std::ios::binary};
        
        if (!is.good()) {
            throw std::runtime_error("Failed to open index file: " + path.string());
        }
        
        // Optimize stream buffer for large file reading
        constexpr size_t BUFFER_SIZE = 1024 * 1024;  // 1MB buffer
        std::vector<char> buffer(BUFFER_SIZE);
        is.rdbuf()->pubsetbuf(buffer.data(), BUFFER_SIZE);
        
        cereal::BinaryInputArchive iarchive{is};
        iarchive(index);
        
        PLOG_INFO << "Index loaded via standard I/O";
    }
    
    // Log index statistics
    PLOG_INFO << "Index loaded with " << index.ibf().bin_count() << " bins and " 
              << index.ibf().bit_size() << " bits";
    PLOG_VERBOSE << "Index parameters: k=" << +index.kmer_size() 
                 << ", w=" << +index.window_size() 
                 << ", fpr=" << index.max_fpr();
}
