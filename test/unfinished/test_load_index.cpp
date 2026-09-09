#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>

#include "load_index.hpp"
#include "index.hpp"
#include "store_index.hpp"

namespace fs = std::filesystem;

// Test data directory
static const std::string TEST_DATA_DIR = "/Users/rmcolq/Work/git/charon/build";

// Helper to create a temporary test index
static void create_test_index(const fs::path& path) {
    // Create a minimal index using the same approach as test_index
    IndexArguments args;
    args.window_size = 41;
    args.kmer_size = 19;
    args.max_fpr = 0.01;
    
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"bin1", "bin2"};
    summary.bin_to_category[0] = "bin1";
    summary.bin_to_category[1] = "bin2";
    
    InputStats stats;
    stats.num_files = 2;
    
    // Create a minimal IBF
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf_uncompressed{
        seqan3::bin_count{2u},
        seqan3::bin_size{128u},
        seqan3::hash_function_count{2u}
    };
    seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> ibf{ibf_uncompressed};
    
    Index index{args, summary, stats, ibf};
    
    // Store the index
    store_index(path, std::move(index));
}

// Helper to clean up temporary files
static void cleanup_test_files(const std::vector<fs::path>& files) {
    for (const auto& file : files) {
        if (fs::exists(file)) {
            fs::remove(file);
        }
    }
}

// ─── Basic Load Tests ────────────────────────────────────────────────────────

TEST_CASE("load_index loads valid index file", "[LoadIndex][Basic]") {
    fs::path test_index_path = fs::path(TEST_DATA_DIR) / "test_index_file.bin";
    
    // Create a test index
    create_test_index(test_index_path);
    REQUIRE(fs::exists(test_index_path));
    
    // Load it back
    Index index;
    REQUIRE_NOTHROW(load_index(index, test_index_path));
    
    // Verify the loaded index has correct properties
    REQUIRE(index.summary().num_bins == 2);
    REQUIRE(index.summary().categories.size() == 2);
    REQUIRE(index.summary().bin_to_category.at(0) == "bin1");
    REQUIRE(index.summary().bin_to_category.at(1) == "bin2");
    
    // Cleanup
    cleanup_test_files({test_index_path});
}

TEST_CASE("load_index throws on non-existent file", "[LoadIndex][ErrorHandling]") {
    fs::path non_existent_path = fs::path(TEST_DATA_DIR) / "non_existent_index.bin";
    
    // Make sure file doesn't exist
    if (fs::exists(non_existent_path)) {
        fs::remove(non_existent_path);
    }
    
    Index index;
    
    // Should throw an exception
    REQUIRE_THROWS_AS(load_index(index, non_existent_path), std::exception);
}

TEST_CASE("load_index throws on corrupted file", "[LoadIndex][ErrorHandling]") {
    fs::path corrupted_path = fs::path(TEST_DATA_DIR) / "corrupted_index.bin";
    
    // Create a corrupted file (just random bytes)
    {
        std::ofstream ofs{corrupted_path, std::ios::binary};
        std::string garbage = "This is not a valid index file";
        ofs.write(garbage.c_str(), garbage.size());
    }
    
    REQUIRE(fs::exists(corrupted_path));
    
    Index index;
    
    // Should throw an exception when trying to deserialize
    REQUIRE_THROWS_AS(load_index(index, corrupted_path), std::exception);
    
    // Cleanup
    cleanup_test_files({corrupted_path});
}

// ─── Index Content Verification ──────────────────────────────────────────────

TEST_CASE("load_index preserves index properties", "[LoadIndex][Verification]") {
    fs::path test_index_path = fs::path(TEST_DATA_DIR) / "test_index_verify.bin";
    
    // Create a test index with specific properties
    create_test_index(test_index_path);
    
    // Load it back
    Index loaded_index;
    load_index(loaded_index, test_index_path);
    
    // Verify all properties are preserved
    REQUIRE(loaded_index.summary().num_bins == 2);
    REQUIRE(loaded_index.stats().num_files == 2);
    REQUIRE(loaded_index.ibf().bin_count() == 2);
    
    // Cleanup
    cleanup_test_files({test_index_path});
}

TEST_CASE("load_index handles multiple loads", "[LoadIndex][Stress]") {
    fs::path test_index_path = fs::path(TEST_DATA_DIR) / "test_index_multi.bin";
    
    // Create a test index
    create_test_index(test_index_path);
    
    // Load it multiple times
    for (int i = 0; i < 3; ++i) {
        Index index;
        REQUIRE_NOTHROW(load_index(index, test_index_path));
        
        // Verify each load works
        REQUIRE(index.summary().num_bins == 2);
    }
    
    // Cleanup
    cleanup_test_files({test_index_path});
}

// ─── Path Handling ───────────────────────────────────────────────────────────

TEST_CASE("load_index handles relative paths", "[LoadIndex][PathHandling]") {
    // Save current directory
    fs::path original_dir = fs::current_path();
    
    try {
        // Change to build directory
        fs::current_path(TEST_DATA_DIR);
        
        // Create index with relative path
        fs::path relative_path = "test_index_relative.bin";
        create_test_index(relative_path);
        REQUIRE(fs::exists(relative_path));
        
        // Load with relative path
        Index index;
        REQUIRE_NOTHROW(load_index(index, relative_path));
        REQUIRE(index.summary().num_bins == 2);
        
        // Cleanup
        cleanup_test_files({relative_path});
        
    } catch (...) {
        // Restore original directory on error
        fs::current_path(original_dir);
        throw;
    }
    
    // Restore original directory
    fs::current_path(original_dir);
}

TEST_CASE("load_index handles absolute paths", "[LoadIndex][PathHandling]") {
    fs::path absolute_path = fs::absolute(fs::path(TEST_DATA_DIR) / "test_index_absolute.bin");
    
    // Create a test index
    create_test_index(absolute_path);
    REQUIRE(fs::exists(absolute_path));
    
    // Load with absolute path
    Index index;
    REQUIRE_NOTHROW(load_index(index, absolute_path));
    REQUIRE(index.summary().num_bins == 2);
    
    // Cleanup
    cleanup_test_files({absolute_path});
}

// ─── Edge Cases ──────────────────────────────────────────────────────────────

TEST_CASE("load_index handles empty directory path", "[LoadIndex][EdgeCases]") {
    // Try to load from a directory path (should fail)
    fs::path dir_path = fs::path(TEST_DATA_DIR);
    
    Index index;
    
    // Should throw because it's a directory, not a file
    REQUIRE_THROWS_AS(load_index(index, dir_path), std::exception);
}

TEST_CASE("load_index handles path with special characters", "[LoadIndex][EdgeCases]") {
    fs::path special_path = fs::path(TEST_DATA_DIR) / "test_index_with spaces & symbols.bin";
    
    // Create a test index
    create_test_index(special_path);
    REQUIRE(fs::exists(special_path));
    
    // Load it back
    Index index;
    REQUIRE_NOTHROW(load_index(index, special_path));
    REQUIRE(index.summary().num_bins == 2);
    
    // Cleanup
    cleanup_test_files({special_path});
}
