#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "utils.hpp"
#include "index_arguments.hpp"

#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include <filesystem>
#include <fstream>
#include <random>

using namespace seqan3::literals;

// ─── make_absolute ───────────────────────────────────────────────────────────

TEST_CASE("make_absolute converts relative path to absolute", "[make_absolute]") {
    const std::filesystem::path relative = "test_file.txt";
    const auto result = make_absolute(relative);
    
    REQUIRE(result.is_absolute());
    REQUIRE(result.filename() == "test_file.txt");
}

TEST_CASE("make_absolute handles already absolute paths", "[make_absolute]") {
    const std::filesystem::path absolute = "/tmp/test_file.txt";
    const auto result = make_absolute(absolute);
    
    REQUIRE(result.is_absolute());
    REQUIRE(result == absolute);
}

TEST_CASE("make_absolute handles nested relative paths", "[make_absolute]") {
    const std::filesystem::path relative = "subdir/nested/file.txt";
    const auto result = make_absolute(relative);
    
    REQUIRE(result.is_absolute());
    REQUIRE(result.filename() == "file.txt");
}

// ─── sequence_to_string ──────────────────────────────────────────────────────

TEST_CASE("sequence_to_string converts empty sequence", "[sequence_to_string]") {
    std::vector<seqan3::dna5> seq;
    REQUIRE(sequence_to_string(seq).empty());
}

TEST_CASE("sequence_to_string converts simple sequence", "[sequence_to_string]") {
    std::vector<seqan3::dna5> seq = {'A'_dna5, 'C'_dna5, 'G'_dna5, 'T'_dna5};
    REQUIRE(sequence_to_string(seq) == "ACGT");
}

TEST_CASE("sequence_to_string handles homopolymer", "[sequence_to_string]") {
    std::vector<seqan3::dna5> seq(10, 'A'_dna5);
    REQUIRE(sequence_to_string(seq) == "AAAAAAAAAA");
}

TEST_CASE("sequence_to_string handles N characters", "[sequence_to_string]") {
    std::vector<seqan3::dna5> seq = {'A'_dna5, 'N'_dna5, 'C'_dna5};
    REQUIRE(sequence_to_string(seq) == "ANC");
}

TEST_CASE("sequence_to_string handles long sequence", "[sequence_to_string]") {
    std::vector<seqan3::dna5> seq(1000, 'G'_dna5);
    const auto result = sequence_to_string(seq);
    REQUIRE(result.size() == 1000);
    REQUIRE(result == std::string(1000, 'G'));
}

// ─── store_hashes and load_hashes (roundtrip tests) ─────────────────────────

namespace {
    // Helper to create temporary directory for tests
    std::filesystem::path create_temp_test_dir() {
        const auto temp_dir = std::filesystem::temp_directory_path() / "charon_test_hashes";
        std::filesystem::create_directories(temp_dir);
        return temp_dir;
    }
    
    // Helper to clean up temporary directory
    void cleanup_temp_test_dir(const std::filesystem::path& temp_dir) {
        if (std::filesystem::exists(temp_dir)) {
            std::filesystem::remove_all(temp_dir);
        }
    }
}

TEST_CASE("store_hashes and load_hashes roundtrip with empty set", "[store_hashes][load_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    ankerl::unordered_dense::set<uint64_t> original;
    REQUIRE_NOTHROW(store_hashes("test_empty", original, temp_dir));
    
    auto loaded = load_hashes("test_empty", temp_dir);
    REQUIRE(loaded.empty());
    
    cleanup_temp_test_dir(temp_dir);
}

TEST_CASE("store_hashes and load_hashes roundtrip with single hash", "[store_hashes][load_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    ankerl::unordered_dense::set<uint64_t> original = {42};
    REQUIRE_NOTHROW(store_hashes("test_single", original, temp_dir));
    
    auto loaded = load_hashes("test_single", temp_dir);
    REQUIRE(loaded.size() == 1);
    REQUIRE(loaded[0] == 42);
    
    cleanup_temp_test_dir(temp_dir);
}

TEST_CASE("store_hashes and load_hashes roundtrip with multiple hashes", "[store_hashes][load_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    ankerl::unordered_dense::set<uint64_t> original = {1, 2, 3, 4, 5, 100, 1000, 10000};
    REQUIRE_NOTHROW(store_hashes("test_multi", original, temp_dir));
    
    auto loaded = load_hashes("test_multi", temp_dir);
    REQUIRE(loaded.size() == original.size());
    
    // Check all values are present (order may differ due to set)
    for (const auto& hash : original) {
        REQUIRE(std::find(loaded.begin(), loaded.end(), hash) != loaded.end());
    }
    
    cleanup_temp_test_dir(temp_dir);
}

TEST_CASE("store_hashes and load_hashes roundtrip with large hash set", "[store_hashes][load_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    ankerl::unordered_dense::set<uint64_t> original;
    for (uint64_t i = 0; i < 10000; ++i) {
        original.insert(i * 7);  // Some pattern
    }
    
    REQUIRE_NOTHROW(store_hashes("test_large", original, temp_dir));
    
    auto loaded = load_hashes("test_large", temp_dir);
    REQUIRE(loaded.size() == original.size());
    
    // Spot check some values
    REQUIRE(std::find(loaded.begin(), loaded.end(), 0) != loaded.end());
    REQUIRE(std::find(loaded.begin(), loaded.end(), 7) != loaded.end());
    REQUIRE(std::find(loaded.begin(), loaded.end(), 69993) != loaded.end());  // 9999 * 7
    
    cleanup_temp_test_dir(temp_dir);
}

TEST_CASE("load_hashes throws on non-existent file", "[load_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    REQUIRE_THROWS_AS(load_hashes("nonexistent", temp_dir), std::runtime_error);
    
    cleanup_temp_test_dir(temp_dir);
}

TEST_CASE("store_hashes creates file with correct name", "[store_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    ankerl::unordered_dense::set<uint64_t> hashes = {1, 2, 3};
    REQUIRE_NOTHROW(store_hashes("mytest", hashes, temp_dir));
    
    const auto expected_file = temp_dir / "mytest.min";
    REQUIRE(std::filesystem::exists(expected_file));
    REQUIRE(std::filesystem::file_size(expected_file) > 0);
    
    cleanup_temp_test_dir(temp_dir);
}

// ─── delete_hashes ───────────────────────────────────────────────────────────

TEST_CASE("delete_hashes removes hash files", "[delete_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    // Create some hash files
    ankerl::unordered_dense::set<uint64_t> hashes1 = {1, 2, 3};
    ankerl::unordered_dense::set<uint64_t> hashes2 = {4, 5, 6};
    store_hashes("0", hashes1, temp_dir);
    store_hashes("1", hashes2, temp_dir);
    
    REQUIRE(std::filesystem::exists(temp_dir / "0.min"));
    REQUIRE(std::filesystem::exists(temp_dir / "1.min"));
    
    // Delete them
    std::vector<uint8_t> targets = {0, 1};
    REQUIRE_NOTHROW(delete_hashes(targets, temp_dir.string()));
    
    REQUIRE_FALSE(std::filesystem::exists(temp_dir / "0.min"));
    REQUIRE_FALSE(std::filesystem::exists(temp_dir / "1.min"));
    
    cleanup_temp_test_dir(temp_dir);
}

TEST_CASE("delete_hashes handles non-existent files gracefully", "[delete_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    std::vector<uint8_t> targets = {0, 1, 2};
    REQUIRE_NOTHROW(delete_hashes(targets, temp_dir.string()));
    
    cleanup_temp_test_dir(temp_dir);
}

TEST_CASE("delete_hashes removes empty directory", "[delete_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    // Create and delete a file
    ankerl::unordered_dense::set<uint64_t> hashes = {1, 2, 3};
    store_hashes("0", hashes, temp_dir);
    
    std::vector<uint8_t> targets = {0};
    delete_hashes(targets, temp_dir.string());
    
    // Directory should be removed if empty
    REQUIRE_FALSE(std::filesystem::exists(temp_dir));
}

TEST_CASE("delete_hashes keeps directory with remaining files", "[delete_hashes]") {
    const auto temp_dir = create_temp_test_dir();
    
    // Create two hash files
    ankerl::unordered_dense::set<uint64_t> hashes1 = {1, 2, 3};
    ankerl::unordered_dense::set<uint64_t> hashes2 = {4, 5, 6};
    store_hashes("0", hashes1, temp_dir);
    store_hashes("1", hashes2, temp_dir);
    
    // Delete only one
    std::vector<uint8_t> targets = {0};
    delete_hashes(targets, temp_dir.string());
    
    // Directory should still exist
    REQUIRE(std::filesystem::exists(temp_dir));
    REQUIRE_FALSE(std::filesystem::exists(temp_dir / "0.min"));
    REQUIRE(std::filesystem::exists(temp_dir / "1.min"));
    
    cleanup_temp_test_dir(temp_dir);
}

// ─── bin_size_in_bits ────────────────────────────────────────────────────────

TEST_CASE("bin_size_in_bits returns reasonable value for typical parameters", "[bin_size_in_bits]") {
    IndexArguments opt;
    opt.num_hash = 3;
    opt.max_fpr = 0.01;
    opt.bits = 64;
    
    const uint64_t num_elements = 1000;
    const auto result = bin_size_in_bits(opt, num_elements);
    
    REQUIRE(result > 0);
    REQUIRE(result <= opt.bits);
}

TEST_CASE("bin_size_in_bits increases with more elements", "[bin_size_in_bits]") {
    IndexArguments opt;
    opt.num_hash = 3;
    opt.max_fpr = 0.01;
    opt.bits = 100000;  // Large enough to not cap results
    
    const auto result1 = bin_size_in_bits(opt, 100);
    const auto result2 = bin_size_in_bits(opt, 1000);
    const auto result3 = bin_size_in_bits(opt, 10000);
    
    REQUIRE(result1 < result2);
    REQUIRE(result2 < result3);
}

TEST_CASE("bin_size_in_bits decreases with lower FPR", "[bin_size_in_bits]") {
    IndexArguments opt;
    opt.num_hash = 3;
    opt.bits = 100000;  // Large enough to not cap results
    
    const uint64_t num_elements = 1000;
    const auto result1 = bin_size_in_bits(opt, num_elements);  // max_fpr = 0.01 (default)
    
    opt.max_fpr = 0.001;
    const auto result2 = bin_size_in_bits(opt, num_elements);
    
    opt.max_fpr = 0.0001;
    const auto result3 = bin_size_in_bits(opt, num_elements);
    
    REQUIRE(result1 < result2);
    REQUIRE(result2 < result3);
}

TEST_CASE("bin_size_in_bits caps at opt.bits", "[bin_size_in_bits]") {
    IndexArguments opt;
    opt.num_hash = 3;
    opt.max_fpr = 0.00001;  // Very strict FPR requires many bits
    opt.bits = 32;
    
    const uint64_t num_elements = 100000;
    const auto result = bin_size_in_bits(opt, num_elements);
    
    REQUIRE(result == opt.bits);
}

TEST_CASE("bin_size_in_bits with single hash function", "[bin_size_in_bits]") {
    IndexArguments opt;
    opt.num_hash = 1;
    opt.max_fpr = 0.01;
    opt.bits = 64;
    
    const uint64_t num_elements = 100;
    const auto result = bin_size_in_bits(opt, num_elements);
    
    REQUIRE(result > 0);
    REQUIRE(result <= opt.bits);
}

// ─── max_num_hashes_for_fpr ──────────────────────────────────────────────────

TEST_CASE("max_num_hashes_for_fpr returns positive value", "[max_num_hashes_for_fpr]") {
    IndexArguments opt;
    opt.num_hash = 3;
    opt.max_fpr = 0.01;
    opt.bits = 64;
    
    const auto result = max_num_hashes_for_fpr(opt);
    
    REQUIRE(result > 0);
}

TEST_CASE("max_num_hashes_for_fpr increases with more bits", "[max_num_hashes_for_fpr]") {
    IndexArguments opt;
    opt.num_hash = 3;
    opt.max_fpr = 0.01;
    
    opt.bits = 32;
    const auto result1 = max_num_hashes_for_fpr(opt);
    
    opt.bits = 64;
    const auto result2 = max_num_hashes_for_fpr(opt);
    
    opt.bits = 128;
    const auto result3 = max_num_hashes_for_fpr(opt);
    
    REQUIRE(result1 < result2);
    REQUIRE(result2 < result3);
}

TEST_CASE("max_num_hashes_for_fpr decreases with stricter FPR", "[max_num_hashes_for_fpr]") {
    IndexArguments opt;
    opt.num_hash = 3;
    opt.bits = 64;
    
    opt.max_fpr = 0.1;
    const auto result1 = max_num_hashes_for_fpr(opt);
    
    opt.max_fpr = 0.01;
    const auto result2 = max_num_hashes_for_fpr(opt);
    
    opt.max_fpr = 0.001;
    const auto result3 = max_num_hashes_for_fpr(opt);
    
    REQUIRE(result1 > result2);
    REQUIRE(result2 > result3);
}

TEST_CASE("max_num_hashes_for_fpr with different hash counts", "[max_num_hashes_for_fpr]") {
    IndexArguments opt;
    opt.max_fpr = 0.01;
    opt.bits = 64;
    
    opt.num_hash = 2;
    const auto result1 = max_num_hashes_for_fpr(opt);
    
    opt.num_hash = 3;
    const auto result2 = max_num_hashes_for_fpr(opt);
    
    opt.num_hash = 4;
    const auto result3 = max_num_hashes_for_fpr(opt);
    
    // More hash functions typically allow more elements for same FPR
    REQUIRE(result1 < result2);
    REQUIRE(result2 < result3);
}
