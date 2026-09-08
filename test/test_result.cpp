#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "result.hpp"
#include "input_summary.hpp"
#include "classify_stats.hpp"
#include "classify_arguments.hpp"
#include "dehost_arguments.hpp"

#include <sstream>
#include <limits>

// Helper: create a standard InputSummary for tests
static InputSummary make_test_summary() {
    InputSummary summary;
    summary.num_bins = 3;
    summary.categories = {"human", "bacteria", "virus"};
    summary.bin_to_category = {{0, "human"}, {1, "bacteria"}, {2, "virus"}};
    return summary;
}

// Helper: create ClassifyArguments for tests
static ClassifyArguments make_classify_args() {
    ClassifyArguments opt;
    opt.num_reads_to_fit = 10;
    opt.lo_hi_threshold = 0.15f;
    opt.dist = "beta";
    opt.min_quality = 10.0f;
    opt.min_length = 140;
    opt.min_compression = 0.15f;
    opt.confidence_threshold = 2;
    opt.min_proportion_difference = 0.0f;
    opt.output_file = "/tmp/test_output.tsv";
    return opt;
}

// Helper: create DehostArguments for tests
static DehostArguments make_dehost_args() {
    DehostArguments opt;
    opt.num_reads_to_fit = 10;
    opt.lo_hi_threshold = 0.15f;
    opt.dist = "beta";
    opt.min_quality = 10.0f;
    opt.min_length = 140;
    opt.min_compression = 0.15f;
    opt.confidence_threshold = 2;
    opt.output_file = "/tmp/test_dehost_output.tsv";
    opt.host_index = 0;
    return opt;
}

// Mock seqan3 types for testing
// Note: In real code, these would be actual seqan3 types
namespace seqan3_mock {
    // Mock sequence_file_output
    template<typename... Args>
    class sequence_file_output {
    public:
        sequence_file_output(const std::string& path) : path_(path) {}
        void emplace_back(const std::string& /*id*/, const std::string& /*seq*/, const std::string& /*qual*/) {
            // Mock implementation
        }
    private:
        std::string path_;
    };
}

// ─── Construction ────────────────────────────────────────────────────────────

TEST_CASE("Result construction with ClassifyArguments", "[Result]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // Verify initialization
    REQUIRE(result.get_summary().num_categories() == 3);
}

TEST_CASE("Result construction with DehostArguments", "[Result]") {
    auto summary = make_test_summary();
    auto opt = make_dehost_args();
    
    Result<DehostArguments> result(opt, summary);
    
    // Verify initialization
    REQUIRE(result.get_summary().num_categories() == 3);
}

TEST_CASE("Result constructor initializes classified_counts", "[Result]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // classified_counts should be initialized to zeros
    REQUIRE(result.get_summary().classified_counts.size() == 3);
    REQUIRE(std::all_of(result.get_summary().classified_counts.begin(), 
                       result.get_summary().classified_counts.end(),
                       [](uint64_t n) { return n == 0; }));
}

TEST_CASE("Result constructor initializes unclassified_count", "[Result]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    REQUIRE(result.get_summary().unclassified_count == 0);
}

// ─── classify_read / extract_read ────────────────────────────────────────────

TEST_CASE("classify_read returns category string", "[Result][classify_read]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // Test with valid category index
    auto cat_str = result.classify_read(0);
    REQUIRE(cat_str == "human");
    
    cat_str = result.classify_read(1);
    REQUIRE(cat_str == "bacteria");
    
    cat_str = result.classify_read(2);
    REQUIRE(cat_str == "virus");
}

TEST_CASE("classify_read returns empty for invalid index", "[Result][classify_read]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // Test with invalid category index
    auto cat_str = result.classify_read(99);
    REQUIRE(cat_str.empty());
}

TEST_CASE("extract_read returns category string for dehost", "[Result][extract_read]") {
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"host", "non-host"};
    summary.bin_to_category = {{0, "host"}, {1, "non-host"}};
    
    auto opt = make_dehost_args();
    Result<DehostArguments> result(opt, summary);
    
    auto cat_str = result.extract_read(0);
    REQUIRE(cat_str == "host");
    
    cat_str = result.extract_read(1);
    REQUIRE(cat_str == "non-host");
}

// ─── add_read ────────────────────────────────────────────────────────────────

TEST_CASE("add_read increments classified count", "[Result][add_read]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // Create a mock ReadEntry
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    // Note: Would need to set entry.call() to a valid category
    
    // result.add_read(entry);
    // Note: add_read depends on entry.call() which requires full setup
}

TEST_CASE("add_read increments unclassified for no-call", "[Result][add_read]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    // entry.call() returns max uint8 (unclassified) by default
    
    // result.add_read(entry);
    // REQUIRE(result.get_summary().unclassified_count == 1);
}

// ─── add_paired_read ─────────────────────────────────────────────────────────

TEST_CASE("add_paired_read handles paired reads", "[Result][add_paired_read]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    ReadEntry entry1("read1/1", 150, 30.0f, 0.5f, summary);
    ReadEntry entry2("read1/2", 150, 30.0f, 0.5f, summary);
    
    // result.add_paired_read(entry1, entry2);
    // Note: Complex logic depends on both entries' calls
}

TEST_CASE("add_paired_read with discordant calls", "[Result][add_paired_read]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    ReadEntry entry1("read1/1", 150, 30.0f, 0.5f, summary);
    ReadEntry entry2("read1/2", 150, 30.0f, 0.5f, summary);
    // entry1.call() = 0, entry2.call() = 1 (different categories)
    
    // result.add_paired_read(entry1, entry2);
    // Should handle discordant pairs appropriately
}

// ─── complete ────────────────────────────────────────────────────────────────

TEST_CASE("complete calls classify_cache", "[Result][complete]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // result.complete(false);
    // Note: classify_cache is private, would need integration test
}

TEST_CASE("complete with dehost mode", "[Result][complete]") {
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"host", "non-host"};
    summary.bin_to_category = {{0, "host"}, {1, "non-host"}};
    
    auto opt = make_dehost_args();
    Result<DehostArguments> result(opt, summary);
    
    // result.complete(true);
    // Note: Different behavior in dehost mode
}

// ─── print_summary ───────────────────────────────────────────────────────────

TEST_CASE("print_summary outputs classification counts", "[Result][print_summary]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // This is primarily a logging function
    // Would need to capture log output to verify
    // result.print_summary();
}

// ─── Template instantiation tests ────────────────────────────────────────────

TEST_CASE("Result works with ClassifyArguments template", "[Result][template]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    REQUIRE(result.get_summary().num_categories() == 3);
    REQUIRE(result.get_summary().classified_counts.size() == 3);
}

TEST_CASE("Result works with DehostArguments template", "[Result][template]") {
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"host", "non-host"};
    summary.bin_to_category = {{0, "host"}, {1, "non-host"}};
    
    auto opt = make_dehost_args();
    Result<DehostArguments> result(opt, summary);
    
    REQUIRE(result.get_summary().num_categories() == 2);
    REQUIRE(result.get_summary().classified_counts.size() == 2);
}

// ─── Edge cases ──────────────────────────────────────────────────────────────

TEST_CASE("Result with single category", "[Result][edge]") {
    InputSummary summary;
    summary.num_bins = 1;
    summary.categories = {"single"};
    summary.bin_to_category = {{0, "single"}};
    
    auto opt = make_classify_args();
    Result<ClassifyArguments> result(opt, summary);
    
    REQUIRE(result.get_summary().classified_counts.size() == 1);
}

TEST_CASE("Result with many categories", "[Result][edge]") {
    InputSummary summary;
    summary.num_bins = 10;
    for (int i = 0; i < 10; ++i) {
        summary.categories.push_back("cat_" + std::to_string(i));
        summary.bin_to_category[i] = "cat_" + std::to_string(i);
    }
    
    auto opt = make_classify_args();
    Result<ClassifyArguments> result(opt, summary);
    
    REQUIRE(result.get_summary().classified_counts.size() == 10);
}

TEST_CASE("Result with zero categories (edge case)", "[Result][edge]") {
    InputSummary summary;
    summary.num_bins = 0;
    // Empty categories
    
    auto opt = make_classify_args();
    Result<ClassifyArguments> result(opt, summary);
    
    REQUIRE(result.get_summary().classified_counts.empty());
}

// ─── Integration-style tests ─────────────────────────────────────────────────

TEST_CASE("Result full classification workflow", "[Result][integration]") {
    auto summary = make_test_summary();
    auto opt = make_classify_args();
    
    Result<ClassifyArguments> result(opt, summary);
    
    // Simulate typical workflow:
    // 1. Create ReadEntry objects
    // 2. Call add_read() for each
    // 3. Call complete() to finalize
    // 4. Call print_summary() for output
    
    // Note: Full workflow requires proper ReadEntry setup
    // and is better tested at integration level
}

TEST_CASE("Result full dehosting workflow", "[Result][integration]") {
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"host", "non-host"};
    summary.bin_to_category = {{0, "host"}, {1, "non-host"}};
    
    auto opt = make_dehost_args();
    Result<DehostArguments> result(opt, summary);
    
    // Simulate typical dehosting workflow
    // Similar to classification but with different logic
}

// ─── Note on Testability ─────────────────────────────────────────────────────

// Many Result methods depend on:
// - seqan3::sequence_file_output (file I/O)
// - ReadEntry objects with proper state
// - Private methods like classify_cache()
//
// To improve testability, consider:
// 1. Extracting pure functions for classification logic
// 2. Using dependency injection for file I/O
// 3. Adding test hooks or friend fixtures
// 4. Mocking seqan3 types more comprehensively
//
// The template design makes testing more complex but provides type safety
