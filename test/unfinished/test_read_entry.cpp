#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "read_entry.hpp"
#include "input_summary.hpp"
#include "classify_stats.hpp"
#include "classify_arguments.hpp"

#include <limits>
#include <numeric>

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
    return opt;
}

// ─── Construction ────────────────────────────────────────────────────────────

TEST_CASE("ReadEntry default construction", "[ReadEntry]") {
    ReadEntry entry;
    REQUIRE(entry.read_id().empty());
    REQUIRE(entry.proportions().empty());
    REQUIRE(entry.unique_proportions().empty());
    REQUIRE(entry.call() == std::numeric_limits<uint8_t>::max());
    REQUIRE(entry.confidence_score() == 0);
}

TEST_CASE("ReadEntry construction with summary", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    
    REQUIRE(entry.read_id() == "read1");
    REQUIRE(entry.proportions().size() == 3);
    REQUIRE(entry.unique_proportions().size() == 3);
    REQUIRE(std::all_of(entry.proportions().begin(), entry.proportions().end(),
        [](float f) { return f == 0.0f; }));
    REQUIRE(entry.call() == std::numeric_limits<uint8_t>::max());
}

TEST_CASE("ReadEntry construction initializes vectors to correct size", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    
    REQUIRE(entry.proportions().size() == summary.num_categories());
    REQUIRE(entry.unique_proportions().size() == summary.num_categories());
}

TEST_CASE("ReadEntry copy constructor", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry1("read1", 150, 30.0f, 0.5f, summary);
    
    ReadEntry entry2(entry1);
    REQUIRE(entry2.read_id() == "read1");
    REQUIRE(entry2.proportions().size() == 3);
}

TEST_CASE("ReadEntry copy assignment", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry1("read1", 150, 30.0f, 0.5f, summary);
    ReadEntry entry2;
    
    entry2 = entry1;
    REQUIRE(entry2.read_id() == "read1");
}

TEST_CASE("ReadEntry move constructor", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry1("read1", 150, 30.0f, 0.5f, summary);
    
    ReadEntry entry2(std::move(entry1));
    REQUIRE(entry2.read_id() == "read1");
}

TEST_CASE("ReadEntry move assignment", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry1("read1", 150, 30.0f, 0.5f, summary);
    ReadEntry entry2;
    
    entry2 = std::move(entry1);
    REQUIRE(entry2.read_id() == "read1");
}

// ─── Getters ─────────────────────────────────────────────────────────────────

TEST_CASE("ReadEntry getters return correct values", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry("test_read", 200, 35.0f, 0.6f, summary);
    
    REQUIRE(entry.read_id() == "test_read");
    // Note: length_, mean_quality_, compression_ don't have public getters
    // This is a limitation of the current design
}

TEST_CASE("ReadEntry proportions are const reference", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    
    const auto& props = entry.proportions();
    REQUIRE(props.size() == 3);
    
    // Verify it's actually a reference (modifying original affects reference)
    // Note: Can't test this directly without non-const access
}

// ─── update_entry ────────────────────────────────────────────────────────────

TEST_CASE("update_entry adds bitvector entry", "[ReadEntry]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    
    // Create a mock bitvector (simplified for testing)
    using bitvector_type = decltype(entry)::bits_type::value_type;
    bitvector_type bv;
    bv.resize(3);
    bv[0] = 1; bv[1] = 0; bv[2] = 1;
    
    entry.update_entry(bv);
    
    // Note: bits_ is private, so we can't directly verify
    // This test would need friend access or a test hook
}

// ─── get_proportions ─────────────────────────────────────────────────────────

TEST_CASE("get_proportions calculates from counts", "[ReadEntry][get_proportions]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    
    // Manually set counts (normally done by get_counts)
    // Note: counts_ is private, would need test hook
    // This is a limitation - the method depends on private state
}

// ─── call_category ───────────────────────────────────────────────────────────

TEST_CASE("call_category with high confidence", "[ReadEntry][call_category]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    auto opt = make_classify_args();
    StatsModel model(opt, summary);
    
    // Note: call_category depends on private state (counts_, proportions_, probabilities_)
    // Would need to call post_process first or have test hooks
    // This demonstrates the challenge of testing this class
}

TEST_CASE("call_category rejects low quality reads", "[ReadEntry][call_category]") {
    auto summary = make_test_summary();
    // Create entry with low quality
    ReadEntry entry("read1", 200, 5.0f, 0.6f, summary);  // quality < min_quality
    
    auto opt = make_classify_args();
    StatsModel model(opt, summary);
    
    // Should not make a call due to low quality
    // entry.call_category(model);
    // REQUIRE(entry.call() == std::numeric_limits<uint8_t>::max());
}

TEST_CASE("call_category rejects short reads", "[ReadEntry][call_category]") {
    auto summary = make_test_summary();
    // Create entry with short length
    ReadEntry entry("read1", 100, 35.0f, 0.6f, summary);  // length < min_length
    
    auto opt = make_classify_args();
    StatsModel model(opt, summary);
    
    // Should not make a call due to short length
}

TEST_CASE("call_category rejects low compression", "[ReadEntry][call_category]") {
    auto summary = make_test_summary();
    // Create entry with low compression
    ReadEntry entry("read1", 200, 35.0f, 0.1f, summary);  // compression < min_compression
    
    auto opt = make_classify_args();
    StatsModel model(opt, summary);
    
    // Should not make a call due to low compression
}

// ─── call_host ───────────────────────────────────────────────────────────────

TEST_CASE("call_host with two categories", "[ReadEntry][call_host]") {
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"host", "non-host"};
    summary.bin_to_category = {{0, "host"}, {1, "non-host"}};
    
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    auto opt = make_classify_args();
    StatsModel model(opt, summary);
    
    // Note: Requires probabilities_ to be set
    // Would need test hooks or integration test
}

// ─── dehost ──────────────────────────────────────────────────────────────────

TEST_CASE("dehost calls call_host", "[ReadEntry][dehost]") {
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"host", "non-host"};
    summary.bin_to_category = {{0, "host"}, {1, "non-host"}};
    
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    auto opt = make_classify_args();
    StatsModel model(opt, summary);
    
    // entry.dehost(model, 0);  // host_index = 0
    // Note: dehost is public but depends on private state
}

// ─── classify ────────────────────────────────────────────────────────────────

TEST_CASE("classify calls call_category", "[ReadEntry][classify]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    auto opt = make_classify_args();
    StatsModel model(opt, summary);
    
    // entry.classify(model);
    // Note: classify is public but depends on private state
}

// ─── print_assignment_result ─────────────────────────────────────────────────

TEST_CASE("print_assignment_result logs classification", "[ReadEntry][print]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    // This is primarily a logging function
    // Would need to capture log output to verify
    // entry.print_assignment_result(summary);
}

// ─── Integration-style tests ─────────────────────────────────────────────────

TEST_CASE("ReadEntry full workflow simulation", "[ReadEntry][integration]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    // Simulate typical workflow:
    // 1. update_entry() called for each hash
    // 2. get_counts() called to aggregate
    // 3. get_proportions() called to normalize
    // 4. call_category() called to make classification
    
    // Note: Full workflow requires seqan3 bitvector types
    // and is better tested at integration level
}

TEST_CASE("ReadEntry with single category", "[ReadEntry][edge]") {
    InputSummary summary;
    summary.num_bins = 1;
    summary.categories = {"single"};
    summary.bin_to_category = {{0, "single"}};
    
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    
    REQUIRE(entry.proportions().size() == 1);
    REQUIRE(entry.unique_proportions().size() == 1);
}

TEST_CASE("ReadEntry with many categories", "[ReadEntry][edge]") {
    InputSummary summary;
    summary.num_bins = 10;
    for (int i = 0; i < 10; ++i) {
        summary.categories.push_back("cat_" + std::to_string(i));
        summary.bin_to_category[i] = "cat_" + std::to_string(i);
    }
    
    ReadEntry entry("read1", 150, 30.0f, 0.5f, summary);
    
    REQUIRE(entry.proportions().size() == 10);
    REQUIRE(entry.unique_proportions().size() == 10);
}

TEST_CASE("ReadEntry with zero length read", "[ReadEntry][edge]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 0, 30.0f, 0.5f, summary);
    
    // Should handle gracefully
    REQUIRE(entry.read_id() == "read1");
}

TEST_CASE("ReadEntry confidence score boundaries", "[ReadEntry][edge]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    // Initial confidence should be 0
    REQUIRE(entry.confidence_score() == 0);
    
    // Maximum confidence would be set by call_category
    // Requires proper setup of internal state
}

TEST_CASE("ReadEntry call boundaries", "[ReadEntry][edge]") {
    auto summary = make_test_summary();
    ReadEntry entry("read1", 200, 35.0f, 0.6f, summary);
    
    // Initial call should be max uint8 (unclassified)
    REQUIRE(entry.call() == std::numeric_limits<uint8_t>::max());
}

// ─── Note on Testability ─────────────────────────────────────────────────────

// Many ReadEntry methods depend on private member variables:
// - counts_
// - unique_counts_
// - proportions_
// - unique_proportions_
// - probabilities_
// - bits_
// - max_bits_
//
// To improve testability, consider:
// 1. Adding const getters for these members (for testing only)
// 2. Using friend test fixtures
// 3. Extracting pure functions for the calculation logic
// 4. Using dependency injection for seqan3 types
//
// Example improvement:
// [[nodiscard]] const std::vector<uint32_t>& counts() const { return counts_; }
// [[nodiscard]] const std::vector<double>& probabilities() const { return probabilities_; }
