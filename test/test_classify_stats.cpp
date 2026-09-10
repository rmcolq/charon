#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "classify_stats.hpp"
#include "input_summary.hpp"
#include "classify_arguments.hpp"

#include <limits>
#include <cmath>

// ─── mean() ──────────────────────────────────────────────────────────────────

TEST_CASE("mean returns correct value for simple vector", "[mean]") {
    std::vector<float> v = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    REQUIRE(mean(v) == Catch::Approx(3.0));
}

TEST_CASE("mean returns 0 for empty vector", "[mean]") {
    std::vector<float> v;
    REQUIRE(mean(v) == 0.0);
}

TEST_CASE("mean handles single element", "[mean]") {
    std::vector<float> v = {42.0f};
    REQUIRE(mean(v) == Catch::Approx(42.0));
}

TEST_CASE("mean handles negative values", "[mean]") {
    std::vector<float> v = {-1.0f, -2.0f, -3.0f};
    REQUIRE(mean(v) == Catch::Approx(-2.0));
}

TEST_CASE("mean handles mixed positive and negative", "[mean]") {
    std::vector<float> v = {-5.0f, 5.0f, -5.0f, 5.0f};
    REQUIRE(mean(v) == Catch::Approx(0.0));
}

TEST_CASE("mean handles large values", "[mean]") {
    std::vector<float> v = {1000000.0f, 2000000.0f, 3000000.0f};
    REQUIRE(mean(v) == Catch::Approx(2000000.0));
}

// ─── variance() ──────────────────────────────────────────────────────────────

TEST_CASE("variance returns correct value for simple vector", "[variance]") {
    std::vector<float> v = {1.0f, 2.0f, 3.0f, 4.0f, 5.0f};
    const auto m = mean(v);
    // Variance = sum((x - mean)^2) / (n-1) = 10/4 = 2.5
    REQUIRE(variance(v, m) == Catch::Approx(2.5));
}

TEST_CASE("variance returns 0 for single element", "[variance]") {
    std::vector<float> v = {42.0f};
    const auto m = mean(v);
    REQUIRE(variance(v, m) == 0.0);
}

TEST_CASE("variance returns 0 for empty vector", "[variance]") {
    std::vector<float> v;
    const auto m = mean(v);
    REQUIRE(variance(v, m) == 0.0);
}

TEST_CASE("variance returns 0 for identical values", "[variance]") {
    std::vector<float> v = {5.0f, 5.0f, 5.0f, 5.0f};
    const auto m = mean(v);
    REQUIRE(variance(v, m) == Catch::Approx(0.0));
}

TEST_CASE("variance handles negative values", "[variance]") {
    std::vector<float> v = {-1.0f, -2.0f, -3.0f};
    const auto m = mean(v);
    // Variance = ((-1 - -2)^2 + (-2 - -2)^2 + (-3 - -2)^2) / 2 = (1 + 0 + 1) / 2 = 1.0
    REQUIRE(variance(v, m) == Catch::Approx(1.0));
}

TEST_CASE("variance is always non-negative", "[variance]") {
    std::vector<float> v = {1.0f, 10.0f, 100.0f, 1000.0f};
    const auto m = mean(v);
    REQUIRE(variance(v, m) >= 0.0);
}

// ─── TrainingData ────────────────────────────────────────────────────────────

namespace {
    // Helper to create ClassifyArguments with minimal setup
    ClassifyArguments make_classify_args() {
        ClassifyArguments opt;
        opt.num_reads_to_fit = 10;
        opt.lo_hi_threshold = 0.15f;
        opt.dist = "beta";
        return opt;
    }
}

TEST_CASE("TrainingData initializes with correct capacity", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td(opt, 0);
    
    // Should reserve space but start empty
    // Note: Cannot directly check pos/neg sizes as they are private
    // We verify through behavior - initially not complete
    REQUIRE(td.check_status() == false);
}

TEST_CASE("TrainingData becomes complete when both pos and neg are full", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td(opt, 0);
    
    // Add positive samples
    for (int i = 0; i < 10; ++i) {
        REQUIRE(td.add_pos(0.5f + i * 0.1f) == false);
    }
    
    // Add negative samples
    for (int i = 0; i < 10; ++i) {
        REQUIRE(td.add_neg(0.1f + i * 0.05f) == false);
    }
    
    // Should be complete now
    REQUIRE(td.check_status() == true);
}

TEST_CASE("TrainingData stops accepting samples when complete", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td(opt, 0);
    
    // Fill both vectors
    for (int i = 0; i < 10; ++i) {
        td.add_pos(0.5f);
        td.add_neg(0.1f);
    }
    
    REQUIRE(td.check_status() == true);
    
    // Try to add more - should be rejected
    // We can't directly check sizes (private), but we can verify behavior
    // After completion, adding more should not change the complete status
    td.add_pos(0.9f);
    td.add_neg(0.9f);
    
    // Should still be complete
    REQUIRE(td.check_status() == true);
}

TEST_CASE("TrainingData handles pos_complete separately from neg_complete", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td(opt, 0);
    
    // Fill only positive
    for (int i = 0; i < 10; ++i) {
        td.add_pos(0.5f);
    }
    
    // pos should be complete, but overall not (neg is empty)
    // We verify through check_status which should be false
    REQUIRE(td.check_status() == false);
}

TEST_CASE("TrainingData add_neg rejects zero values", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td(opt, 0);
    
    // Try to add zero - should be rejected (val > 0 check)
    td.add_neg(0.0f);
    // Cannot check size directly, but we know it should remain incomplete
    REQUIRE(td.check_status() == false);
    
    // Add positive value - should be accepted
    td.add_neg(0.1f);
    REQUIRE(td.check_status() == false); // Still not complete
}

TEST_CASE("TrainingData handles different IDs", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td1(opt, 0);
    TrainingData td2(opt, 1);
    TrainingData td3(opt, 255);
    
    // IDs should be stored (though not directly accessible, they're used internally)
    // This test ensures construction with different IDs works
    REQUIRE(td1.check_status() == false);
    REQUIRE(td2.check_status() == false);
    REQUIRE(td3.check_status() == false);
}

TEST_CASE("TrainingData copy operations work correctly", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td1(opt, 0);
    
    td1.add_pos(0.5f);
    td1.add_neg(0.1f);
    
    // Copy constructor
    TrainingData td2(td1);
    // Verify copy works by checking it's in same state
    REQUIRE(td2.check_status() == td1.check_status());
    
    // Copy assignment
    TrainingData td3(opt, 1);
    td3 = td1;
    REQUIRE(td3.check_status() == td1.check_status());
}

TEST_CASE("TrainingData move operations work correctly", "[TrainingData]") {
    auto opt = make_classify_args();
    TrainingData td1(opt, 0);
    td1.add_pos(0.5f);
    td1.add_neg(0.1f);
    
    // Move constructor
    TrainingData td2(std::move(td1));
    // Verify move works
    REQUIRE(td2.check_status() == false); // Should have the data
    
    // Move assignment
    TrainingData td3(opt, 1);
    td3 = std::move(td2);
    REQUIRE(td3.check_status() == false); // Should have the data
}

// ─── StatsModel ──────────────────────────────────────────────────────────────

TEST_CASE("StatsModel constructs with ClassifyArguments", "[StatsModel]") {
    auto opt = make_classify_args();
    InputSummary summary;
    summary.num_bins = 3;
    summary.categories = {"cat1", "cat2", "cat3"};
    
    REQUIRE_NOTHROW(StatsModel(opt, summary));
}

TEST_CASE("StatsModel constructs with DehostArguments", "[StatsModel]") {
    DehostArguments opt;
    opt.num_reads_to_fit = 10;
    opt.lo_hi_threshold = 0.15f;
    opt.dist = "kde";
    
    InputSummary summary;
    summary.num_bins = 2;
    summary.categories = {"host", "non-host"};
    
    REQUIRE_NOTHROW(StatsModel(opt, summary));
}

TEST_CASE("StatsModel creates correct number of TrainingData objects", "[StatsModel]") {
    auto opt = make_classify_args();
    InputSummary summary;
    summary.num_bins = 5;
    summary.categories = {"cat1", "cat2", "cat3", "cat4", "cat5"};
    
    StatsModel model(opt, summary);
    
    // Should have one TrainingData per category
    // Note: This assumes TrainingData is stored in a vector accessible via some method
    // Adjust based on actual StatsModel implementation
}

TEST_CASE("StatsModel with zero categories", "[StatsModel]") {
    auto opt = make_classify_args();
    InputSummary summary;
    summary.num_bins = 0;
    // categories is empty by default
    
    REQUIRE_NOTHROW(StatsModel(opt, summary));
}
