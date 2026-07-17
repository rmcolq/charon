#include <catch2/catch_test_macros.hpp>

#include "input_summary.hpp"

// Helper: build a standard InputSummary for reuse across tests
static InputSummary make_summary() {
    InputSummary s;
    s.num_bins = 3;
    s.categories = {"human", "bacteria", "virus"};
    s.bin_to_category = {{0, "human"}, {1, "bacteria"}, {2, "virus"}};
    return s;
}

// ─── num_categories ──────────────────────────────────────────────────────────

TEST_CASE("num_categories returns correct count", "[InputSummary][num_categories]") {
    REQUIRE(make_summary().num_categories() == 3);
}

TEST_CASE("num_categories returns zero for empty summary", "[InputSummary][num_categories]") {
    InputSummary s;
    REQUIRE(s.num_categories() == 0);
}

// ─── category_index ──────────────────────────────────────────────────────────

TEST_CASE("category_index returns correct index for known category", "[InputSummary][category_index]") {
    auto s = make_summary();
    REQUIRE(s.category_index("human") == 0);
    REQUIRE(s.category_index("bacteria") == 1);
    REQUIRE(s.category_index("virus") == 2);
}

TEST_CASE("category_index returns max uint8 for unknown category", "[InputSummary][category_index]") {
    auto s = make_summary();
    REQUIRE(s.category_index("fungus") == std::numeric_limits<uint8_t>::max());
}

TEST_CASE("category_index returns max uint8 for empty categories", "[InputSummary][category_index]") {
    InputSummary s;
    REQUIRE(s.category_index("human") == std::numeric_limits<uint8_t>::max());
}

// ─── host_category_index ─────────────────────────────────────────────────────

TEST_CASE("host_category_index finds 'human' category", "[InputSummary][host_category_index]") {
    auto s = make_summary(); // categories = {"human", "bacteria", "virus"}
    REQUIRE(s.host_category_index() == 0);
}

TEST_CASE("host_category_index finds 'host' category", "[InputSummary][host_category_index]") {
    InputSummary s;
    s.categories = {"bacteria", "host", "virus"};
    REQUIRE(s.host_category_index() == 1);
}

TEST_CASE("host_category_index returns minimum index when both human and host present",
          "[InputSummary][host_category_index]") {
    InputSummary s;
    s.categories = {"host", "human", "virus"}; // host=0, human=1 → min=0
    REQUIRE(s.host_category_index() == 0);
}

// ─── category_name ───────────────────────────────────────────────────────────

TEST_CASE("category_name returns correct name for valid index", "[InputSummary][category_name]") {
    auto s = make_summary();
    REQUIRE(s.category_name(0) == "human");
    REQUIRE(s.category_name(1) == "bacteria");
    REQUIRE(s.category_name(2) == "virus");
}

TEST_CASE("category_name returns empty string for out-of-bounds index",
          "[InputSummary][category_name]") {
    auto s = make_summary(); // 3 categories, valid indices 0-2
    REQUIRE(s.category_name(3) == "");  // was a bug: > vs >= meant this called .at(3) and threw
    REQUIRE(s.category_name(255) == "");
}

TEST_CASE("category_name returns empty string for empty categories",
          "[InputSummary][category_name]") {
    InputSummary s;
    REQUIRE(s.category_name(0) == "");
}
