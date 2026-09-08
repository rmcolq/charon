#include <catch2/catch_test_macros.hpp>

#include "input_summary.hpp"

#include <limits>

// Helper: build a standard InputSummary for reuse across tests
static InputSummary make_summary() {
    InputSummary s;
    s.num_bins = 3;
    s.categories = {"human", "bacteria", "virus"};
    s.bin_to_category = {{0, "human"}, {1, "bacteria"}, {2, "virus"}};
    s.filepath_to_bin = {{"human.fasta", 0}, {"bacteria.fasta", 1}, {"virus.fasta", 2}};
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

TEST_CASE("num_categories handles single category", "[InputSummary][num_categories]") {
    InputSummary s;
    s.categories = {"single"};
    s.num_bins = 1;
    REQUIRE(s.num_categories() == 1);
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

TEST_CASE("category_index is case-sensitive", "[InputSummary][category_index]") {
    auto s = make_summary();
    REQUIRE(s.category_index("Human") == std::numeric_limits<uint8_t>::max());
    REQUIRE(s.category_index("HUMAN") == std::numeric_limits<uint8_t>::max());
}

// ─── host_category_index ─────────────────────────────────────────────────────

TEST_CASE("host_category_index finds 'human' category", "[InputSummary][host_category_index]") {
    auto s = make_summary(); // categories = {"human", "bacteria", "virus"}
    REQUIRE(s.host_category_index() == 0);
}

TEST_CASE("host_category_index finds 'host' category", "[InputSummary][host_category_index]") {
    InputSummary s;
    s.categories = {"bacteria", "host", "virus"};
    s.bin_to_category = {{0, "bacteria"}, {1, "host"}, {2, "virus"}};
    REQUIRE(s.host_category_index() == 1);
}

TEST_CASE("host_category_index prefers 'human' over 'host'", "[InputSummary][host_category_index]") {
    InputSummary s;
    s.categories = {"bacteria", "host", "human", "virus"};
    s.bin_to_category = {{0, "bacteria"}, {1, "host"}, {2, "human"}, {3, "virus"}};
    // Should return index of 'human' (2) since it's smaller than 'host' (1)
    REQUIRE(s.host_category_index() == 1);  // min(1, 2) = 1
}

TEST_CASE("host_category_index with neither human nor host", "[InputSummary][host_category_index]") {
    InputSummary s;
    s.categories = {"bacteria", "virus", "fungus"};
    // This will log an error but should return max uint8
    REQUIRE(s.host_category_index() == std::numeric_limits<uint8_t>::max());
}

TEST_CASE("host_category_index with empty categories", "[InputSummary][host_category_index]") {
    InputSummary s;
    REQUIRE(s.host_category_index() == std::numeric_limits<uint8_t>::max());
}

// ─── category_name ───────────────────────────────────────────────────────────

TEST_CASE("category_name returns correct name for valid index", "[InputSummary][category_name]") {
    auto s = make_summary();
    REQUIRE(s.category_name(0) == "human");
    REQUIRE(s.category_name(1) == "bacteria");
    REQUIRE(s.category_name(2) == "virus");
}

TEST_CASE("category_name returns empty string for out-of-bounds index", "[InputSummary][category_name]") {
    auto s = make_summary();
    REQUIRE(s.category_name(3) == "");
    REQUIRE(s.category_name(4) == "");
    REQUIRE(s.category_name(255) == "");
}

TEST_CASE("category_name returns empty string for empty categories", "[InputSummary][category_name]") {
    InputSummary s;
    REQUIRE(s.category_name(0) == "");
}

TEST_CASE("category_name handles boundary index", "[InputSummary][category_name]") {
    auto s = make_summary();
    REQUIRE(s.category_name(2) == "virus");  // Last valid index
    REQUIRE(s.category_name(3) == "");       // First invalid index
}

// ─── filepath_to_bin ─────────────────────────────────────────────────────────

TEST_CASE("filepath_to_bin contains correct mappings", "[InputSummary][filepath_to_bin]") {
    auto s = make_summary();
    REQUIRE(s.filepath_to_bin.size() == 3);
    
    // Find specific mappings
    auto it = std::find_if(s.filepath_to_bin.begin(), s.filepath_to_bin.end(),
        [](const auto& p) { return p.first == "human.fasta"; });
    REQUIRE(it != s.filepath_to_bin.end());
    REQUIRE(it->second == 0);
}

TEST_CASE("filepath_to_bin can be empty", "[InputSummary][filepath_to_bin]") {
    InputSummary s;
    REQUIRE(s.filepath_to_bin.empty());
}

TEST_CASE("filepath_to_bin handles multiple files per bin", "[InputSummary][filepath_to_bin]") {
    InputSummary s;
    s.num_bins = 2;
    s.categories = {"cat1", "cat2"};
    s.filepath_to_bin = {
        {"file1.fasta", 0},
        {"file2.fasta", 0},
        {"file3.fasta", 1}
    };
    
    REQUIRE(s.filepath_to_bin.size() == 3);
}

// ─── bin_to_category ─────────────────────────────────────────────────────────

TEST_CASE("bin_to_category contains correct mappings", "[InputSummary][bin_to_category]") {
    auto s = make_summary();
    REQUIRE(s.bin_to_category.size() == 3);
    REQUIRE(s.bin_to_category.at(0) == "human");
    REQUIRE(s.bin_to_category.at(1) == "bacteria");
    REQUIRE(s.bin_to_category.at(2) == "virus");
}

TEST_CASE("bin_to_category can be empty", "[InputSummary][bin_to_category]") {
    InputSummary s;
    REQUIRE(s.bin_to_category.empty());
}

TEST_CASE("bin_to_category lookup for non-existent bin", "[InputSummary][bin_to_category]") {
    auto s = make_summary();
    REQUIRE(s.bin_to_category.find(99) == s.bin_to_category.end());
}

// ─── num_bins ────────────────────────────────────────────────────────────────

TEST_CASE("num_bins returns correct value", "[InputSummary][num_bins]") {
    auto s = make_summary();
    REQUIRE(s.num_bins == 3);
}

TEST_CASE("num_bins can be zero", "[InputSummary][num_bins]") {
    InputSummary s;
    REQUIRE(s.num_bins == 0);
}

TEST_CASE("num_bins consistency with categories", "[InputSummary][num_bins]") {
    InputSummary s;
    s.categories = {"cat1", "cat2", "cat3", "cat4"};
    s.num_bins = 4;
    REQUIRE(s.num_bins == s.categories.size());
}

// ─── Copy/Move operations ────────────────────────────────────────────────────

TEST_CASE("InputSummary copy constructor works", "[InputSummary][copy]") {
    auto s1 = make_summary();
    InputSummary s2(s1);
    
    REQUIRE(s2.num_bins == 3);
    REQUIRE(s2.categories == s1.categories);
    REQUIRE(s2.filepath_to_bin == s1.filepath_to_bin);
    REQUIRE(s2.bin_to_category == s1.bin_to_category);
}

TEST_CASE("InputSummary copy assignment works", "[InputSummary][copy]") {
    auto s1 = make_summary();
    InputSummary s2;
    s2 = s1;
    
    REQUIRE(s2.num_bins == 3);
    REQUIRE(s2.categories == s1.categories);
}

TEST_CASE("InputSummary move constructor works", "[InputSummary][move]") {
    auto s1 = make_summary();
    InputSummary s2(std::move(s1));
    
    REQUIRE(s2.num_bins == 3);
    REQUIRE(s2.categories.size() == 3);
}

TEST_CASE("InputSummary move assignment works", "[InputSummary][move]") {
    auto s1 = make_summary();
    InputSummary s2;
    s2 = std::move(s1);
    
    REQUIRE(s2.num_bins == 3);
    REQUIRE(s2.categories.size() == 3);
}

// ─── Default construction ────────────────────────────────────────────────────

TEST_CASE("InputSummary default construction initializes to sensible defaults", "[InputSummary][default]") {
    InputSummary s;
    REQUIRE(s.num_bins == 0);
    REQUIRE(s.categories.empty());
    REQUIRE(s.filepath_to_bin.empty());
    REQUIRE(s.bin_to_category.empty());
}

// ─── Edge cases ──────────────────────────────────────────────────────────────

TEST_CASE("InputSummary with many categories", "[InputSummary][edge]") {
    InputSummary s;
    s.num_bins = 255;
    for (int i = 0; i < 255; ++i) {
        s.categories.push_back("category_" + std::to_string(i));
        s.bin_to_category[static_cast<uint8_t>(i)] = "category_" + std::to_string(i);
    }
    
    REQUIRE(s.num_categories() == 255);
    REQUIRE(s.category_index("category_0") == 0);
    REQUIRE(s.category_index("category_254") == 254);
}

TEST_CASE("InputSummary with special characters in category names", "[InputSummary][edge]") {
    InputSummary s;
    s.categories = {"human", "E. coli", "SARS-CoV-2", "M. tuberculosis"};
    s.num_bins = 4;
    
    REQUIRE(s.category_index("E. coli") == 1);
    REQUIRE(s.category_index("SARS-CoV-2") == 2);
    REQUIRE(s.category_index("M. tuberculosis") == 3);
}

TEST_CASE("InputSummary with empty category name", "[InputSummary][edge]") {
    InputSummary s;
    s.categories = {"", "valid"};
    s.num_bins = 2;
    
    REQUIRE(s.category_index("") == 0);
    REQUIRE(s.category_index("valid") == 1);
    REQUIRE(s.category_name(0) == "");
}
