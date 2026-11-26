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

// ─── shorten_filepaths ───────────────────────────────────────────────────────

static std::vector<std::string> shortened(const std::vector<std::string>& paths) {
    InputSummary s;
    for (size_t i = 0; i < paths.size(); ++i)
        s.filepath_to_bin.emplace_back(paths[i], static_cast<uint8_t>(i));
    s.shorten_filepaths();
    std::vector<std::string> result;
    for (const auto& [path, bin] : s.filepath_to_bin)
        result.push_back(path);
    return result;
}

TEST_CASE("shorten_filepaths keeps only file names when they are unique", "[InputSummary][shorten_filepaths]") {
    REQUIRE(shortened({"/home/me/db/human_chr1.fa.gz", "/home/me/db/viruses_01.fa.gz", "relative/dir/fungi.fa"}) ==
            std::vector<std::string>{"human_chr1.fa.gz", "viruses_01.fa.gz", "fungi.fa"});
}

TEST_CASE("shorten_filepaths adds directories to tell apart files with the same name",
          "[InputSummary][shorten_filepaths]") {
    REQUIRE(shortened({"/data/a/ref.fa", "/data/b/ref.fa", "/data/c/other.fa"}) ==
            std::vector<std::string>{"a/ref.fa", "b/ref.fa", "other.fa"});
    REQUIRE(shortened({"/x/same/ref.fa", "/y/same/ref.fa", "/z/ref.fa"}) ==
            std::vector<std::string>{"x/same/ref.fa", "y/same/ref.fa", "z/ref.fa"});
}

TEST_CASE("shorten_filepaths handles one path ending another and duplicate paths",
          "[InputSummary][shorten_filepaths]") {
    REQUIRE(shortened({"/x/f.fa", "/y/x/f.fa"}) == std::vector<std::string>{"x/f.fa", "y/x/f.fa"});
    REQUIRE(shortened({"/d/f.fa", "/d/f.fa"}) == std::vector<std::string>{"d/f.fa", "d/f.fa"});
}

TEST_CASE("shorten_filepaths keeps bins and leaves no root", "[InputSummary][shorten_filepaths]") {
    InputSummary s;
    s.filepath_to_bin = {{"/only.fa", 7}};
    s.shorten_filepaths();
    REQUIRE(s.filepath_to_bin.front() == std::pair<std::string, uint8_t>{"only.fa", 7});
}
