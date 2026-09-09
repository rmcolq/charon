#include <catch2/catch_test_macros.hpp>
#include <cinttypes>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>

#include "index.hpp"

// ─── Test fixture helpers ─────────────────────────────────────────────────────

static InputSummary make_summary(const std::vector<std::string>& categories) {
    InputSummary s;
    s.num_bins = static_cast<uint8_t>(categories.size());
    s.categories = categories;
    for (uint8_t i = 0; i < categories.size(); ++i) {
        s.bin_to_category[i] = categories[i];
    }
    return s;
}

static InputStats make_stats(uint32_t num_files = 3) {
    InputStats s;
    s.num_files = num_files;
    return s;
}

// Build a minimal compressed IBF with n bins
static seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>
make_ibf(uint8_t num_bins) {
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf{
        seqan3::bin_count{static_cast<size_t>(num_bins)},
        seqan3::bin_size{128u},
        seqan3::hash_function_count{2u}
    };
    return seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>{ibf};
}

// Build a minimal IndexArguments-like struct for Index construction.
// We only need the fields the Index constructor uses.
// Uses IndexArguments struct defaults: window_size=41, kmer_size=19, max_fpr=0.01
static IndexArguments make_arguments() {
    return IndexArguments{};
}

// Build a complete minimal Index for reuse across tests
static Index make_index(const std::vector<std::string>& categories = {"human", "bacteria", "virus"}) {
    auto args = make_arguments();
    auto summary = make_summary(categories);
    auto stats = make_stats();
    auto ibf = make_ibf(static_cast<uint8_t>(categories.size()));
    return Index{args, summary, stats, ibf};
}

// ─── basic accessors ─────────────────────────────────────────────────────────

TEST_CASE("Index stores and returns window_size", "[Index][accessors]") {
    auto idx = make_index();
    REQUIRE(idx.window_size() == 41);
}

TEST_CASE("Index stores and returns kmer_size", "[Index][accessors]") {
    auto idx = make_index();
    REQUIRE(idx.kmer_size() == 19);
}

TEST_CASE("Index stores and returns max_fpr", "[Index][accessors]") {
    auto idx = make_index();
    REQUIRE(idx.max_fpr() == 0.01);
}

TEST_CASE("Index returns correct num_bins", "[Index][accessors]") {
    auto idx = make_index();
    REQUIRE(idx.num_bins() == 3);
}

TEST_CASE("Index returns correct num_categories", "[Index][accessors]") {
    auto idx = make_index();
    REQUIRE(idx.num_categories() == 3);
}

TEST_CASE("Index returns correct categories", "[Index][accessors]") {
    auto idx = make_index();
    auto cats = idx.categories();
    REQUIRE(cats.size() == 3);
    REQUIRE(cats[0] == "human");
    REQUIRE(cats[1] == "bacteria");
    REQUIRE(cats[2] == "virus");
}

// ─── get_host_index ──────────────────────────────────────────────────────────

TEST_CASE("get_host_index finds 'human' category", "[Index][get_host_index]") {
    auto idx = make_index({"human", "bacteria", "virus"});
    REQUIRE(idx.get_host_index() == 0);
}

TEST_CASE("get_host_index finds 'host' category", "[Index][get_host_index]") {
    auto idx = make_index({"bacteria", "host", "virus"});
    REQUIRE(idx.get_host_index() == 1);
}

TEST_CASE("get_host_index returns minimum when both human and host present",
          "[Index][get_host_index]") {
    auto idx = make_index({"host", "human", "virus"}); // host=0, human=1
    REQUIRE(idx.get_host_index() == 0);
}

// ─── get_category_index ──────────────────────────────────────────────────────

TEST_CASE("get_category_index returns correct index for known category",
          "[Index][get_category_index]") {
    auto idx = make_index({"human", "bacteria", "virus"});
    REQUIRE(idx.get_category_index("human") == 0);
    REQUIRE(idx.get_category_index("bacteria") == 1);
    REQUIRE(idx.get_category_index("virus") == 2);
}

// ─── bin_to_category ─────────────────────────────────────────────────────────

TEST_CASE("bin_to_category map is correctly stored and returned", "[Index][bin_to_category]") {
    auto idx = make_index({"human", "bacteria", "virus"});
    auto map = idx.bin_to_category();
    REQUIRE(map.at(0) == "human");
    REQUIRE(map.at(1) == "bacteria");
    REQUIRE(map.at(2) == "virus");
}

// ─── default construction ────────────────────────────────────────────────────

TEST_CASE("Default-constructed Index has zero window and kmer size", "[Index][default]") {
    Index idx;
    REQUIRE(idx.window_size() == 0);
    REQUIRE(idx.kmer_size() == 0);
}
