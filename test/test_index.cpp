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

// ─── Layout (balanced technical bins) ────────────────────────────────────────

#include <filesystem>
#include <fstream>
#include <cereal/archives/binary.hpp>

#include "layout.hpp"
#include "load_index.hpp"
#include "store_index.hpp"

TEST_CASE("Identity layout maps each group to its own technical bin", "[Layout]") {
    const auto layout = Layout::identity(5);
    REQUIRE(layout.num_groups() == 5);
    REQUIRE(layout.num_technical_bins == 5);
    REQUIRE(layout.is_identity());
    for (size_t group = 0; group < 5; ++group)
        REQUIRE(layout.technical_bin(group, 0xABCDEFULL * group) == group);
}

TEST_CASE("Split groups map hashes into their own technical bins", "[Layout]") {
    Layout layout;
    layout.first_bin = {0, 1, 4};
    layout.num_parts = {1, 3, 1};
    layout.num_technical_bins = 5;
    REQUIRE_FALSE(layout.is_identity());
    for (uint64_t hash = 0; hash < 1000; ++hash) {
        REQUIRE(layout.technical_bin(0, hash) == 0);
        const auto technical_bin = layout.technical_bin(1, hash);
        REQUIRE(technical_bin >= 1);
        REQUIRE(technical_bin < 4);
        REQUIRE(technical_bin == layout.technical_bin(1, hash));  // deterministic
        REQUIRE(layout.technical_bin(2, hash) == 4);
    }
}

TEST_CASE("Layout partitions k-mer-like hashes evenly", "[Layout]") {
    // seqan3 minimiser hashes are k-mer codes XOR a seed, so consecutive codes are a hard case
    constexpr size_t parts = 7;
    constexpr uint64_t seed = 0x8F3F73B5CF1C9ADEULL;
    std::vector<size_t> counts(parts, 0);
    constexpr size_t n = 70000;
    for (uint64_t code = 0; code < n; ++code)
        counts[Layout::part(code ^ seed, parts)] += 1;
    for (const auto count: counts) {
        REQUIRE(count > n / parts * 0.9);
        REQUIRE(count < n / parts * 1.1);
    }
}

static std::filesystem::path temp_index_path(const std::string &name) {
    return std::filesystem::temp_directory_path() / ("charon_test_" + name + ".idx");
}

TEST_CASE("Index with a balanced layout round-trips through store/load", "[Index][Layout][io]") {
    auto args = make_arguments();
    auto summary = make_summary({"human", "bacteria"});
    Layout layout;
    layout.first_bin = {0, 1};
    layout.num_parts = {1, 2};
    layout.num_technical_bins = 3;
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf{
            seqan3::bin_count{3u}, seqan3::bin_size{1024u}, seqan3::hash_function_count{2u}};
    ibf.emplace(42u, seqan3::bin_index{layout.technical_bin(1, 42u)});

    const auto path = temp_index_path("v4");
    store_index(path, Index{args, summary, make_stats(), ibf, layout});
    Index loaded;
    load_index(loaded, path);
    std::filesystem::remove(path);

    REQUIRE_FALSE(loaded.is_compressed());
    REQUIRE(loaded.technical_bins() == 3);
    REQUIRE(loaded.num_bins() == 2);
    REQUIRE(loaded.layout().first_bin == layout.first_bin);
    REQUIRE(loaded.layout().num_parts == layout.num_parts);
    bool found = false;
    loaded.with_agent([&](auto &agent) {
        found = agent.bulk_contains(42u)[layout.technical_bin(1, 42u)];
    });
    REQUIRE(found);
}

TEST_CASE("Index written in the pre-v4 format still loads, compressed", "[Index][io]") {
    auto args = make_arguments();
    auto summary = make_summary({"human", "bacteria", "virus"});
    auto stats = make_stats();
    seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed> ibf{
            seqan3::bin_count{3u}, seqan3::bin_size{1024u}, seqan3::hash_function_count{2u}};
    ibf.emplace(7u, seqan3::bin_index{2u});
    const seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed> compressed{ibf};

    // the old format: the cereal archive straight away, no magic, no layout, compressed IBF
    const auto path = temp_index_path("v3");
    {
        std::ofstream os{path, std::ios::binary};
        cereal::BinaryOutputArchive archive{os};
        archive(args.window_size, args.kmer_size, args.max_fpr, summary, stats, compressed);
    }
    Index loaded;
    load_index(loaded, path);
    std::filesystem::remove(path);

    REQUIRE(loaded.is_compressed());
    REQUIRE(loaded.window_size() == args.window_size);
    REQUIRE(loaded.num_bins() == 3);
    REQUIRE(loaded.layout().is_identity());
    bool found = false;
    loaded.with_agent([&](auto &agent) { found = agent.bulk_contains(7u)[2]; });
    REQUIRE(found);
}
