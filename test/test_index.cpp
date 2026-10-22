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

// ─── ConcurrentIbfInserter ────────────────────────────────────────────────────

#include <random>

#include "concurrent_ibf.hpp"
#include "hyperloglog.hpp"

TEST_CASE("ConcurrentIbfInserter sets the same bits as emplace", "[ConcurrentIbfInserter]") {
    for (const size_t num_bins: {1, 63, 64, 70, 300}) {
        for (const size_t num_hash: {1, 2, 3, 5}) {
            const size_t num_bits = 10007;
            Index::ibf_uncompressed expected{seqan3::bin_count{num_bins}, seqan3::bin_size{num_bits},
                                             seqan3::hash_function_count{num_hash}};
            Index::ibf_uncompressed actual{seqan3::bin_count{num_bins}, seqan3::bin_size{num_bits},
                                           seqan3::hash_function_count{num_hash}};
            std::mt19937_64 rng{num_bins * 10 + num_hash};
            std::vector<std::pair<uint64_t, size_t>> values;
            for (size_t i = 0; i < 20000; ++i)
                values.emplace_back(rng(), rng() % num_bins);
            for (const auto &[value, bin]: values)
                expected.emplace(value, seqan3::bin_index{bin});
            const ConcurrentIbfInserter inserter{actual};
#pragma omp parallel for num_threads(4)
            for (size_t i = 0; i < values.size(); ++i)
                inserter.emplace(values[i].first, values[i].second);
            REQUIRE(actual == expected);
        }
    }
}

// ─── HyperLogLog ──────────────────────────────────────────────────────────────

TEST_CASE("HyperLogLog estimates distinct values within 1 %", "[HyperLogLog]") {
    for (const uint64_t n: {0ULL, 1ULL, 100ULL, 10000ULL, 300000ULL, 3000000ULL}) {
        HyperLogLog sketch;
        // Small consecutive integers, like k-mer codes, and each added twice
        for (uint64_t i = 0; i < n; ++i) {
            sketch.add(i);
            sketch.add(i);
        }
        const auto estimate = static_cast<double>(sketch.estimate());
        if (n == 0)
            REQUIRE(estimate == 0);
        else
            REQUIRE(std::abs(estimate - static_cast<double>(n)) <= std::max(1.0, 0.01 * static_cast<double>(n)));
    }
}

// ─── for_each_fasta_chunk ─────────────────────────────────────────────────────

#include <set>

#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

#include "fasta_chunks.hpp"

TEST_CASE("for_each_fasta_chunk yields the same minimisers as whole records", "[for_each_fasta_chunk]") {
    std::mt19937_64 rng{42};
    std::string fasta;
    for (size_t record = 0; record < 6; ++record) {
        fasta += ">record " + std::to_string(record) + " with a description\n";
        const size_t length = record == 3 ? 10 : 500 + 731 * record;  // one record shorter than the window
        for (size_t i = 0; i < length; ++i) {
            fasta += "ACGTacgtNRY"[rng() % (i % 97 == 0 ? 11 : 8)];
            if (i % 60 == 59)
                fasta += i % 120 == 119 ? "\r\n" : " 12\n";
        }
        fasta += "\n";
    }
    const auto dir = std::filesystem::temp_directory_path() / "charon_test_fasta_chunks";
    std::filesystem::create_directories(dir);
    const auto plain = dir / "test.fa";
    const auto compressed = dir / "test.fa.gz";
    std::ofstream{plain} << fasta;
    gzFile out = gzopen(compressed.c_str(), "wb");
    gzwrite(out, fasta.data(), static_cast<unsigned>(fasta.size()));
    gzclose(out);

    const uint8_t window = 23;
    const auto hash_adaptor = seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{15}},
                                                            seqan3::window_size{window});
    std::set<uint64_t> expected;
    uint64_t expected_records = 0;
    seqan3::sequence_file_input fin{plain};
    for (const auto &record: fin) {
        for (const auto value: record.sequence() | hash_adaptor)
            expected.insert(value);
        ++expected_records;
    }

    REQUIRE(is_chunkable_fasta(compressed));
    REQUIRE(not is_chunkable_fasta(dir / "test.fq.gz"));
    for (const auto &file: {plain, compressed}) {
        for (const size_t chunk_size: {1, 7, 100, 1000000}) {
            std::set<uint64_t> actual;
            const auto records = for_each_fasta_chunk(file, chunk_size, window - 1,
                                                      [&](const std::vector<seqan3::dna5> &chunk) {
                                                          for (const auto value: chunk | hash_adaptor)
                                                              actual.insert(value);
                                                      });
            REQUIRE(records == expected_records);
            REQUIRE(actual == expected);
        }
    }
    std::filesystem::remove_all(dir);
}
