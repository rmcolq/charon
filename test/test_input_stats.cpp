#include <catch2/catch_test_macros.hpp>

#include "input_stats.hpp"

// Helper: build an InputStats with known bin sizes
static InputStats make_stats() {
    InputStats s;
    s.num_files = 3;
    s.records_per_bin = {{0, 100}, {1, 200}, {2, 50}};
    s.hashes_per_bin  = {{0, 1000}, {1, 5000}, {2, 300}};
    return s;
}

// ─── bins_by_size ────────────────────────────────────────────────────────────

TEST_CASE("bins_by_size returns bins sorted ascending by hash count", "[InputStats][bins_by_size]") {
    auto s = make_stats();
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.size() == 3);
    REQUIRE(sorted[0].second <= sorted[1].second);
    REQUIRE(sorted[1].second <= sorted[2].second);
}

TEST_CASE("bins_by_size smallest bin is correct", "[InputStats][bins_by_size]") {
    auto s = make_stats(); // bin 2 has 300 hashes — smallest
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.front().first == 2);
    REQUIRE(sorted.front().second == 300);
}

TEST_CASE("bins_by_size largest bin is correct", "[InputStats][bins_by_size]") {
    auto s = make_stats(); // bin 1 has 5000 hashes — largest
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.back().first == 1);
    REQUIRE(sorted.back().second == 5000);
}

TEST_CASE("bins_by_size returns empty vector for empty stats", "[InputStats][bins_by_size]") {
    InputStats s;
    REQUIRE(s.bins_by_size().empty());
}

TEST_CASE("bins_by_size handles single bin", "[InputStats][bins_by_size]") {
    InputStats s;
    s.hashes_per_bin = {{0, 42}};
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.size() == 1);
    REQUIRE(sorted[0].first == 0);
    REQUIRE(sorted[0].second == 42);
}

// ─── max_num_hashes ──────────────────────────────────────────────────────────

TEST_CASE("max_num_hashes returns largest hash count", "[InputStats][max_num_hashes]") {
    auto s = make_stats(); // bin 1 has 5000 hashes
    REQUIRE(s.max_num_hashes() == 5000);
}

TEST_CASE("max_num_hashes returns 0 for empty stats", "[InputStats][max_num_hashes]") {
    // Fixed: original implementation called .back() on empty vector — undefined behaviour
    InputStats s;
    REQUIRE(s.max_num_hashes() == 0);
}

TEST_CASE("max_num_hashes handles single bin", "[InputStats][max_num_hashes]") {
    InputStats s;
    s.hashes_per_bin = {{0, 99}};
    REQUIRE(s.max_num_hashes() == 99);
}

TEST_CASE("max_num_hashes handles bins with equal hash counts", "[InputStats][max_num_hashes]") {
    InputStats s;
    s.hashes_per_bin = {{0, 500}, {1, 500}, {2, 500}};
    REQUIRE(s.max_num_hashes() == 500);
}
