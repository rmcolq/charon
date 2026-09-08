#include <catch2/catch_test_macros.hpp>

#include "input_stats.hpp"

#include <limits>

// Helper: build an InputStats with known bin sizes
static InputStats make_stats() {
    InputStats s;
    s.num_files = 3;
    s.records_per_bin = {{0, 100}, {1, 200}, {2, 50}};
    s.hashes_per_bin  = {{0, 1000}, {1, 5000}, {2, 300}};
    return s;
}

// ─── num_files ───────────────────────────────────────────────────────────────

TEST_CASE("num_files tracks total file count", "[InputStats][num_files]") {
    auto s = make_stats();
    REQUIRE(s.num_files == 3);
}

TEST_CASE("num_files can be zero", "[InputStats][num_files]") {
    InputStats s;
    REQUIRE(s.num_files == 0);
}

TEST_CASE("num_files accumulates correctly", "[InputStats][num_files]") {
    InputStats s;
    s.num_files = 5;
    s.num_files += 3;
    REQUIRE(s.num_files == 8);
}

// ─── records_per_bin ─────────────────────────────────────────────────────────

TEST_CASE("records_per_bin stores correct counts", "[InputStats][records_per_bin]") {
    auto s = make_stats();
    REQUIRE(s.records_per_bin.at(0) == 100);
    REQUIRE(s.records_per_bin.at(1) == 200);
    REQUIRE(s.records_per_bin.at(2) == 50);
}

TEST_CASE("records_per_bin can be empty", "[InputStats][records_per_bin]") {
    InputStats s;
    REQUIRE(s.records_per_bin.empty());
}

TEST_CASE("records_per_bin handles accumulation", "[InputStats][records_per_bin]") {
    InputStats s;
    s.records_per_bin[0] = 100;
    s.records_per_bin[0] += 50;
    REQUIRE(s.records_per_bin[0] == 150);
}

TEST_CASE("records_per_bin with many bins", "[InputStats][records_per_bin]") {
    InputStats s;
    for (uint8_t i = 0; i < 100; ++i) {
        s.records_per_bin[i] = i * 10;
    }
    REQUIRE(s.records_per_bin.size() == 100);
    REQUIRE(s.records_per_bin.at(50) == 500);
}

// ─── hashes_per_bin ──────────────────────────────────────────────────────────

TEST_CASE("hashes_per_bin stores correct counts", "[InputStats][hashes_per_bin]") {
    auto s = make_stats();
    REQUIRE(s.hashes_per_bin.at(0) == 1000);
    REQUIRE(s.hashes_per_bin.at(1) == 5000);
    REQUIRE(s.hashes_per_bin.at(2) == 300);
}

TEST_CASE("hashes_per_bin can be empty", "[InputStats][hashes_per_bin]") {
    InputStats s;
    REQUIRE(s.hashes_per_bin.empty());
}

TEST_CASE("hashes_per_bin handles large values", "[InputStats][hashes_per_bin]") {
    InputStats s;
    s.hashes_per_bin[0] = 1000000000;  // 1 billion
    REQUIRE(s.hashes_per_bin[0] == 1000000000);
}

TEST_CASE("hashes_per_bin handles accumulation", "[InputStats][hashes_per_bin]") {
    InputStats s;
    s.hashes_per_bin[0] = 1000;
    s.hashes_per_bin[0] += 500;
    REQUIRE(s.hashes_per_bin[0] == 1500);
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

TEST_CASE("bins_by_size handles bins with equal hash counts", "[InputStats][bins_by_size]") {
    InputStats s;
    s.hashes_per_bin = {{0, 500}, {1, 500}, {2, 500}};
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.size() == 3);
    // All should have same count
    for (const auto& [bin, count] : sorted) {
        REQUIRE(count == 500);
    }
}

TEST_CASE("bins_by_size preserves bin IDs correctly", "[InputStats][bins_by_size]") {
    InputStats s;
    s.hashes_per_bin = {{5, 100}, {10, 200}, {15, 50}};
    auto sorted = s.bins_by_size();
    
    REQUIRE(sorted[0].first == 15);  // bin 15 has 50 (smallest)
    REQUIRE(sorted[0].second == 50);
    REQUIRE(sorted[1].first == 5);   // bin 5 has 100
    REQUIRE(sorted[1].second == 100);
    REQUIRE(sorted[2].first == 10);  // bin 10 has 200 (largest)
    REQUIRE(sorted[2].second == 200);
}

TEST_CASE("bins_by_size with many bins", "[InputStats][bins_by_size]") {
    InputStats s;
    for (uint8_t i = 0; i < 50; ++i) {
        s.hashes_per_bin[i] = (50 - i) * 100;  // Reverse order
    }
    
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.size() == 50);
    
    // Should be sorted ascending
    for (size_t i = 1; i < sorted.size(); ++i) {
        REQUIRE(sorted[i-1].second <= sorted[i].second);
    }
    
    // Smallest should be last
    REQUIRE(sorted.front().second == 100);  // bin 49
    REQUIRE(sorted.back().second == 5000);  // bin 0
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

TEST_CASE("max_num_hashes with large values", "[InputStats][max_num_hashes]") {
    InputStats s;
    s.hashes_per_bin = {{0, 1000000}, {1, 5000000}, {2, 3000000}};
    REQUIRE(s.max_num_hashes() == 5000000);
}

TEST_CASE("max_num_handles non-sequential bin IDs", "[InputStats][max_num_hashes]") {
    InputStats s;
    s.hashes_per_bin = {{10, 100}, {50, 500}, {100, 50}};
    REQUIRE(s.max_num_hashes() == 500);  // bin 50 has max
}

TEST_CASE("max_num_hashes is const-correct", "[InputStats][max_num_hashes]") {
    const InputStats s = make_stats();
    REQUIRE(s.max_num_hashes() == 5000);
}

// ─── Copy/Move operations ────────────────────────────────────────────────────

TEST_CASE("InputStats copy constructor works", "[InputStats][copy]") {
    auto s1 = make_stats();
    InputStats s2(s1);
    
    REQUIRE(s2.num_files == 3);
    REQUIRE(s2.records_per_bin == s1.records_per_bin);
    REQUIRE(s2.hashes_per_bin == s1.hashes_per_bin);
}

TEST_CASE("InputStats copy assignment works", "[InputStats][copy]") {
    auto s1 = make_stats();
    InputStats s2;
    s2 = s1;
    
    REQUIRE(s2.num_files == 3);
    REQUIRE(s2.hashes_per_bin.at(1) == 5000);
}

TEST_CASE("InputStats move constructor works", "[InputStats][move]") {
    auto s1 = make_stats();
    InputStats s2(std::move(s1));
    
    REQUIRE(s2.num_files == 3);
    REQUIRE(s2.hashes_per_bin.size() == 3);
}

TEST_CASE("InputStats move assignment works", "[InputStats][move]") {
    auto s1 = make_stats();
    InputStats s2;
    s2 = std::move(s1);
    
    REQUIRE(s2.num_files == 3);
    REQUIRE(s2.hashes_per_bin.size() == 3);
}

// ─── Default construction ────────────────────────────────────────────────────

TEST_CASE("InputStats default construction initializes to sensible defaults", "[InputStats][default]") {
    InputStats s;
    REQUIRE(s.num_files == 0);
    REQUIRE(s.records_per_bin.empty());
    REQUIRE(s.hashes_per_bin.empty());
}

// ─── Edge cases ──────────────────────────────────────────────────────────────

TEST_CASE("InputStats with zero records in bin", "[InputStats][edge]") {
    InputStats s;
    s.records_per_bin[0] = 0;
    s.hashes_per_bin[0] = 0;
    
    REQUIRE(s.records_per_bin.at(0) == 0);
    REQUIRE(s.hashes_per_bin.at(0) == 0);
    REQUIRE(s.max_num_hashes() == 0);
}

TEST_CASE("InputStats with very large number of files", "[InputStats][edge]") {
    InputStats s;
    s.num_files = 10000;
    
    for (uint8_t i = 0; i < 100; ++i) {
        s.records_per_bin[i] = 100;
        s.hashes_per_bin[i] = 1000;
    }
    
    REQUIRE(s.num_files == 10000);
    REQUIRE(s.records_per_bin.size() == 100);
}

TEST_CASE("InputStats bins_by_size with single element", "[InputStats][edge]") {
    InputStats s;
    s.hashes_per_bin[42] = 12345;
    
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.size() == 1);
    REQUIRE(sorted[0].first == 42);
    REQUIRE(sorted[0].second == 12345);
}

TEST_CASE("InputStats with maximum uint8 bin ID", "[InputStats][edge]") {
    InputStats s;
    s.hashes_per_bin[255] = 999;
    
    REQUIRE(s.hashes_per_bin.at(255) == 999);
    REQUIRE(s.max_num_hashes() == 999);
    
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.size() == 1);
    REQUIRE(sorted[0].first == 255);
}

// ─── Consistency checks ──────────────────────────────────────────────────────

TEST_CASE("InputStats records_per_bin and hashes_per_bin can have different sizes", "[InputStats][consistency]") {
    InputStats s;
    s.records_per_bin = {{0, 100}, {1, 200}};
    s.hashes_per_bin = {{0, 1000}};
    
    REQUIRE(s.records_per_bin.size() == 2);
    REQUIRE(s.hashes_per_bin.size() == 1);
    
    // bins_by_size should only consider hashes_per_bin
    auto sorted = s.bins_by_size();
    REQUIRE(sorted.size() == 1);
}

TEST_CASE("InputStats num_files independent of bin counts", "[InputStats][consistency]") {
    InputStats s;
    s.num_files = 10;
    s.records_per_bin = {{0, 100}};
    s.hashes_per_bin = {{0, 1000}};
    
    REQUIRE(s.num_files == 10);
    REQUIRE(s.records_per_bin.size() == 1);
}
