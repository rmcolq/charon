#include <catch2/catch_test_macros.hpp>

#include "memory_manager.hpp"

namespace {
constexpr size_t MB = 1024 * 1024;
constexpr size_t GB = 1024 * MB;
}

TEST_CASE("explicit memory limit is used as given", "[memory_manager]") {
    MemoryManager mgr(8 * 1024);
    REQUIRE(mgr.max_allowed() == 8 * GB);
}

TEST_CASE("index and overhead come out of the working memory", "[memory_manager]") {
    MemoryManager mgr(10 * 1024);
    const auto without_index = mgr.working_memory(4);
    mgr.set_index_bytes(4 * GB);
    const auto with_index = mgr.working_memory(4);
    REQUIRE(with_index < without_index);
    // the 4GB index must be taken off before the safety margin is applied
    REQUIRE(without_index - with_index >= 4 * GB * 0.75);
}

TEST_CASE("more threads leave less working memory", "[memory_manager]") {
    MemoryManager mgr(10 * 1024);
    mgr.set_index_bytes(4 * GB);
    REQUIRE(mgr.working_memory(8) < mgr.working_memory(1));
}

TEST_CASE("an index larger than the limit is infeasible", "[memory_manager]") {
    MemoryManager mgr(4 * 1024);
    mgr.set_index_bytes(6 * GB);
    REQUIRE(mgr.working_memory(1) == 0);
    REQUIRE_FALSE(mgr.is_feasible(1));
}

TEST_CASE("an index that leaves almost nothing is infeasible", "[memory_manager]") {
    MemoryManager mgr(4 * 1024);
    mgr.set_index_bytes(3 * GB + 800 * MB);
    REQUIRE_FALSE(mgr.is_feasible(1));
}

TEST_CASE("a comfortable limit is feasible", "[memory_manager]") {
    MemoryManager mgr(16 * 1024);
    mgr.set_index_bytes(6 * GB);
    REQUIRE(mgr.is_feasible(6));
}

TEST_CASE("chunk and cache budgets fit within the working memory", "[memory_manager]") {
    MemoryManager mgr(12 * 1024);
    mgr.set_index_bytes(6 * GB);
    const auto working = mgr.working_memory(6);
    // buffered reads cost more than 1 byte per base, so the chunk budget in bytes is below working memory
    REQUIRE(mgr.chunk_base_budget(6) * 2 + mgr.cache_byte_budget(6) <= working);
}

TEST_CASE("budgets shrink as the index grows", "[memory_manager]") {
    MemoryManager small_index(12 * 1024);
    small_index.set_index_bytes(2 * GB);
    MemoryManager big_index(12 * 1024);
    big_index.set_index_bytes(8 * GB);
    REQUIRE(big_index.cache_byte_budget(4) < small_index.cache_byte_budget(4));
}

TEST_CASE("chunk budget does not depend on read length", "[memory_manager]") {
    // chunks are sized in bases, so the budget is a single number
    MemoryManager mgr(12 * 1024);
    mgr.set_index_bytes(4 * GB);
    REQUIRE(mgr.chunk_base_budget(4) >= 1000000);
    REQUIRE(mgr.chunk_base_budget(4) <= 256ULL * 1000 * 1000 * 2);
}

TEST_CASE("cached read cost grows with read length and mates", "[memory_manager]") {
    const auto short_read = MemoryManager::cached_read_bytes(1000, 20, 2, 100);
    const auto long_read = MemoryManager::cached_read_bytes(100000, 20, 2, 10000);
    REQUIRE(long_read > short_read * 50);
    const auto paired = MemoryManager::cached_read_bytes(1000, 20, 2, 100, true);
    REQUIRE(paired > short_read);
}

TEST_CASE("resident memory can be measured on supported platforms", "[memory_manager]") {
#if defined(__linux__) || defined(__APPLE__)
    REQUIRE(MemoryManager::current_rss() > 0);
#else
    SUCCEED();
#endif
}

TEST_CASE("peak memory is tracked", "[memory_manager]") {
    MemoryManager mgr(64 * 1024);
    REQUIRE(mgr.within_limit());
    REQUIRE(mgr.peak_rss() == MemoryManager::current_rss());
}

TEST_CASE("within_limit is false once usage passes the fraction of the limit", "[memory_manager]") {
    MemoryManager mgr(1);  // 1 MB limit; the test process is far above it
    if (MemoryManager::current_rss() > 0) {
        REQUIRE_FALSE(mgr.within_limit(0.9));
    }
}

TEST_CASE("format_size uses sensible units", "[memory_manager]") {
    REQUIRE(MemoryManager::format_size(512) == "512.00 B");
    REQUIRE(MemoryManager::format_size(2 * GB) == "2.00 GB");
}
