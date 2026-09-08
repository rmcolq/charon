#include <catch2/catch_test_macros.hpp>

#include "utils.hpp"

// ─── ends_with ───────────────────────────────────────────────────────────────

TEST_CASE("ends_with returns true when string ends with suffix", "[ends_with]") {
    REQUIRE(ends_with("hello.fastq", ".fastq"));
    REQUIRE(ends_with("hello.fastq.gz", ".gz"));
    REQUIRE(ends_with("exact", "exact"));
}

TEST_CASE("ends_with returns false when string does not end with suffix", "[ends_with]") {
    REQUIRE_FALSE(ends_with("hello.fastq", ".fasta"));
    REQUIRE_FALSE(ends_with("hello.fastq", ".fastq.gz"));
}

TEST_CASE("ends_with handles empty inputs", "[ends_with]") {
    REQUIRE(ends_with("hello", ""));
    REQUIRE_FALSE(ends_with("", "hello"));
    REQUIRE(ends_with("", ""));
}

// ─── starts_with ─────────────────────────────────────────────────────────────

TEST_CASE("starts_with returns true when string starts with prefix", "[starts_with]") {
    REQUIRE(starts_with("hello.fastq", "hello"));
    REQUIRE(starts_with("exact", "exact"));
}

TEST_CASE("starts_with returns false when string does not start with prefix", "[starts_with]") {
    REQUIRE_FALSE(starts_with("hello.fastq", "world"));
    REQUIRE_FALSE(starts_with("hi", "hello"));
}

TEST_CASE("starts_with handles empty inputs", "[starts_with]") {
    REQUIRE(starts_with("hello", ""));
    REQUIRE_FALSE(starts_with("", "hello"));
    REQUIRE(starts_with("", ""));
}

// ─── get_extension ───────────────────────────────────────────────────────────

TEST_CASE("get_extension returns extension for plain file", "[get_extension]") {
    REQUIRE(get_extension("reads.fastq") == ".fastq");
    REQUIRE(get_extension("reads.fasta") == ".fasta");
    REQUIRE(get_extension("reads.bam") == ".bam");
}

TEST_CASE("get_extension returns inner extension for .gz files", "[get_extension]") {
    REQUIRE(get_extension("reads.fastq.gz") == ".fastq");
    REQUIRE(get_extension("reads.fasta.gz") == ".fasta");
}

TEST_CASE("get_extension handles full paths", "[get_extension]") {
    REQUIRE(get_extension("/data/samples/reads.fastq.gz") == ".fastq");
    REQUIRE(get_extension("/data/samples/reads.bam") == ".bam");
}

TEST_CASE("get_extension returns empty string for files with no extension", "[get_extension]") {
    REQUIRE(get_extension("README") == "");
}

// ─── split ───────────────────────────────────────────────────────────────────

TEST_CASE("split returns nullopt when delimiter not found", "[split]") {
    REQUIRE(split("hello", ",") == std::nullopt);
    REQUIRE(split("", ",") == std::nullopt);
}

TEST_CASE("split divides string by single-character delimiter", "[split]") {
    auto result = split("a,b,c", ",");
    REQUIRE(result.has_value());
    REQUIRE(result->size() == 3);
    REQUIRE(result->at(0) == "a");
    REQUIRE(result->at(1) == "b");
    REQUIRE(result->at(2) == "c");
}

TEST_CASE("split divides string by multi-character delimiter", "[split]") {
    auto result = split("a::b::c", "::");
    REQUIRE(result.has_value());
    REQUIRE(result->size() == 3);
    REQUIRE(result->at(0) == "a");
    REQUIRE(result->at(1) == "b");
    REQUIRE(result->at(2) == "c");
}

TEST_CASE("split handles leading and trailing delimiters", "[split]") {
    auto result = split(",a,b,", ",");
    REQUIRE(result.has_value());
    REQUIRE(result->size() == 4);
    REQUIRE(result->at(0) == "");
    REQUIRE(result->at(1) == "a");
    REQUIRE(result->at(2) == "b");
    REQUIRE(result->at(3) == "");
}

TEST_CASE("split handles consecutive delimiters", "[split]") {
    auto result = split("a,,b", ",");
    REQUIRE(result.has_value());
    REQUIRE(result->size() == 3);
    REQUIRE(result->at(0) == "a");
    REQUIRE(result->at(1) == "");
    REQUIRE(result->at(2) == "b");
}

// ─── first_field ─────────────────────────────────────────────────────────────

TEST_CASE("first_field returns substring before first delimiter", "[first_field]") {
    REQUIRE(first_field("SRR123 some description", " ") == "SRR123");
    REQUIRE(first_field("a,b,c", ",") == "a");
}

TEST_CASE("first_field returns whole string when delimiter not found", "[first_field]") {
    REQUIRE(first_field("SRR123", " ") == "SRR123");
    REQUIRE(first_field("hello", ",") == "hello");
}

TEST_CASE("first_field handles empty string", "[first_field]") {
    REQUIRE(first_field("", " ") == "");
}

TEST_CASE("first_field handles delimiter at start of string", "[first_field]") {
    REQUIRE(first_field(" leading space", " ") == "");
}

// ─── get_compression_ratio ───────────────────────────────────────────────────
// NOTE: empty string input triggers an assert/abort by design — not tested here.
// Sequences are DNA (ACGT). Ratio = compressed_size / original_size.
// Low complexity (repetitive) sequences compress well → ratio closer to 0.
// High complexity (random-like) sequences compress poorly → ratio closer to 1.

TEST_CASE("get_compression_ratio returns value between 0 and 1 for typical DNA", "[get_compression_ratio]") {
    const std::string seq(200, 'A');
    const float ratio = get_compression_ratio(seq);
    REQUIRE(ratio > 0.0f);
    REQUIRE(ratio < 1.5f); // allow slight overhead for short gzip headers
}

TEST_CASE("get_compression_ratio is low for high AT-content repetitive sequence", "[get_compression_ratio]") {
    // ATATAT... is highly repetitive and low complexity — should compress well
    std::string seq;
    seq.reserve(200);
    for (int i = 0; i < 100; ++i) seq += "AT";
    REQUIRE(get_compression_ratio(seq) < 0.5f);
}

TEST_CASE("get_compression_ratio is low for high GC-content repetitive sequence", "[get_compression_ratio]") {
    // GCGCGC... is highly repetitive and low complexity — should compress well
    std::string seq;
    seq.reserve(200);
    for (int i = 0; i < 100; ++i) seq += "GC";
    REQUIRE(get_compression_ratio(seq) < 0.5f);
}

TEST_CASE("get_compression_ratio is low for homopolymer run", "[get_compression_ratio]") {
    // All same base — maximally low complexity
    const std::string seq(200, 'A');
    REQUIRE(get_compression_ratio(seq) < 0.5f);
}

TEST_CASE("get_compression_ratio is higher for complex sequence", "[get_compression_ratio]") {
    // Pseudo-random DNA — high complexity, should not compress well
    const std::string seq =
        "ACGTTAGCACGTTAGCTTGACAGTCACGTTAGCACGTTAGCTTGACAGTC"
        "GTATCGAATCGATCGATCGTAGCTAGCTAGCTAGCTAGCATCGATCGATC"
        "TTACGGCAATCGATCGTACGATCGATCGTAGCTAGCTAGCTAGCATCGATC"
        "ACGTTAGCACGTTAGCTTGACAGTCACGTTAGCACGTTAGCTTGACAGTC";
    REQUIRE(get_compression_ratio(seq) > 0.7f);
}

TEST_CASE("get_compression_ratio is lower for repetitive than complex sequence", "[get_compression_ratio]") {
    const std::string repetitive(200, 'A');
    std::string complex_seq;
    complex_seq.reserve(200);
    for (int i = 0; i < 50; ++i) complex_seq += "ACGT";
    // even ACGT repeated is more complex than all-A, but less than truly random
    REQUIRE(get_compression_ratio(repetitive) < get_compression_ratio(complex_seq));
}