#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <string>
#include <vector>

#include "dehost_arguments.hpp"

// These tests check how defaults are chosen, not the values themselves, so that tuning a default needs no test change.

TEST_CASE("paired reads take the paired default for every option not given", "[defaults]") {
    DehostArguments opt;
    opt.apply_paired_defaults();
    REQUIRE(opt.min_length == DehostArguments::paired_min_length);
    REQUIRE(opt.min_compression == DehostArguments::paired_min_compression);
    REQUIRE(opt.host_unique_prop_lo_threshold == DehostArguments::paired_host_unique_prop_lo_threshold);
    REQUIRE(opt.min_proportion_difference == DehostArguments::paired_min_proportion_difference);
}

TEST_CASE("paired defaults leave options without a paired default unchanged", "[defaults]") {
    const DehostArguments single;
    DehostArguments opt;
    opt.apply_paired_defaults();
    REQUIRE(opt.min_quality == single.min_quality);
    REQUIRE(opt.confidence_threshold == single.confidence_threshold);
    REQUIRE(opt.min_prob_difference == single.min_prob_difference);
}

TEST_CASE("an option given explicitly is not replaced by a paired default", "[defaults]") {
    // explicit values that differ from both the single-end and the paired defaults
    const DehostArguments single;
    const uint32_t min_length = single.min_length + DehostArguments::paired_min_length + 1;
    const float min_compression = single.min_compression + DehostArguments::paired_min_compression + 0.01f;
    const float lo = single.host_unique_prop_lo_threshold + DehostArguments::paired_host_unique_prop_lo_threshold + 0.01f;
    const float mpd = single.min_proportion_difference + DehostArguments::paired_min_proportion_difference + 0.01f;

    DehostArguments opt;
    opt.min_length = min_length;
    opt.min_length_set = true;
    opt.min_compression = min_compression;
    opt.min_compression_set = true;
    opt.host_unique_prop_lo_threshold = lo;
    opt.host_unique_prop_lo_threshold_set = true;
    opt.min_proportion_difference = mpd;
    opt.min_proportion_difference_set = true;
    opt.apply_paired_defaults();
    REQUIRE(opt.min_length == min_length);
    REQUIRE(opt.min_compression == min_compression);
    REQUIRE(opt.host_unique_prop_lo_threshold == lo);
    REQUIRE(opt.min_proportion_difference == mpd);
}

TEST_CASE("only the options given explicitly keep their value", "[defaults]") {
    const DehostArguments single;
    DehostArguments opt;
    opt.min_length = single.min_length + DehostArguments::paired_min_length + 1;
    opt.min_length_set = true;
    opt.apply_paired_defaults();
    REQUIRE(opt.min_length == single.min_length + DehostArguments::paired_min_length + 1);
    REQUIRE(opt.min_proportion_difference == DehostArguments::paired_min_proportion_difference);
}
