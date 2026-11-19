#include <catch2/catch_test_macros.hpp>

#include <sstream>
#include <string>
#include <vector>

#include "read_entry.hpp"

namespace {
std::vector<std::string> split_tabs(const std::string &line) {
    std::vector<std::string> fields;
    std::stringstream ss(line);
    std::string field;
    while (std::getline(ss, field, '\t'))
        fields.push_back(field);
    return fields;
}
}

TEST_CASE("a filtered read has a row in the same columns as a classified read", "[filtered]") {
    const auto line = filtered_assignment_line("read1", 120, 17.5f, FilterReason::too_short);
    REQUIRE(line.back() == '\n');
    // status, read_id, call, length, num_hashes, mean_quality, confidence, compression, details
    const auto fields = split_tabs(line.substr(0, line.size() - 1));
    REQUIRE(fields.size() == 9);
    REQUIRE(fields[0] == "F");
    REQUIRE(fields[1] == "read1");
    REQUIRE(fields[2] == "");          // no call
    REQUIRE(fields[3] == "120");
    REQUIRE(fields[4] == "0");         // no hashes were computed
    REQUIRE(fields[5] == "17.5");
    REQUIRE(fields[6] == "0");         // no confidence
    REQUIRE(fields[7] == "NA");        // compression was not computed
    REQUIRE(fields[8] == "filter:too_short");
}

TEST_CASE("the reason distinguishes length from quality", "[filtered]") {
    const auto short_line = filtered_assignment_line("a", 10, 30.0f, FilterReason::too_short);
    const auto quality_line = filtered_assignment_line("a", 10, 3.0f, FilterReason::low_quality);
    REQUIRE(short_line.find("filter:too_short") != std::string::npos);
    REQUIRE(quality_line.find("filter:low_quality") != std::string::npos);
    REQUIRE(quality_line.find("filter:too_short") == std::string::npos);
}

TEST_CASE("a filtered row is not mistaken for classified or unclassified", "[filtered]") {
    const auto line = filtered_assignment_line("a", 10, 30.0f, FilterReason::too_short);
    REQUIRE(line.rfind("F\t", 0) == 0);
    REQUIRE(line.rfind("U\t", 0) != 0);
    REQUIRE(line.rfind("C\t", 0) != 0);
}
