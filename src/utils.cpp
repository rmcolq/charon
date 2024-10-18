#include "utils.hpp"

#include <string_view>

std::filesystem::path make_absolute(std::filesystem::path path) { return std::filesystem::absolute(path); }

std::vector<std::string> split(const std::string& s, const std::string& delimiter){
    std::vector<std::string> substrings;

    int start, end = -1*delimiter.size();
    do {
        start = end + delimiter.size();
        end = s.find(delimiter, start);
        substrings.push_back(s.substr(start, end - start));
    } while (end != -1);

    return substrings;
}

/*static bool ends_with(std::string_view str, std::string_view suffix)
{
    return str.size() >= suffix.size() && str.compare(str.size()-suffix.size(), suffix.size(), suffix) == 0;
}

static bool starts_with(std::string_view str, std::string_view prefix)
{
    return str.size() >= prefix.size() && str.compare(0, prefix.size(), prefix) == 0;
}*/