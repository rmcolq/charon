#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#include <filesystem>
#include <vector>
#include <string>
#include <string_view>

// used to transform paths to absolute paths - designed to be used with CLI11 transform
std::filesystem::path make_absolute(std::filesystem::path);

std::vector<std::string> split(const std::string&, const std::string&);

//static bool ends_with(std::string_view str, std::string_view suffix);
//static bool starts_with(std::string_view str, std::string_view prefix);

#endif
