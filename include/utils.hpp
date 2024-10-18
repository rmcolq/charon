#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

#include <filesystem>
#include <vector>
#include <string>

// used to transform paths to absolute paths - designed to be used with CLI11 transform
std::filesystem::path make_absolute(std::filesystem::path);

std::vector<std::string> split(const std::string&, const std::string&);

bool ends_with(std::string str, std::string suffix);
bool starts_with(std::string str, std::string prefix);

#endif
