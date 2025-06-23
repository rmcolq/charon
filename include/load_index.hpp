#ifndef CHARON_LOAD_INDEX_H
#define CHARON_LOAD_INDEX_H

#pragma once

#include <filesystem>
#include <index.hpp>

void load_index(Index &index, std::filesystem::path const &path);

#endif // CHARON_LOAD_INDEX_MAIN_H