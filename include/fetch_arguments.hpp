#ifndef CHARON_FETCH_ARGUMENTS_H
#define CHARON_FETCH_ARGUMENTS_H

#pragma once

#include <string>

/// Collection of all options of the fetch subcommand.
struct FetchArguments {
    std::string tech{"nanopore"};  // "nanopore" or "illumina"; empty means every technology (--list only)
    std::string size{"big"};     // "small" or "big"
    bool host_only{false};
    std::string outdir{};        // empty: $CHARON_DB, else the current directory
    std::string record{};        // pin a Zenodo record id instead of following the latest version
    std::string api_url{};       // Zenodo API base; empty: $CHARON_ZENODO_API, else https://zenodo.org/api
    bool list{false};
    bool keep_gz{false};
    bool force{false};
};

#endif // CHARON_FETCH_ARGUMENTS_H
