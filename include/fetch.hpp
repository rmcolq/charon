#ifndef CHARON_FETCH_H
#define CHARON_FETCH_H

#pragma once

// `charon fetch`: download a published Charon index from Zenodo, check its md5 and decompress it.

#include <cstdint>
#include <filesystem>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "CLI11.hpp"
#include "fetch_arguments.hpp"

class Md5;

namespace fetch {

/// A file in a Zenodo record.
struct RemoteFile {
    std::string key;   // file name
    uint64_t size{0};  // bytes
    std::string md5;   // lower-case hex
    std::string url;   // download url
};

/// The parts of a Zenodo record that fetch uses.
struct RecordInfo {
    std::string id;
    std::string doi;
    std::string version;
    std::string publication_date;
    std::vector<RemoteFile> files;
};

/// What a fetched index was downloaded from, kept next to it in <index>.source.tsv.
struct Source {
    std::string file;
    std::string record;
    std::string doi;
    std::string version;
    std::string md5;  // md5 of the .gz on Zenodo
};

/// An HTTP (or file://) request that failed; status is the HTTP status, 404 for a missing file:// path, or 0.
struct RequestError : std::runtime_error {
    long status;

    RequestError(const std::string &what, long status_) : std::runtime_error(what), status(status_) {}
};

inline constexpr std::string_view default_api_url{"https://zenodo.org/api"};

/// Technologies with a published index, in --list order.
const std::vector<std::string> &technologies();

/// A Zenodo record id for this technology's index. Fetch follows it to the latest version unless --record is given.
std::string record_for(const std::string &tech);

/// --api-url, else $CHARON_ZENODO_API, else https://zenodo.org/api (without a trailing slash).
std::string api_base(const FetchArguments &opt);

/// The API url of a record, or of the latest version of the record's concept.
std::string record_url(const std::string &base, const std::string &record_id, bool latest);

/// -o, else $CHARON_DB, else the current directory.
std::filesystem::path output_dir(const FetchArguments &opt);

/// Parse a Zenodo API record (JSON). Throws std::runtime_error if it is not one.
RecordInfo parse_record(std::string_view json);

/// charon_<size>[_host_only].<tech>.idx.gz
std::string index_file_name(const std::string &tech, const std::string &size, bool host_only);

/// The record's file with this name. Throws std::runtime_error listing the indexes the record does have.
const RemoteFile &select_file(const RecordInfo &record, const std::string &name);

/// The command-line options that select this index file, e.g. "--tech nanopore --size big --host-only".
/// Empty if the name does not follow the charon_<size>[_host_only].<tech>.idx.gz pattern.
std::string selection_flags(const std::string &name);

/// e.g. "6.26 GB"
std::string human_size(uint64_t bytes);

/// Lower-case hex md5 of a file (bytes [0, length) if length is given).
std::string file_md5(const std::filesystem::path &path, std::optional<uint64_t> length = std::nullopt);

/// Feed the first `length` bytes of a file to an Md5.
void hash_file_prefix(Md5 &md5, const std::filesystem::path &path, uint64_t length);

/// Throws std::runtime_error if these leading bytes of a decompressed file are not an index this Charon can load.
/// Returns a short description of the format, e.g. "v4".
std::string check_index_start(std::string_view first_bytes);

/// Decompress a gzip file to `out`, checking at the start that it is an index this Charon can load.
/// Writes to `out` directly; the caller renames it into place.
void gunzip_index(const std::filesystem::path &gz, const std::filesystem::path &out);

std::filesystem::path source_path(const std::filesystem::path &index);

void write_source(const std::filesystem::path &index, const Source &source);

std::optional<Source> read_source(const std::filesystem::path &index);

// ---- network (libcurl; file:// urls work too, which the tests use) ----

/// GET a url into a string. Throws RequestError.
std::string http_get(const std::string &url);

/// Download a gzipped index's first bytes, decompress them and run check_index_start on them.
std::string check_remote_index(const std::string &url);

/// Download `url` into `part`, resuming from what `part` already holds, retrying failed transfers,
/// and return the md5 of the whole file.
std::string download(const std::string &url, const std::filesystem::path &part, uint64_t expected_size);

void list_indexes(const FetchArguments &opt);

/// Returns the process exit code.
int fetch_index(const FetchArguments &opt);

} // namespace fetch

void setup_fetch_subcommand(CLI::App &app);

#endif // CHARON_FETCH_H
