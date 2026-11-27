#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <thread>

#include <unistd.h>

#include <curl/curl.h>
#include <zlib.h>

#include <cereal/external/rapidjson/document.h>

#include "fetch.hpp"
#include "index.hpp"
#include "md5.hpp"
#include "version.h"

namespace fs = std::filesystem;

namespace fetch {

namespace {

// Any version of each record works as a starting point: fetch follows it to the latest version.
const std::map<std::string, std::string> &records() {
    static const std::map<std::string, std::string> r{
            {"nanopore", "24559698"},  // 10.5281/zenodo.24559698, concept 10.5281/zenodo.15398094
            {"illumina", "24531993"},  // 10.5281/zenodo.24531993
    };
    return r;
}

std::string user_agent() {
    return std::string{"charon/"} + SOFTWARE_VERSION + " (https://github.com/rmcolq/charon)";
}

struct CurlHandle {
    CURL *curl{curl_easy_init()};

    CurlHandle() {
        if (!curl)
            throw std::runtime_error("Could not initialise libcurl");
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        curl_easy_setopt(curl, CURLOPT_FAILONERROR, 1L);
        curl_easy_setopt(curl, CURLOPT_USERAGENT, user_agent().c_str());
        curl_easy_setopt(curl, CURLOPT_CONNECTTIMEOUT, 60L);
        // treat a transfer that stalls (under 1 KB/s for 2 minutes) as failed, so it can be retried
        curl_easy_setopt(curl, CURLOPT_LOW_SPEED_LIMIT, 1024L);
        curl_easy_setopt(curl, CURLOPT_LOW_SPEED_TIME, 120L);
    }

    ~CurlHandle() { curl_easy_cleanup(curl); }

    CurlHandle(const CurlHandle &) = delete;

    CurlHandle &operator=(const CurlHandle &) = delete;
};

long status_of(CURL *curl, CURLcode code) {
    long status{0};
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &status);
    if (status == 0 && code == CURLE_FILE_COULDNT_READ_FILE)
        status = 404;
    return status;
}

[[noreturn]] void throw_request_error(CURL *curl, CURLcode code, const std::string &url, const char *errbuf) {
    const long status = status_of(curl, code);
    std::string what = "Request to " + url + " failed: ";
    what += errbuf[0] ? errbuf : curl_easy_strerror(code);
    throw RequestError(what, status);
}

size_t append_to_string(char *ptr, size_t size, size_t nmemb, void *userdata) {
    static_cast<std::string *>(userdata)->append(ptr, size * nmemb);
    return size * nmemb;
}

std::string get_string(const rapidjson::Value &object, const char *key) {
    if (!object.IsObject())
        return {};
    auto it = object.FindMember(key);
    if (it == object.MemberEnd())
        return {};
    if (it->value.IsString())
        return it->value.GetString();
    if (it->value.IsInt64())
        return std::to_string(it->value.GetInt64());
    return {};
}

} // namespace

const std::vector<std::string> &technologies() {
    static const std::vector<std::string> t{"nanopore", "illumina"};
    return t;
}

std::string record_for(const std::string &tech) {
    auto it = records().find(tech);
    if (it == records().end())
        throw std::runtime_error("No index record is known for technology '" + tech + "'");
    return it->second;
}

std::string api_base(const FetchArguments &opt) {
    std::string base{opt.api_url};
    if (base.empty()) {
        const char *env = std::getenv("CHARON_ZENODO_API");
        base = (env && *env) ? env : std::string{default_api_url};
    }
    while (!base.empty() && base.back() == '/')
        base.pop_back();
    return base;
}

std::string record_url(const std::string &base, const std::string &record_id, bool latest) {
    return base + "/records/" + record_id + (latest ? "/versions/latest" : "");
}

fs::path output_dir(const FetchArguments &opt) {
    if (!opt.outdir.empty())
        return opt.outdir;
    const char *env = std::getenv("CHARON_DB");
    if (env && *env)
        return env;
    return fs::current_path();
}

RecordInfo parse_record(std::string_view json) {
    rapidjson::Document doc;
    doc.Parse(json.data(), json.size());
    if (doc.HasParseError() || !doc.IsObject())
        throw std::runtime_error("Zenodo did not return a valid record (could not parse JSON)");

    RecordInfo record;
    record.id = get_string(doc, "id");
    record.doi = get_string(doc, "doi");
    if (doc.HasMember("metadata")) {
        record.version = get_string(doc["metadata"], "version");
        record.publication_date = get_string(doc["metadata"], "publication_date");
    }
    if (record.id.empty() || !doc.HasMember("files") || !doc["files"].IsArray())
        throw std::runtime_error("Zenodo did not return a valid record (no id or file list)");

    for (const auto &f: doc["files"].GetArray()) {
        RemoteFile file;
        file.key = get_string(f, "key");
        if (f.HasMember("size") && f["size"].IsUint64())
            file.size = f["size"].GetUint64();
        std::string checksum = get_string(f, "checksum");
        if (checksum.rfind("md5:", 0) == 0)
            file.md5 = checksum.substr(4);
        if (f.HasMember("links"))
            file.url = get_string(f["links"], "self");
        if (!file.key.empty())
            record.files.push_back(std::move(file));
    }
    return record;
}

std::string index_file_name(const std::string &tech, const std::string &size, bool host_only) {
    return "charon_" + size + (host_only ? "_host_only" : "") + "." + tech + ".idx.gz";
}

const RemoteFile &select_file(const RecordInfo &record, const std::string &name) {
    for (const auto &f: record.files)
        if (f.key == name) {
            if (f.url.empty() || f.md5.empty())
                throw std::runtime_error("Zenodo lists " + name + " without a download link or md5");
            return f;
        }
    std::string what = "Record " + record.id + " has no file " + name + ". Indexes it has:";
    for (const auto &f: record.files)
        if (f.key.ends_with(".idx.gz"))
            what += "\n  " + f.key + "  (" + selection_flags(f.key) + ")";
    throw std::runtime_error(what);
}

std::string selection_flags(const std::string &name) {
    static constexpr std::string_view prefix{"charon_"}, suffix{".idx.gz"};
    if (!name.starts_with(prefix) || !name.ends_with(suffix))
        return {};
    std::string middle = name.substr(prefix.size(), name.size() - prefix.size() - suffix.size());
    const auto dot = middle.find('.');
    if (dot == std::string::npos)
        return {};
    std::string size = middle.substr(0, dot);
    const std::string tech = middle.substr(dot + 1);
    bool host_only = false;
    if (size.ends_with("_host_only")) {
        host_only = true;
        size.resize(size.size() - std::string_view{"_host_only"}.size());
    }
    return "--tech " + tech + " --size " + size + (host_only ? " --host-only" : "");
}

std::string human_size(uint64_t bytes) {
    static const char *units[] = {"B", "KB", "MB", "GB", "TB"};
    double value = static_cast<double>(bytes);
    int unit = 0;
    while (value >= 1000 && unit < 4) {
        value /= 1000;
        ++unit;
    }
    std::ostringstream out;
    out << std::fixed << std::setprecision(unit == 0 ? 0 : 2) << value << " " << units[unit];
    return out.str();
}

void hash_file_prefix(Md5 &md5, const fs::path &path, uint64_t length) {
    std::ifstream in{path, std::ios::binary};
    if (!in)
        throw std::runtime_error("Could not read " + path.string());
    std::vector<char> buffer(1 << 20);
    while (length > 0) {
        const auto want = static_cast<std::streamsize>(std::min<uint64_t>(buffer.size(), length));
        in.read(buffer.data(), want);
        const auto got = in.gcount();
        if (got <= 0)
            throw std::runtime_error("Unexpected end of " + path.string());
        md5.update(buffer.data(), static_cast<size_t>(got));
        length -= static_cast<uint64_t>(got);
    }
}

std::string file_md5(const fs::path &path, std::optional<uint64_t> length) {
    Md5 md5;
    hash_file_prefix(md5, path, length.value_or(fs::file_size(path)));
    return md5.hexdigest();
}

std::string check_index_start(std::string_view first_bytes) {
    const std::string_view magic{Index::format_magic.data(), Index::format_magic.size()};
    if (first_bytes.substr(0, magic.size()) == magic)
        return std::string{magic.substr(magic.size() - 2)};
    if (first_bytes.starts_with("CHARONv"))
        throw std::runtime_error("This index has format " + std::string{first_bytes.substr(6, 2)} +
                                 ", which Charon " + SOFTWARE_VERSION +
                                 " cannot read. Update Charon, or use --record to choose an older index version.");
    // Indexes written by Charon <= v1.1.4 (format v3) start straight with the archive and have no magic.
    return "v3 or older (no format marker)";
}

namespace {

/// Streaming gzip decompression (multi-member) that hands decompressed bytes to a sink.
class Gunzip {
public:
    Gunzip() {
        if (inflateInit2(&stream_, 15 + 32) != Z_OK)
            throw std::runtime_error("Could not initialise zlib");
    }

    ~Gunzip() { inflateEnd(&stream_); }

    Gunzip(const Gunzip &) = delete;

    Gunzip &operator=(const Gunzip &) = delete;

    /// Returns false once the sink asks to stop.
    template<typename Sink>
    bool feed(const char *data, size_t len, Sink &&sink) {
        stream_.next_in = reinterpret_cast<const Bytef *>(data);
        stream_.avail_in = static_cast<uInt>(len);
        while (stream_.avail_in > 0) {
            if (finished_) {
                // another gzip member follows
                inflateReset(&stream_);
                finished_ = false;
            }
            stream_.next_out = reinterpret_cast<Bytef *>(out_.data());
            stream_.avail_out = static_cast<uInt>(out_.size());
            const int ret = inflate(&stream_, Z_NO_FLUSH);
            if (ret != Z_OK && ret != Z_STREAM_END && ret != Z_BUF_ERROR)
                throw std::runtime_error(std::string{"Corrupt gzip data: "} + (stream_.msg ? stream_.msg : "error"));
            const size_t produced = out_.size() - stream_.avail_out;
            if (produced && !sink(out_.data(), produced))
                return false;
            if (ret == Z_STREAM_END)
                finished_ = true;
            else if (ret == Z_BUF_ERROR && produced == 0)
                break;
        }
        return true;
    }

    /// True if the last gzip member was complete (zlib has checked its CRC and length).
    [[nodiscard]] bool complete() const { return finished_; }

private:
    z_stream stream_{};
    std::vector<char> out_ = std::vector<char>(1 << 20);
    bool finished_{false};
};

} // namespace

void gunzip_index(const fs::path &gz, const fs::path &out) {
    std::ifstream in{gz, std::ios::binary};
    if (!in)
        throw std::runtime_error("Could not read " + gz.string());
    std::ofstream os{out, std::ios::binary | std::ios::trunc};
    if (!os)
        throw std::runtime_error("Could not write " + out.string());

    Gunzip gunzip;
    std::string head;
    bool checked = false;
    auto sink = [&](const char *data, size_t len) {
        if (!checked) {
            head.append(data, std::min(len, Index::format_magic.size() - head.size()));
            if (head.size() == Index::format_magic.size()) {
                check_index_start(head);
                checked = true;
            }
        }
        os.write(data, static_cast<std::streamsize>(len));
        if (!os)
            throw std::runtime_error("Could not write " + out.string() + " (disk full?)");
        return true;
    };

    std::vector<char> buffer(1 << 20);
    while (in) {
        in.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        if (in.gcount() > 0)
            gunzip.feed(buffer.data(), static_cast<size_t>(in.gcount()), sink);
    }
    if (!gunzip.complete())
        throw std::runtime_error(gz.string() + " ends before the end of its gzip data (truncated?)");
    if (!checked)
        check_index_start(head);
    os.close();
    if (!os)
        throw std::runtime_error("Could not write " + out.string());
}

fs::path source_path(const fs::path &index) {
    return fs::path{index.string() + ".source.tsv"};
}

void write_source(const fs::path &index, const Source &source) {
    std::ofstream os{source_path(index)};
    os << "file\t" << source.file << "\n"
       << "record\t" << source.record << "\n"
       << "doi\t" << source.doi << "\n"
       << "version\t" << source.version << "\n"
       << "md5_gz\t" << source.md5 << "\n";
    if (!os)
        throw std::runtime_error("Could not write " + source_path(index).string());
}

std::optional<Source> read_source(const fs::path &index) {
    std::ifstream in{source_path(index)};
    if (!in)
        return std::nullopt;
    Source source;
    std::string line;
    while (std::getline(in, line)) {
        const auto tab = line.find('\t');
        if (tab == std::string::npos)
            continue;
        const std::string key = line.substr(0, tab), value = line.substr(tab + 1);
        if (key == "file") source.file = value;
        else if (key == "record") source.record = value;
        else if (key == "doi") source.doi = value;
        else if (key == "version") source.version = value;
        else if (key == "md5_gz") source.md5 = value;
    }
    if (source.md5.empty())
        return std::nullopt;
    return source;
}

std::string http_get(const std::string &url) {
    CurlHandle handle;
    std::string body;
    char errbuf[CURL_ERROR_SIZE] = "";
    curl_easy_setopt(handle.curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(handle.curl, CURLOPT_ERRORBUFFER, errbuf);
    curl_easy_setopt(handle.curl, CURLOPT_WRITEFUNCTION, append_to_string);
    curl_easy_setopt(handle.curl, CURLOPT_WRITEDATA, &body);
    const CURLcode code = curl_easy_perform(handle.curl);
    if (code != CURLE_OK)
        throw_request_error(handle.curl, code, url, errbuf);
    return body;
}

std::string check_remote_index(const std::string &url) {
    struct State {
        Gunzip gunzip;
        std::string head;
        std::string error;
    } state;

    auto write = [](char *ptr, size_t size, size_t nmemb, void *userdata) -> size_t {
        auto &s = *static_cast<State *>(userdata);
        try {
            s.gunzip.feed(ptr, size * nmemb, [&](const char *data, size_t len) {
                s.head.append(data, std::min(len, Index::format_magic.size() - s.head.size()));
                return s.head.size() < Index::format_magic.size();
            });
        } catch (const std::exception &e) {
            s.error = e.what();
            return 0;
        }
        // stop as soon as the format marker has been read
        return s.head.size() < Index::format_magic.size() ? size * nmemb : 0;
    };

    CurlHandle handle;
    char errbuf[CURL_ERROR_SIZE] = "";
    curl_easy_setopt(handle.curl, CURLOPT_URL, url.c_str());
    curl_easy_setopt(handle.curl, CURLOPT_ERRORBUFFER, errbuf);
    curl_easy_setopt(handle.curl, CURLOPT_RANGE, "0-65535");
    curl_easy_setopt(handle.curl, CURLOPT_WRITEFUNCTION, static_cast<curl_write_callback>(write));
    curl_easy_setopt(handle.curl, CURLOPT_WRITEDATA, &state);
    const CURLcode code = curl_easy_perform(handle.curl);
    if (!state.error.empty())
        throw std::runtime_error("Could not read the start of " + url + ": " + state.error);
    if (code != CURLE_OK && code != CURLE_WRITE_ERROR)
        throw_request_error(handle.curl, code, url, errbuf);
    return check_index_start(state.head);
}

namespace {

struct Transfer {
    std::FILE *file{nullptr};
    Md5 *md5{nullptr};
    uint64_t written{0};       // bytes in the part file and in md5
    uint64_t expected{0};
    bool tty{false};
    int last_decile{-1};
    std::chrono::steady_clock::time_point last_print{};
    std::chrono::steady_clock::time_point start{std::chrono::steady_clock::now()};
    uint64_t start_bytes{0};
};

size_t write_part(char *ptr, size_t size, size_t nmemb, void *userdata) {
    auto &t = *static_cast<Transfer *>(userdata);
    const size_t len = size * nmemb;
    if (std::fwrite(ptr, 1, len, t.file) != len)
        return 0;
    t.md5->update(ptr, len);
    t.written += len;
    return len;
}

void print_progress(Transfer &t, bool final) {
    const auto now = std::chrono::steady_clock::now();
    const double pct = t.expected ? 100.0 * static_cast<double>(t.written) / static_cast<double>(t.expected) : 0;
    const double secs = std::chrono::duration<double>(now - t.start).count();
    const double rate = secs > 0 ? static_cast<double>(t.written - t.start_bytes) / secs : 0;
    std::ostringstream line;
    line << "  " << human_size(t.written) << " / " << human_size(t.expected) << " (" << std::fixed
         << std::setprecision(1) << pct << "%) " << human_size(static_cast<uint64_t>(rate)) << "/s";
    if (t.tty) {
        if (final || now - t.last_print > std::chrono::milliseconds(500)) {
            std::cerr << "\r" << line.str() << "    " << (final ? "\n" : "") << std::flush;
            t.last_print = now;
        }
    } else {
        const int decile = static_cast<int>(pct / 10);
        if (decile > t.last_decile || (final && t.last_decile < 10)) {
            std::cerr << line.str() << "\n";
            t.last_decile = decile;
        }
    }
}

int progress(void *userdata, curl_off_t, curl_off_t, curl_off_t, curl_off_t) {
    print_progress(*static_cast<Transfer *>(userdata), false);
    return 0;
}

} // namespace

std::string download(const std::string &url, const fs::path &part, uint64_t expected_size) {
    uint64_t have = fs::exists(part) ? fs::file_size(part) : 0;
    if (have > expected_size) {
        std::cerr << "Discarding " << part.string() << ": it is larger than the file on Zenodo\n";
        fs::remove(part);
        have = 0;
    }
    Md5 md5;
    if (have > 0) {
        std::cerr << "Resuming from " << human_size(have) << " already downloaded\n";
        hash_file_prefix(md5, part, have);
    }

    Transfer t;
    t.md5 = &md5;
    t.written = have;
    t.expected = expected_size;
    t.tty = isatty(fileno(stderr));
    t.start_bytes = have;

    constexpr int max_attempts = 10;
    for (int attempt = 1; t.written < expected_size; ++attempt) {
        t.file = std::fopen(part.c_str(), "ab");
        if (!t.file)
            throw std::runtime_error("Could not write " + part.string());

        CurlHandle handle;
        char errbuf[CURL_ERROR_SIZE] = "";
        curl_easy_setopt(handle.curl, CURLOPT_URL, url.c_str());
        curl_easy_setopt(handle.curl, CURLOPT_ERRORBUFFER, errbuf);
        curl_easy_setopt(handle.curl, CURLOPT_RESUME_FROM_LARGE, static_cast<curl_off_t>(t.written));
        curl_easy_setopt(handle.curl, CURLOPT_WRITEFUNCTION, write_part);
        curl_easy_setopt(handle.curl, CURLOPT_WRITEDATA, &t);
        curl_easy_setopt(handle.curl, CURLOPT_NOPROGRESS, 0L);
        curl_easy_setopt(handle.curl, CURLOPT_XFERINFOFUNCTION, progress);
        curl_easy_setopt(handle.curl, CURLOPT_XFERINFODATA, &t);
        const CURLcode code = curl_easy_perform(handle.curl);
        const bool closed = std::fclose(t.file) == 0;
        t.file = nullptr;
        if (!closed)
            throw std::runtime_error("Could not write " + part.string() + " (disk full?)");

        if (code == CURLE_OK)
            break;
        if (code == CURLE_WRITE_ERROR)
            throw std::runtime_error("Could not write " + part.string() + " (disk full?)");
        if (code == CURLE_RANGE_ERROR) {
            // the server would not resume: start again
            std::cerr << "\nThe server cannot resume this download; starting again\n";
            fs::remove(part);
            md5 = Md5{};
            t.written = 0;
            t.start_bytes = 0;
        }
        const long status = status_of(handle.curl, code);
        const bool permanent = status >= 400 && status < 500 && status != 408 && status != 429;
        if (permanent || attempt == max_attempts)
            throw_request_error(handle.curl, code, url, errbuf);
        // keep the part file consistent with what has been hashed
        if (fs::exists(part) && fs::file_size(part) != t.written)
            fs::resize_file(part, t.written);
        const auto wait = std::chrono::seconds(std::min(60, 5 * attempt));
        std::cerr << "\nDownload interrupted (" << (errbuf[0] ? errbuf : curl_easy_strerror(code))
                  << "); retrying in " << wait.count() << " s (attempt " << attempt + 1 << " of " << max_attempts
                  << ")\n";
        std::this_thread::sleep_for(wait);
    }
    print_progress(t, true);
    if (t.written != expected_size)
        throw std::runtime_error("Downloaded " + std::to_string(t.written) + " bytes but Zenodo lists " +
                                 std::to_string(expected_size));
    return md5.hexdigest();
}

namespace {

/// Fetch a record, following it to its latest version unless `pinned`.
RecordInfo get_record(const std::string &base, const std::string &record_id, bool pinned, const std::string &tech) {
    const std::string url = record_url(base, record_id, !pinned);
    try {
        return parse_record(http_get(url));
    } catch (const RequestError &e) {
        if (e.status == 404 || e.status == 410)
            throw std::runtime_error("The " + tech + " index record (" + record_id +
                                     ") is not available on Zenodo: it may not be published yet.\n  " + url);
        throw;
    }
}

void print_record(const RecordInfo &record) {
    std::cout << "Zenodo record " << record.id;
    if (!record.version.empty())
        std::cout << ", version " << record.version;
    if (!record.publication_date.empty())
        std::cout << ", published " << record.publication_date;
    if (!record.doi.empty())
        std::cout << " (doi:" << record.doi << ")";
    std::cout << "\n";
}

} // namespace

void list_indexes(const FetchArguments &opt) {
    const std::string base = api_base(opt);
    std::vector<std::string> techs;
    if (opt.tech.empty())
        techs = technologies();
    else
        techs.push_back(opt.tech);

    for (const auto &tech: techs) {
        std::cout << tech << ": ";
        RecordInfo record;
        try {
            record = get_record(base, opt.record.empty() ? record_for(tech) : opt.record, !opt.record.empty(), tech);
        } catch (const std::runtime_error &e) {
            std::cout << e.what() << "\n\n";
            continue;
        }
        print_record(record);
        for (const auto &f: record.files) {
            if (!f.key.ends_with(".idx.gz"))
                continue;
            std::cout << "  " << std::left << std::setw(40) << f.key << std::right << std::setw(10)
                      << human_size(f.size) << "  " << selection_flags(f.key) << "\n";
        }
        std::cout << "\n";
    }
    std::cout << "Sizes are compressed; the decompressed index is larger. Download one with: charon fetch --tech <tech> "
                 "[--size small|big] [--host-only]\n";
}

int fetch_index(const FetchArguments &opt) {
    const std::string base = api_base(opt);
    const std::string record_id = opt.record.empty() ? record_for(opt.tech) : opt.record;
    const RecordInfo record = get_record(base, record_id, !opt.record.empty(), opt.tech);
    const std::string name = index_file_name(opt.tech, opt.size, opt.host_only);
    const RemoteFile &remote = select_file(record, name);

    const fs::path dir = output_dir(opt);
    fs::create_directories(dir);
    const fs::path gz = dir / remote.key;
    const fs::path index = dir / remote.key.substr(0, remote.key.size() - 3);  // strip .gz
    const fs::path part = fs::path{gz.string() + "." + remote.md5.substr(0, 8) + ".part"};
    const fs::path tmp = fs::path{index.string() + ".tmp"};

    std::cout << "Index:  " << remote.key << " (" << opt.tech << ", " << opt.size
              << (opt.host_only ? ", host only" : "") << ")\n"
              << "From:   ";
    print_record(record);
    std::cout << "Size:   " << human_size(remote.size) << " to download (gzip); md5 " << remote.md5 << "\n"
              << "To:     " << index.string() << "\n";

    if (fs::exists(index) && !opt.force) {
        const auto source = read_source(index);
        if (source && source->md5 == remote.md5) {
            std::cout << index.string() << " is already up to date.\n";
            return 0;
        }
        if (source)
            std::cerr << index.string() << " exists, from record " << source->record << " (version " << source->version
                      << "). Record " << record.id << " has a different version of it: rerun with --force to replace "
                      << "it.\n";
        else
            std::cerr << index.string() << " exists, but there is no " << source_path(index).filename().string()
                      << " to say where it came from. Rerun with --force to replace it.\n";
        return 1;
    }

    std::cout << "Checking the index format... " << std::flush;
    std::cout << check_remote_index(remote.url) << ", which this Charon can read\n";

    const bool have_gz = fs::exists(gz) && fs::file_size(gz) == remote.size && file_md5(gz) == remote.md5;
    if (have_gz) {
        std::cout << "Using " << gz.string() << ", which is already downloaded\n";
    } else {
        const uint64_t have = fs::exists(part) ? fs::file_size(part) : 0;
        const uint64_t remaining = remote.size - std::min(have, remote.size);
        const auto space = fs::space(dir);
        std::cout << "Downloading " << human_size(remaining) << (have ? " (resuming)" : "") << "; "
                  << human_size(space.available) << " free in " << dir.string() << "\n";
        // the decompressed index is a little larger than the .gz, and both exist while decompressing
        if (space.available < remaining + remote.size)
            std::cerr << "Warning: there may not be enough disk space. Downloading and decompressing need about "
                      << human_size(remaining + remote.size) << " more than the decompressed index needs alone; "
                      << "use -o to choose another directory.\n";

        const std::string md5 = download(remote.url, part, remote.size);
        if (md5 != remote.md5) {
            fs::remove(part);
            throw std::runtime_error("md5 of the download (" + md5 + ") does not match Zenodo's (" + remote.md5 +
                                     "). The partial file has been removed; please try again.");
        }
        std::cout << "md5 matches Zenodo\n";
        fs::rename(part, gz);
    }

    std::cout << "Decompressing to " << index.string() << "\n";
    try {
        gunzip_index(gz, tmp);
    } catch (...) {
        fs::remove(tmp);
        throw;
    }
    fs::rename(tmp, index);
    write_source(index, Source{remote.key, record.id, record.doi, record.version, remote.md5});
    if (!opt.keep_gz)
        fs::remove(gz);

    std::cout << "Done: " << index.string() << " (" << human_size(fs::file_size(index)) << ")\n";
    return 0;
}

} // namespace fetch

void setup_fetch_subcommand(CLI::App &app) {
    auto opt = std::make_shared<FetchArguments>();
    auto *sub = app.add_subcommand("fetch", "Download a published Charon index from Zenodo.");

    sub->add_option("--tech", opt->tech, "Read technology the index is for (required unless --list).")
            ->check(CLI::IsMember(fetch::technologies()))
            ->type_name("STRING");
    sub->add_option("--size", opt->size, "Index size: big (200 human assemblies) or small (4).")
            ->check(CLI::IsMember({"small", "big"}))
            ->type_name("STRING")
            ->capture_default_str();
    sub->add_flag("--host-only", opt->host_only, "The index with only the human category (big only), for dehost.");
    sub->add_option("-o,--outdir", opt->outdir, "Directory for the index [default: $CHARON_DB if set, else the "
                                                "current directory].")
            ->type_name("DIR");
    sub->add_flag("--list", opt->list, "List the indexes available to fetch, with their sizes, and exit.");
    sub->add_option("--record", opt->record, "Use this Zenodo record (version) instead of the latest.")
            ->type_name("ID");
    sub->add_flag("--keep-gz", opt->keep_gz, "Keep the downloaded .gz file after decompressing it.");
    sub->add_flag("--force", opt->force, "Replace an existing index of the same name.");
    sub->add_option("--api-url", opt->api_url, "Zenodo API base url [default: $CHARON_ZENODO_API if set, else "
                                               "https://zenodo.org/api].")
            ->group("");  // hidden: for testing and mirrors

    sub->callback([opt]() {
        if (!opt->list && opt->tech.empty())
            throw CLI::RequiredError("--tech (or use --list)");
        if (opt->host_only && opt->size == "small")
            throw CLI::ValidationError("--host-only", "there is no small host-only index; use --size big");

        curl_global_init(CURL_GLOBAL_DEFAULT);
        int code = 0;
        try {
            if (opt->list)
                fetch::list_indexes(*opt);
            else
                code = fetch::fetch_index(*opt);
        } catch (const std::exception &e) {
            std::cerr << "\nError: " << e.what() << "\n";
            code = 1;
        }
        curl_global_cleanup();
        if (code != 0)
            std::exit(code);
    });
}
