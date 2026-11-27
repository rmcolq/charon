// Tests for `charon fetch`. Everything is served from file:// urls in a temporary directory, so no network
// access is needed: a fake index (the v4 format marker followed by random bytes) is gzipped and described by a
// Zenodo-style record JSON.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

#include <curl/curl.h>
#include <zlib.h>

#include "fetch.hpp"
#include "md5.hpp"

namespace fs = std::filesystem;
using Catch::Matchers::ContainsSubstring;

namespace {

std::string md5_of(std::string_view s) {
    Md5 md5;
    md5.update(s);
    return md5.hexdigest();
}

std::string read_file(const fs::path &p) {
    std::ifstream in{p, std::ios::binary};
    return {std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
}

void write_file(const fs::path &p, std::string_view content) {
    fs::create_directories(p.parent_path());
    std::ofstream os{p, std::ios::binary};
    os.write(content.data(), static_cast<std::streamsize>(content.size()));
}

std::string gzip(std::string_view content, const fs::path &scratch) {
    gzFile f = gzopen(scratch.c_str(), "wb1");
    gzwrite(f, content.data(), static_cast<unsigned>(content.size()));
    gzclose(f);
    std::string out = read_file(scratch);
    fs::remove(scratch);
    return out;
}

std::string random_bytes(size_t n, unsigned seed) {
    std::mt19937 rng{seed};
    std::string s(n, '\0');
    for (auto &c: s)
        c = static_cast<char>(rng() % 7);  // compressible, like an index
    return s;
}

std::string file_url(const fs::path &p) { return "file://" + p.string(); }

struct TempDir {
    fs::path path;

    TempDir() {
        static int counter = 0;
        path = fs::temp_directory_path() / ("charon_test_fetch_" + std::to_string(::getpid()) + "_" +
                                            std::to_string(counter++));
        fs::remove_all(path);
        fs::create_directories(path);
    }

    ~TempDir() { fs::remove_all(path); }
};

/// A fake Zenodo API under <dir>/api with one nanopore record holding charon_big.nanopore.idx.gz
/// (and charon_big_host_only.nanopore.idx.gz with a wrong md5).
struct FakeZenodo {
    TempDir dir;
    std::string index;  // the uncompressed index
    std::string gz;     // what is served

    FakeZenodo() {
        index = std::string{"CHARONv4"} + random_bytes(1'500'000, 1);
        gz = gzip(index, dir.path / "tmp.gz");
        const fs::path files = dir.path / "files";
        write_file(files / "charon_big.nanopore.idx.gz", gz);
        write_file(files / "charon_big_host_only.nanopore.idx.gz", gz);

        const std::string json = record_json("24559698", "v2", {
                {"charon_big.nanopore.idx.gz", md5_of(gz)},
                {"charon_big_host_only.nanopore.idx.gz", std::string(32, '0')},
        });
        write_file(dir.path / "api/records/24559698/versions/latest", json);
        write_file(dir.path / "api/records/24559698", json);
    }

    std::string record_json(const std::string &id, const std::string &version,
                            const std::vector<std::pair<std::string, std::string>> &entries) const {
        std::string files;
        for (const auto &[name, md5]: entries) {
            if (!files.empty())
                files += ",";
            files += R"({"key": ")" + name + R"(", "size": )" + std::to_string(gz.size()) +
                     R"(, "checksum": "md5:)" + md5 + R"(", "links": {"self": ")" +
                     file_url(dir.path / "files" / name) + R"("}})";
        }
        files += R"(, {"key": "README.md", "size": 5, "checksum": "md5:00", "links": {"self": "x"}})";
        return R"({"id": )" + id + R"(, "doi": "10.5281/zenodo.)" + id + R"(", "metadata": {"version": ")" +
               version + R"(", "publication_date": "2026-11-26"}, "files": [)" + files + "]}";
    }

    [[nodiscard]] std::string api() const { return file_url(dir.path / "api"); }
};

struct Run {
    int exit_code;
    std::string output;
};

Run run_charon(const std::string &args) {
    const std::string command = std::string{CHARON_BIN} + " " + args + " 2>&1";
    Run result{-1, {}};
    FILE *pipe = popen(command.c_str(), "r");
    if (!pipe)
        return result;
    char buffer[512];
    while (fgets(buffer, sizeof(buffer), pipe))
        result.output += buffer;
    const int status = pclose(pipe);
    result.exit_code = WEXITSTATUS(status);
    return result;
}

struct CurlGlobal {
    CurlGlobal() { curl_global_init(CURL_GLOBAL_DEFAULT); }

    ~CurlGlobal() { curl_global_cleanup(); }
} curl_global;

} // namespace

TEST_CASE("md5 matches the RFC 1321 test vectors", "[fetch][md5]") {
    CHECK(md5_of("") == "d41d8cd98f00b204e9800998ecf8427e");
    CHECK(md5_of("abc") == "900150983cd24fb0d6963f7d28e17f72");
    CHECK(md5_of("message digest") == "f96b697d7cb7938d525a2f31aaf161d0");
    CHECK(md5_of("12345678901234567890123456789012345678901234567890123456789012345678901234567890") ==
          "57edf4a22be3c955ac49da2e2107b67a");
}

TEST_CASE("md5 does not depend on how the input is split", "[fetch][md5]") {
    const std::string data = random_bytes(10'000, 2);
    const std::string whole = md5_of(data);
    for (size_t chunk: {1u, 3u, 63u, 64u, 65u, 1000u}) {
        Md5 md5;
        for (size_t i = 0; i < data.size(); i += chunk)
            md5.update(std::string_view{data}.substr(i, chunk));
        CHECK(md5.hexdigest() == whole);
    }
}

TEST_CASE("file_md5 hashes a whole file or a prefix", "[fetch][md5]") {
    TempDir dir;
    write_file(dir.path / "f", "abcdef");
    CHECK(fetch::file_md5(dir.path / "f") == md5_of("abcdef"));
    CHECK(fetch::file_md5(dir.path / "f", 3) == md5_of("abc"));
}

TEST_CASE("index file names and the options that select them", "[fetch]") {
    CHECK(fetch::index_file_name("nanopore", "big", false) == "charon_big.nanopore.idx.gz");
    CHECK(fetch::index_file_name("illumina", "big", true) == "charon_big_host_only.illumina.idx.gz");
    CHECK(fetch::selection_flags("charon_big_host_only.illumina.idx.gz") ==
          "--tech illumina --size big --host-only");
    CHECK(fetch::selection_flags("charon_small.nanopore.idx.gz") == "--tech nanopore --size small");
    CHECK(fetch::selection_flags("README.md").empty());
}

TEST_CASE("human_size", "[fetch]") {
    CHECK(fetch::human_size(512) == "512 B");
    CHECK(fetch::human_size(616661317) == "616.66 MB");
    CHECK(fetch::human_size(6255801888) == "6.26 GB");
}

TEST_CASE("api base and record urls", "[fetch]") {
    FetchArguments opt;
    opt.api_url = "https://example.org/api/";
    CHECK(fetch::api_base(opt) == "https://example.org/api");
    CHECK(fetch::record_url("https://zenodo.org/api", "24559698", true) ==
          "https://zenodo.org/api/records/24559698/versions/latest");
    CHECK(fetch::record_url("https://zenodo.org/api", "15398095", false) ==
          "https://zenodo.org/api/records/15398095");
    CHECK_FALSE(fetch::record_for("nanopore").empty());
    CHECK_FALSE(fetch::record_for("illumina").empty());
    CHECK_THROWS(fetch::record_for("pacbio"));
}

TEST_CASE("parse a Zenodo record and select a file", "[fetch]") {
    const std::string json = R"({"id": 24559698, "doi": "10.5281/zenodo.24559698",
        "metadata": {"version": "v2", "publication_date": "2026-11-26"},
        "files": [
          {"key": "charon_small.nanopore.idx.gz", "size": 1251009181,
           "checksum": "md5:c114b1b6a425ed581843766f39e9a618",
           "links": {"self": "https://zenodo.org/api/records/24559698/files/charon_small.nanopore.idx.gz/content"}},
          {"key": "README.md", "size": 6357, "checksum": "md5:f9fb6c1e1cf10709be93c9e64953c202",
           "links": {"self": "https://zenodo.org/x"}}]})";
    const auto record = fetch::parse_record(json);
    CHECK(record.id == "24559698");
    CHECK(record.doi == "10.5281/zenodo.24559698");
    CHECK(record.version == "v2");
    REQUIRE(record.files.size() == 2);

    const auto &f = fetch::select_file(record, "charon_small.nanopore.idx.gz");
    CHECK(f.size == 1251009181);
    CHECK(f.md5 == "c114b1b6a425ed581843766f39e9a618");
    CHECK_THAT(f.url, ContainsSubstring("/content"));

    try {
        fetch::select_file(record, "charon_small_host_only.nanopore.idx.gz");
        FAIL("expected an error");
    } catch (const std::runtime_error &e) {
        CHECK_THAT(e.what(), ContainsSubstring("charon_small.nanopore.idx.gz  (--tech nanopore --size small)"));
        CHECK_THAT(e.what(), !ContainsSubstring("README"));
    }

    CHECK_THROWS(fetch::parse_record("not json"));
    CHECK_THROWS(fetch::parse_record(R"({"status": 404, "message": "not found"})"));
}

TEST_CASE("the index format marker is checked", "[fetch]") {
    CHECK(fetch::check_index_start("CHARONv4") == "v4");
    CHECK_THROWS_WITH(fetch::check_index_start("CHARONv5"), ContainsSubstring("cannot read"));
    CHECK_THAT(fetch::check_index_start(std::string(8, '\x01')), ContainsSubstring("v3"));
}

TEST_CASE("gunzip_index", "[fetch][gzip]") {
    TempDir dir;
    const std::string index = std::string{"CHARONv4"} + random_bytes(2'500'000, 3);

    SECTION("round trip") {
        write_file(dir.path / "a.gz", gzip(index, dir.path / "t.gz"));
        fetch::gunzip_index(dir.path / "a.gz", dir.path / "a.idx");
        CHECK(read_file(dir.path / "a.idx") == index);
    }
    SECTION("several gzip members") {
        const std::string half1 = index.substr(0, 1'000'000), half2 = index.substr(1'000'000);
        write_file(dir.path / "a.gz", gzip(half1, dir.path / "t.gz") + gzip(half2, dir.path / "t.gz"));
        fetch::gunzip_index(dir.path / "a.gz", dir.path / "a.idx");
        CHECK(read_file(dir.path / "a.idx") == index);
    }
    SECTION("a truncated file is an error") {
        const std::string gz = gzip(index, dir.path / "t.gz");
        write_file(dir.path / "a.gz", gz.substr(0, gz.size() / 2));
        CHECK_THROWS_WITH(fetch::gunzip_index(dir.path / "a.gz", dir.path / "a.idx"),
                          ContainsSubstring("truncated"));
    }
    SECTION("an index format this Charon cannot read is an error") {
        write_file(dir.path / "a.gz", gzip("CHARONv9" + index.substr(8), dir.path / "t.gz"));
        CHECK_THROWS_WITH(fetch::gunzip_index(dir.path / "a.gz", dir.path / "a.idx"),
                          ContainsSubstring("cannot read"));
    }
}

TEST_CASE("source file round trip", "[fetch]") {
    TempDir dir;
    const fs::path index = dir.path / "x.idx";
    CHECK_FALSE(fetch::read_source(index).has_value());
    fetch::write_source(index, {"x.idx.gz", "24559698", "10.5281/zenodo.24559698", "v2", "abc"});
    const auto source = fetch::read_source(index);
    REQUIRE(source.has_value());
    CHECK(source->record == "24559698");
    CHECK(source->version == "v2");
    CHECK(source->md5 == "abc");
}

TEST_CASE("requests over file:// urls", "[fetch][download]") {
    FakeZenodo zenodo;
    const fs::path served = zenodo.dir.path / "files/charon_big.nanopore.idx.gz";
    const fs::path part = zenodo.dir.path / "out/x.part";
    fs::create_directories(part.parent_path());

    SECTION("a missing file is a 404") {
        try {
            fetch::http_get(file_url(zenodo.dir.path / "nothing"));
            FAIL("expected an error");
        } catch (const fetch::RequestError &e) {
            CHECK(e.status == 404);
        }
    }
    SECTION("the format marker is read from the start of the remote file") {
        CHECK(fetch::check_remote_index(file_url(served)) == "v4");
    }
    SECTION("full download") {
        CHECK(fetch::download(file_url(served), part, zenodo.gz.size()) == md5_of(zenodo.gz));
        CHECK(read_file(part) == zenodo.gz);
    }
    SECTION("resume a partial download") {
        write_file(part, zenodo.gz.substr(0, zenodo.gz.size() / 3));
        CHECK(fetch::download(file_url(served), part, zenodo.gz.size()) == md5_of(zenodo.gz));
        CHECK(read_file(part) == zenodo.gz);
    }
    SECTION("a part file larger than the remote file is discarded") {
        write_file(part, zenodo.gz + "extra");
        CHECK(fetch::download(file_url(served), part, zenodo.gz.size()) == md5_of(zenodo.gz));
        CHECK(read_file(part) == zenodo.gz);
    }
}

TEST_CASE("charon fetch from a fake Zenodo", "[fetch][cli]") {
    FakeZenodo zenodo;
    const fs::path out = zenodo.dir.path / "out";
    const std::string common = "--api-url " + zenodo.api() + " -o " + out.string();
    const fs::path index = out / "charon_big.nanopore.idx";

    SECTION("download, then up to date, then a different version") {
        auto run = run_charon("fetch --tech nanopore " + common);
        INFO(run.output);
        REQUIRE(run.exit_code == 0);
        // the size is reported before the download starts
        const auto size_line = run.output.find("Size:");
        REQUIRE(size_line != std::string::npos);
        CHECK(size_line < run.output.find("Downloading"));
        CHECK(read_file(index) == zenodo.index);
        CHECK_FALSE(fs::exists(out / "charon_big.nanopore.idx.gz"));
        REQUIRE(fetch::read_source(index).has_value());
        CHECK(fetch::read_source(index)->md5 == md5_of(zenodo.gz));

        run = run_charon("fetch --tech nanopore " + common);
        CHECK(run.exit_code == 0);
        CHECK_THAT(run.output, ContainsSubstring("already up to date"));

        // a new version of the record replaces the file
        const std::string newer = "CHARONv4" + random_bytes(100'000, 9);
        zenodo.gz = gzip(newer, zenodo.dir.path / "t.gz");
        write_file(zenodo.dir.path / "files/charon_big.nanopore.idx.gz", zenodo.gz);
        write_file(zenodo.dir.path / "api/records/24559698/versions/latest",
                   zenodo.record_json("24559699", "v3", {{"charon_big.nanopore.idx.gz", md5_of(zenodo.gz)}}));
        run = run_charon("fetch --tech nanopore " + common);
        CHECK(run.exit_code == 1);
        CHECK_THAT(run.output, ContainsSubstring("--force"));
        run = run_charon("fetch --tech nanopore --force " + common);
        INFO(run.output);
        CHECK(run.exit_code == 0);
        CHECK(read_file(index) == newer);
        CHECK(fetch::read_source(index)->record == "24559699");
    }
    SECTION("--keep-gz keeps the download") {
        REQUIRE(run_charon("fetch --tech nanopore --keep-gz " + common).exit_code == 0);
        CHECK(read_file(out / "charon_big.nanopore.idx.gz") == zenodo.gz);
    }
    SECTION("an md5 mismatch is an error and leaves no files") {
        auto run = run_charon("fetch --tech nanopore --host-only " + common);
        CHECK(run.exit_code == 1);
        CHECK_THAT(run.output, ContainsSubstring("does not match"));
        for (const auto &entry: fs::directory_iterator(out))
            CHECK_THAT(entry.path().filename().string(), !ContainsSubstring("host_only"));
    }
    SECTION("an index the record does not have") {
        auto run = run_charon("fetch --tech nanopore --size small " + common);
        CHECK(run.exit_code == 1);
        CHECK_THAT(run.output, ContainsSubstring("--tech nanopore --size big --host-only"));
    }
    SECTION("an unpublished record") {
        auto run = run_charon("fetch --tech illumina " + common);
        CHECK(run.exit_code == 1);
        CHECK_THAT(run.output, ContainsSubstring("not available on Zenodo"));
    }
    SECTION("there is no small host-only index") {
        auto run = run_charon("fetch --tech nanopore --size small --host-only " + common);
        CHECK(run.exit_code != 0);
        CHECK_THAT(run.output, ContainsSubstring("no small host-only index"));
    }
    SECTION("--tech defaults to nanopore") {
        auto run = run_charon("fetch " + common);
        INFO(run.output);
        CHECK(run.exit_code == 0);
        CHECK(read_file(index) == zenodo.index);
    }
    SECTION("--list --tech lists one technology") {
        auto run = run_charon("fetch --list --tech nanopore --api-url " + zenodo.api());
        CHECK(run.exit_code == 0);
        CHECK_THAT(run.output, ContainsSubstring("charon_big.nanopore.idx.gz"));
        CHECK_THAT(run.output, !ContainsSubstring("illumina:"));
    }
    SECTION("--list") {
        auto run = run_charon("fetch --list --api-url " + zenodo.api());
        CHECK(run.exit_code == 0);
        CHECK_THAT(run.output, ContainsSubstring("charon_big.nanopore.idx.gz"));
        CHECK_THAT(run.output, ContainsSubstring("--tech nanopore --size big"));
        CHECK_THAT(run.output, !ContainsSubstring("README"));
        CHECK_THAT(run.output, ContainsSubstring("illumina: The illumina index record"));
    }
}
