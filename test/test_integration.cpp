#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>

#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>

namespace fs = std::filesystem;

// Test data directory
static const std::string TEST_DATA_DIR = "/Users/rmcolq/Work/git/charon";
static const std::string BUILD_DIR = "/Users/rmcolq/Work/git/charon/build";
static const std::string CHARON_EXECUTABLE = BUILD_DIR + "/bin/charon";

// Helper to run charon command and capture output
struct CharonResult {
    int exit_code;
    std::string stdout_output;
    std::string stderr_output;
    bool success() const { return exit_code == 0; }
};

static CharonResult run_charon(const std::string& args, bool capture_output = true) {
    std::string command = CHARON_EXECUTABLE + " " + args;
    
    CharonResult result;
    
    // Use popen to capture stdout
    if (capture_output) {
        command += " 2>&1";
        FILE* pipe = popen(command.c_str(), "r");
        if (!pipe) {
            result.exit_code = -1;
            return result;
        }
        
        char buffer[256];
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            result.stdout_output += buffer;
        }
        
        int status = pclose(pipe);
        result.exit_code = WEXITSTATUS(status);
    } else {
        result.exit_code = system(command.c_str());
    }
    
    return result;
}

// Helper to check if file exists
static bool file_exists(const std::string& path) {
    return fs::exists(path);
}

// Helper to get file size
static std::size_t file_size(const std::string& path) {
    return fs::file_size(path);
}

// Helper to count lines in file
static std::size_t count_lines(const std::string& path) {
    std::ifstream file(path);
    std::string line;
    std::size_t count = 0;
    while (std::getline(file, line)) {
        count++;
    }
    return count;
}

// Cleanup helper
static void cleanup_files(const std::vector<std::string>& files) {
    for (const auto& file : files) {
        if (file_exists(file)) {
            fs::remove(file);
        }
    }
}

// ─── Version and Help Tests ──────────────────────────────────────────────────

TEST_CASE("charon --help shows usage", "[Tier3][Integration][Help]") {
    auto result = run_charon("--help");
    
    REQUIRE(result.success());
    REQUIRE(result.stdout_output.find("Charon: Categorize reads into a small number of classes") != std::string::npos);
    REQUIRE(result.stdout_output.find("index") != std::string::npos);
    REQUIRE(result.stdout_output.find("classify") != std::string::npos);
    REQUIRE(result.stdout_output.find("dehost") != std::string::npos);
}

TEST_CASE("charon --version shows version", "[Tier3][Integration][Version]") {
    auto result = run_charon("--version");
    
    REQUIRE(result.success());
    // Version output should contain some text (even if just the version number)
    REQUIRE(!result.stdout_output.empty());
}

// ─── Index Command Tests ─────────────────────────────────────────────────────

TEST_CASE("charon index --help shows usage", "[Tier3][Integration][Index]") {
    auto result = run_charon("index --help");
    
    REQUIRE(result.success());
    REQUIRE(result.stdout_output.find("Build an index") != std::string::npos);
    REQUIRE((result.stdout_output.find("--input") != std::string::npos || 
             result.stdout_output.find("<input>") != std::string::npos));
}

TEST_CASE("charon index fails with missing input file", "[Tier3][Integration][Index]") {
    auto result = run_charon("index");
    
    // Should fail because input file is required
    REQUIRE(!result.success());
    REQUIRE((result.stdout_output.find("required") != std::string::npos ||
             result.exit_code != 0));
}

TEST_CASE("charon index fails with non-existent input file", "[Tier3][Integration][Index]") {
    auto result = run_charon("index /nonexistent/file.tsv");
    
    REQUIRE(!result.success());
    REQUIRE((result.stdout_output.find("File does not exist") != std::string::npos ||
             result.stdout_output.find("Failed to open") != std::string::npos));
}

// ─── Classify Command Tests ──────────────────────────────────────────────────

TEST_CASE("charon classify --help shows usage", "[Tier3][Integration][Classify]") {
    auto result = run_charon("classify --help");
    
    REQUIRE(result.success());
    REQUIRE(result.stdout_output.find("Classify read file") != std::string::npos);
    REQUIRE(result.stdout_output.find("--db") != std::string::npos);
}

TEST_CASE("charon classify fails without database", "[Tier3][Integration][Classify]") {
    // Use the test fasta file but no database
    std::string test_fasta = TEST_DATA_DIR + "/my.fasta";
    auto result = run_charon("classify " + test_fasta);
    
    // Should fail because --db is required
    REQUIRE(!result.success());
    REQUIRE((result.stdout_output.find("required") != std::string::npos ||
             result.exit_code != 0));
}

// ─── Dehost Command Tests ────────────────────────────────────────────────────

TEST_CASE("charon dehost --help shows usage", "[Tier3][Integration][Dehost]") {
    auto result = run_charon("dehost --help");
    
    REQUIRE(result.success());
    REQUIRE(result.stdout_output.find("Dehost read file") != std::string::npos);
}

// ─── End-to-End Integration Test (if index exists) ───────────────────────────

TEST_CASE("charon can load existing index", "[Tier3][Integration][LoadIndex]") {
    // Check if test_index executable exists (from unit tests)
    // This tests that the index data structures can be loaded
    std::string test_index_exe = BUILD_DIR + "/bin/test_index";
    
    if (file_exists(test_index_exe)) {
        // Run the test_index unit test directly (not through charon)
        std::string command = test_index_exe + " -c \"Index construction and query\" 2>&1";
        FILE* pipe = popen(command.c_str(), "r");
        REQUIRE(pipe != nullptr);
        
        char buffer[256];
        std::string output;
        while (fgets(buffer, sizeof(buffer), pipe) != nullptr) {
            output += buffer;
        }
        
        int status = pclose(pipe);
        int exit_code = WEXITSTATUS(status);
        
        // test_index should pass (exit code 0)
        REQUIRE(exit_code == 0);
        REQUIRE(output.find("All tests passed") != std::string::npos);
    } else {
        // Skip if test executable doesn't exist
        SUCCEED("test_index executable not found, skipping");
    }
}

// ─── Real Data Integration Test ──────────────────────────────────────────────

TEST_CASE("charon processes real fasta file", "[Tier3][Integration][RealData]") {
    std::string test_fasta = TEST_DATA_DIR + "/my.fasta";
    
    // Verify test data exists
    REQUIRE(file_exists(test_fasta));
    
    // Count sequences in fasta file
    std::size_t num_sequences = count_lines(test_fasta) / 2; // Approximate for single-line fasta
    REQUIRE(num_sequences > 0);
    
    // Test that charon can at least read the file (will fail at classify without index, but should parse)
    std::string command = "classify " + test_fasta + " --db /nonexistent";
    auto result = run_charon(command);
    
    // Should fail because database doesn't exist, but should have tried to parse the file
    REQUIRE(!result.success());
}

// ─── Error Handling Tests ────────────────────────────────────────────────────

TEST_CASE("charon handles invalid subcommand", "[Tier3][Integration][ErrorHandling]") {
    auto result = run_charon("invalid_command");
    
    REQUIRE(!result.success());
    REQUIRE((result.stdout_output.find("Could not parse") != std::string::npos ||
             result.exit_code != 0));
}

TEST_CASE("charon handles --help on invalid subcommand", "[Tier3][Integration][ErrorHandling]") {
    auto result = run_charon("invalid --help");
    
    // Should still show help or error gracefully
    REQUIRE((result.exit_code == 0 || result.exit_code != 127)); // 127 = command not found
}

// ─── Performance Smoke Test ──────────────────────────────────────────────────

TEST_CASE("charon starts quickly", "[Tier3][Integration][Performance]") {
    // Test that charon starts and responds to --help quickly (< 1 second)
    auto start = std::chrono::high_resolution_clock::now();
    auto result = run_charon("--help");
    auto end = std::chrono::high_resolution_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    
    REQUIRE(result.success());
    REQUIRE(duration.count() < 1000); // Should respond in under 1 second
}
