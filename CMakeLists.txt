# =============================================================================
# charon
# =============================================================================
cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(charon LANGUAGES CXX)

# -----------------------------------------------------------------------------
# C++ standard
# -----------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# output directories
# -----------------------------------------------------------------------------
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_BINDIR "bin")

# -----------------------------------------------------------------------------
# version header
# -----------------------------------------------------------------------------
add_custom_target(version
    ${CMAKE_COMMAND}
    -D SRC=${CMAKE_SOURCE_DIR}/version.h.in
    -D DST=${CMAKE_SOURCE_DIR}/include/version.h
    -P ${CMAKE_SOURCE_DIR}/cmake/GenerateVersion.cmake
)

# -----------------------------------------------------------------------------
# dependencies
# -----------------------------------------------------------------------------
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)

set(CPM_INDENT "CMake Package Manager CPM: ")
include(${PROJECT_SOURCE_DIR}/cmake/CPM.cmake)
CPMUsePackageLock(${PROJECT_SOURCE_DIR}/cmake/package-lock.cmake)

CPMGetPackage(seqan3)
CPMGetPackage(plog)
CPMGetPackage(gcem)
CPMGetPackage(statslib)
CPMGetPackage(bzip2)
CPMGetPackage(zlib)
CPMGetPackage(gzip)
CPMGetPackage(Catch2)
CPMGetPackage(unordered_dense)

# -----------------------------------------------------------------------------
# main executable
# -----------------------------------------------------------------------------
# List source files explicitly — do not use GLOB_RECURSE
# (CMake will not detect new files without re-running cmake if GLOB is used)
set(CHARON_SOURCES
    src/main.cpp
    src/classify_main.cpp
    src/dehost_main.cpp
    src/index_main.cpp
    src/load_index.cpp
    src/main.cpp
    src/utils.cpp
    # add new .cpp files here explicitly
)

add_executable(charon ${CHARON_SOURCES})

add_dependencies(charon version)

target_include_directories(charon
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/lib/include
        ${gcem_SOURCE_DIR}/include
        ${statslib_SOURCE_DIR}/include
        ${bzip2_SOURCE_DIR}
        ${zlib_SOURCE_DIR}
        ${gzip_SOURCE_DIR}/include
        ${unordered_dense_SOURCE_DIR}/include
)

target_compile_definitions(charon PRIVATE ZLIB_CONST)

target_link_libraries(charon
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
)

install(TARGETS charon RUNTIME DESTINATION bin)

# -----------------------------------------------------------------------------
# tests
# -----------------------------------------------------------------------------
enable_testing()

# Utility functions under test — list only the src files each test needs,
# not the full executable sources, to keep compilation fast and dependencies minimal
add_executable(test_utils
    test/test_utils.cpp
    src/utils.cpp
)

target_include_directories(test_utils
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/lib/include
        ${gcem_SOURCE_DIR}/include
        ${statslib_SOURCE_DIR}/include
        ${bzip2_SOURCE_DIR}
        ${zlib_SOURCE_DIR}
        ${gzip_SOURCE_DIR}/include
        ${unordered_dense_SOURCE_DIR}/include
)

target_compile_definitions(test_utils PRIVATE ZLIB_CONST)

target_link_libraries(test_utils
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
        Catch2::Catch2WithMain
)

add_test(NAME test_utils COMMAND test_utils)

# Test InputSummary
add_executable(test_input_summary
    test/test_input_summary.cpp
)
target_include_directories(test_input_summary
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_link_libraries(test_input_summary
    PRIVATE
        seqan3::seqan3
        plog::plog
        Catch2::Catch2WithMain
)
add_test(NAME test_input_summary COMMAND test_input_summary)

add_executable(test_input_stats
    test/test_input_stats.cpp
)
target_include_directories(test_input_stats
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_link_libraries(test_input_stats
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
        Catch2::Catch2WithMain
)
add_test(NAME test_input_stats COMMAND test_input_stats)

# Test Index
add_executable(test_index
    test/test_index.cpp
)
target_include_directories(test_index
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/lib/include
        ${gcem_SOURCE_DIR}/include
        ${statslib_SOURCE_DIR}/include
        ${bzip2_SOURCE_DIR}
        ${zlib_SOURCE_DIR}
        ${gzip_SOURCE_DIR}/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_compile_definitions(test_index PRIVATE ZLIB_CONST)
target_link_libraries(test_index
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
        Catch2::Catch2WithMain
)
add_test(NAME test_index COMMAND test_index)

# Test classify_stats (Tier 1 - pure functions)
add_executable(test_classify_stats
    test/test_classify_stats.cpp
    test/model_defaults_stub.cpp
)
target_include_directories(test_classify_stats
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/include
        ${gcem_SOURCE_DIR}/include
        ${statslib_SOURCE_DIR}/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_link_libraries(test_classify_stats
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
        Catch2::Catch2WithMain
)
add_test(NAME test_classify_stats COMMAND test_classify_stats)

# Test utils I/O functions (Tier 1 - file operations)
add_executable(test_utils_io
    test/test_utils_io.cpp
    src/utils.cpp
)
target_include_directories(test_utils_io
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/lib/include
        ${gcem_SOURCE_DIR}/include
        ${statslib_SOURCE_DIR}/include
        ${bzip2_SOURCE_DIR}
        ${zlib_SOURCE_DIR}
        ${gzip_SOURCE_DIR}/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_compile_definitions(test_utils_io PRIVATE ZLIB_CONST)
target_link_libraries(test_utils_io
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
        Catch2::Catch2WithMain
)
add_test(NAME test_utils_io COMMAND test_utils_io)

# Test InputSummary extended (Tier 1 - data structure)
add_executable(test_input_summary_extended
    test/test_input_summary_extended.cpp
)
target_include_directories(test_input_summary_extended
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_link_libraries(test_input_summary_extended
    PRIVATE
        seqan3::seqan3
        plog::plog
        Catch2::Catch2WithMain
)
add_test(NAME test_input_summary_extended COMMAND test_input_summary_extended)

# Test InputStats extended (Tier 1 - data structure)
add_executable(test_input_stats_extended
    test/test_input_stats_extended.cpp
)
target_include_directories(test_input_stats_extended
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_link_libraries(test_input_stats_extended
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
        Catch2::Catch2WithMain
)
add_test(NAME test_input_stats_extended COMMAND test_input_stats_extended)

# Note: test_read_entry moved to test/unfinished/ - seqan3 API compatibility issues
# Test ReadEntry (Tier 2 - classification logic)
# add_executable(test_read_entry
#     test/test_read_entry.cpp
# )
# target_include_directories(test_read_entry
#     PRIVATE
#         ${PROJECT_SOURCE_DIR}/include
#         ${PROJECT_SOURCE_DIR}/lib/include
#         ${gcem_SOURCE_DIR}/include
#         ${statslib_SOURCE_DIR}/include
#         ${unordered_dense_SOURCE_DIR}/include
# )
# target_link_libraries(test_read_entry
#     PRIVATE
#         seqan3::seqan3
#         plog::plog
#         OpenMP::OpenMP_CXX
#         Catch2::Catch2WithMain
# )
# add_test(NAME test_read_entry COMMAND test_read_entry)

# Note: test_result moved to test/unfinished/ - depends on test_read_entry fixes
# Test Result (Tier 2 - template class with I/O)
# add_executable(test_result
#     test/test_result.cpp
# )
# target_include_directories(test_result
#     PRIVATE
#         ${PROJECT_SOURCE_DIR}/include
#         ${PROJECT_SOURCE_DIR}/lib/include
#         ${gcem_SOURCE_DIR}/include
#         ${statslib_SOURCE_DIR}/include
#         ${unordered_dense_SOURCE_DIR}/include
# )
# target_link_libraries(test_result
#     PRIVATE
#         seqan3::seqan3
#         plog::plog
#         OpenMP::OpenMP_CXX
#         Catch2::Catch2WithMain
# )
# add_test(NAME test_result COMMAND test_result)

# Test load_index (critical I/O function)
add_executable(test_load_index
    test/test_load_index.cpp
    src/load_index.cpp
)
target_include_directories(test_load_index
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/include
        ${unordered_dense_SOURCE_DIR}/include
)
target_link_libraries(test_load_index
    PRIVATE
        seqan3::seqan3
        plog::plog
        OpenMP::OpenMP_CXX
        Catch2::Catch2WithMain
)
add_test(NAME test_load_index COMMAND test_load_index)

# Test Integration Tests (end-to-end CLI testing)
add_executable(test_integration
    test/test_integration.cpp
)
target_include_directories(test_integration
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/lib/include
        ${plog_SOURCE_DIR}/include
)
target_link_libraries(test_integration
    PRIVATE
        plog::plog
        Catch2::Catch2WithMain
)
add_test(NAME test_integration COMMAND test_integration)