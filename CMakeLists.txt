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