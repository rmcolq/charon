# =============================================================================
# charon
# =============================================================================
set(PROJECT_NAME_STR charon)
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

add_custom_target(version
        ${CMAKE_COMMAND} -D SRC=${CMAKE_SOURCE_DIR}/version.h.in
        -D DST=${CMAKE_SOURCE_DIR}/include/version.h
        -P ${CMAKE_SOURCE_DIR}/cmake/GenerateVersion.cmake
)

project(${PROJECT_NAME_STR} LANGUAGES CXX)

# Use C++20
set(CMAKE_CXX_STANDARD 20)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# dependencies and 3rd party libraries
# -----------------------------------------------------------------------------
set(PROGRAM_SUBMODULES_DIR
        "${CMAKE_CURRENT_LIST_DIR}/lib"
        CACHE STRING "Directory containing submodules."
)

# Specify the directories where to store the built archives, libraries and executables
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_BINDIR "bin")

# Require OPENMP
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

add_definitions(-D ZLIB_CONST)

# -----------------------------------------------------------------------------
# install
# -----------------------------------------------------------------------------

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/lib/include
        ${gcem_SOURCE_DIR}/include
        ${statslib_SOURCE_DIR}/include
        ${bzip2_SOURCE_DIR}
        ${zlib_SOURCE_DIR}
        ${gzip_SOURCE_DIR}/include
)

file(GLOB_RECURSE SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
)

add_executable(${PROJECT_NAME} ${SRC_FILES}
        src/main.cpp
)

add_dependencies(${PROJECT_NAME} version)

target_link_libraries(${PROJECT_NAME} PRIVATE seqan3::seqan3 plog::plog OpenMP::OpenMP_CXX)

install(TARGETS ${PROJECT_NAME} charon RUNTIME DESTINATION bin)
