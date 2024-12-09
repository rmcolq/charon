# =============================================================================
# sifter
# =============================================================================
set(PROJECT_NAME_STR sifter)
cmake_minimum_required( VERSION 3.9 FATAL_ERROR )
project( ${PROJECT_NAME_STR} VERSION 0.0.1 LANGUAGES CXX )

# Use C++20
set(CMAKE_CXX_STANDARD 20)
# Require (at least) it
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Don't use e.g. GNU extension (like -std=gnu++11) for portability
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# dependencies and 3rd party libraries
# -----------------------------------------------------------------------------
set (PROGRAM_SUBMODULES_DIR
     "${CMAKE_CURRENT_LIST_DIR}/lib"
     CACHE STRING "Directory containing submodules."
)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set (CMAKE_INSTALL_BINDIR "bin")

set (CPM_INDENT "CMake Package Manager CPM: ")
include (${PROJECT_SOURCE_DIR}/cmake/CPM.cmake)
CPMUsePackageLock (${PROJECT_SOURCE_DIR}/cmake/package-lock.cmake)

CPMGetPackage (seqan3)
CPMGetPackage (plog)

# Require OPENMP
find_package(OpenMP REQUIRED)

# -----------------------------------------------------------------------------
# install
# -----------------------------------------------------------------------------

include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_SOURCE_DIR}/lib/include
)

file(GLOB_RECURSE SRC_FILES
        ${PROJECT_SOURCE_DIR}/src/*.cpp
)

add_executable(${PROJECT_NAME} ${SRC_FILES}
        src/main.cpp
)

target_link_libraries(${PROJECT_NAME} seqan3::seqan3 plog::plog OpenMP::OpenMP_CXX)

include_directories(
        ${PROJECT_SOURCE_DIR}/lib/include
)
file(GLOB_RECURSE TEST_SRC_FILES
        ${PROJECT_SOURCE_DIR}/scratch/my_app.cpp
)
add_executable(test_sifter ${TEST_SRC_FILES}
        scratch/my_app.cpp
)
target_link_libraries(test_sifter seqan3::seqan3)

#target_link_libraries(${PROJECT_NAME} seqan3::seqan3 plog::plog)

# enable_testing()
# add_subdirectory(test)

install(TARGETS ${PROJECT_NAME} test_sifter RUNTIME DESTINATION bin)
