# =============================================================================
# sifter
# =============================================================================
set(PROJECT_NAME_STR sifter)
cmake_minimum_required( VERSION 3.4 FATAL_ERROR )
project( ${PROJECT_NAME_STR} VERSION 0.0.1 LANGUAGES CXX )

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

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_link_libraries(${PROJECT_NAME} seqan3::seqan3 plog::plog)

# enable_testing()
# add_subdirectory(test)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
