# NB follows example detailed in https://www.marcusfolkesson.se/blog/git-version-in-cmake/

find_package(Git)

if(GIT_EXECUTABLE)
    get_filename_component(WORKING_DIR ${SRC} DIRECTORY)
    execute_process(
            COMMAND ${GIT_EXECUTABLE} describe --tags --dirty --always
            WORKING_DIRECTORY ${WORKING_DIR}
            OUTPUT_VARIABLE CHARON_VERSION
            RESULT_VARIABLE ERROR_CODE
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )
endif()

if(CHARON_VERSION STREQUAL "")
    set(CHARON_VERSION 0.0.0-unknown)
    message(WARNING "Failed to determine version from Git tags. Using default version \"${CHARON_VERSION}\".")
endif()

configure_file(${SRC} ${DST} @ONLY)