# CPM Package Lock
# This file should be committed to version control

# seqan3
set(SEQAN3_VERSION 30bdf8d0a5d59b79342b472504b95ae50c33da6d)
CPMDeclarePackage(seqan3
    NAME seqan3
    GIT_TAG ${SEQAN3_VERSION}
    GITHUB_REPOSITORY seqan/seqan3
    SYSTEM TRUE
    EXCLUDE_FROM_ALL TRUE
    OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# plog
set(PLOG_VERSION e21baecd4753f14da64ede979c5a19302618b752)
CPMDeclarePackage(plog
    NAME plog
    GIT_TAG ${PLOG_VERSION}
    GITHUB_REPOSITORY SergiusTheBest/plog
    SYSTEM TRUE
    EXCLUDE_FROM_ALL TRUE
    OPTIONS "INSTALL_PLOG OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# statslib
CPMDeclarePackage(statslib
    NAME statslib
    GITHUB_REPOSITORY kthohr/stats
    VERSION 3.4.0
    GIT_SHALLOW TRUE
    DOWNLOAD_ONLY TRUE
)

# gcem
CPMDeclarePackage(gcem
    NAME gcem
    GITHUB_REPOSITORY kthohr/gcem
    VERSION 1.18.0
    GIT_SHALLOW TRUE
    DOWNLOAD_ONLY TRUE
)

# bzip2
CPMDeclarePackage(bzip2
    NAME bzip2
    GITHUB_REPOSITORY libarchive/bzip2
    GIT_TAG 6a8690f
    GIT_SHALLOW TRUE
    DOWNLOAD_ONLY TRUE
)

# zlib
CPMDeclarePackage(zlib
    NAME zlib
    GITHUB_REPOSITORY madler/zlib
    VERSION 1.3.1
    GIT_SHALLOW TRUE
    DOWNLOAD_ONLY TRUE
)

# ankerl::unordered_dense - high-performance hash maps
CPMDeclarePackage(unordered_dense
    NAME unordered_dense
    GITHUB_REPOSITORY martinus/unordered_dense
    VERSION 4.4.0
    GIT_SHALLOW TRUE
    DOWNLOAD_ONLY TRUE
)

# gzip-hpp
CPMDeclarePackage(gzip
    NAME gzip
    GITHUB_REPOSITORY mapbox/gzip-hpp
    GIT_TAG 7546b35
    GIT_SHALLOW TRUE
    DOWNLOAD_ONLY TRUE
)

# Catch2
CPMDeclarePackage(Catch2
    NAME Catch2
    VERSION 3.6.0
    GITHUB_REPOSITORY catchorg/Catch2
    SYSTEM TRUE
    EXCLUDE_FROM_ALL TRUE
    OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)
