# CPM Package Lock
# This file should be committed to version control

# seqan3
set (SEQAN3_VERSION 30bdf8d0a5d59b79342b472504b95ae50c33da6d)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${SEQAN3_VERSION}
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)
set (PLOG_VERSION e21baecd4753f14da64ede979c5a19302618b752)
CPMDeclarePackage (plog
                   NAME plog
                   GIT_TAG ${PLOG_VERSION}
                   GITHUB_REPOSITORY SergiusTheBest/plog
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_PLOG OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)
CPMDeclarePackage( statslib
        NAME statslib
        GITHUB_REPOSITORY kthohr/stats
        VERSION 3.4.0
        GIT_SHALLOW TRUE
        DOWNLOAD_ONLY TRUE
)
CPMDeclarePackage ( gcem
        NAME gcem
        GITHUB_REPOSITORY kthohr/gcem
        VERSION 1.18.0
        GIT_SHALLOW TRUE
        DOWNLOAD_ONLY TRUE
)
CPMDeclarePackage ( bzip2
        NAME bzip2
        GITHUB_REPOSITORY libarchive/bzip2
        GIT_TAG 6a8690f
        GIT_SHALLOW TRUE
        DOWNLOAD_ONLY TRUE
)
CPMDeclarePackage ( zlib
        NAME zlib
        GITHUB_REPOSITORY madler/zlib
        VERSION 1.3.1
        GIT_SHALLOW TRUE
        DOWNLOAD_ONLY TRUE
)
CPMDeclarePackage( gzip
        NAME gzip
        GITHUB_REPOSITORY mapbox/gzip-hpp
        GIT_TAG 7546b35
        GIT_SHALLOW TRUE
        DOWNLOAD_ONLY TRUE
)

# googletest
set (GOOGLETEST_VERSION 1.15.2)
CPMDeclarePackage (googletest
                   NAME googletest
                   VERSION ${GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
                           "CMAKE_CXX_STANDARD 20"
)
