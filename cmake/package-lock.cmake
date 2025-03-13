# SPDX-FileCopyrightText: 2006-2025, Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025, Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0

# CPM Package Lock
# This file should be committed to version control

# hibf
set (HIBF_VERSION af26c24dfbd8489760166acdb68b84921866ce3f)
CPMDeclarePackage (hibf
                   NAME hibf
                   GIT_TAG ${HIBF_VERSION} # main
                   GITHUB_REPOSITORY seqan/hibf
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_HIBF OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# sharg
set (SHARG_VERSION e9bc14ba8818f980727221dc936dbe6361eb87fd)
CPMDeclarePackage (sharg
                   NAME sharg
                   GIT_TAG ${SHARG_VERSION} # main
                   GITHUB_REPOSITORY seqan/sharg-parser
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SHARG OFF" "INSTALL_TDL OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING" "SHARG_NO_TDL ON"
)

# seqan3
set (SEQAN3_VERSION 7e0d88d15fc82b8b8a5548f6eebea8602faf6446)
CPMDeclarePackage (seqan3
                   NAME seqan3
                   GIT_TAG ${SEQAN3_VERSION} # main
                   GITHUB_REPOSITORY seqan/seqan3
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "INSTALL_SEQAN3 OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# libjst
set (LIBJST_VERSION bb3da768840026bf62a5d9869c04da6c8460a34f)
CPMDeclarePackage (libjst
                   NAME libjst
                   GIT_TAG ${LIBJST_VERSION}
                   GITHUB_REPOSITORY rrahn/libjst
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "LIBJST_DEVELOPER_MODE OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# libspm
set (LIBSPM_VERSION 094efd745853f617bc60456ca73262018d3e71d5)
CPMDeclarePackage (libspm
                   NAME libspm
                   GIT_TAG ${LIBSPM_VERSION}
                   GITHUB_REPOSITORY rrahn/libspm
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "LIBSPM_DEVELOPER_MODE OFF" "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# googletest
set (GOOGLETEST_VERSION 1.16.0)
CPMDeclarePackage (googletest
                   NAME GTest
                   VERSION ${GOOGLETEST_VERSION}
                   GITHUB_REPOSITORY google/googletest
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BUILD_GMOCK OFF" "INSTALL_GTEST OFF" "CMAKE_CXX_STANDARD 20"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

#googlebenchmark
set (GOOGLEBENCHMARK_VERSION 1.9.1)
CPMDeclarePackage (googlebenchmark
                   NAME benchmark
                   VERSION ${GOOGLEBENCHMARK_VERSION}
                   GITHUB_REPOSITORY google/benchmark
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "BENCHMARK_ENABLE_TESTING OFF" "BENCHMARK_ENABLE_WERROR OFF" "BENCHMARK_ENABLE_INSTALL OFF"
                           "CMAKE_MESSAGE_LOG_LEVEL WARNING" "CMAKE_CXX_STANDARD 20"
)

# catch2
set (CATCH2_VERSION 3.8.0)
CPMDeclarePackage (catch2
                   NAME Catch2
                   VERSION ${CATCH2_VERSION}
                   GITHUB_REPOSITORY catchorg/Catch2
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
                   OPTIONS "CMAKE_MESSAGE_LOG_LEVEL WARNING"
)

# use_ccache
set (USE_CCACHE_VERSION d2a54ef555b6fc2d496a4c9506dbeb7cf899ce37)
CPMDeclarePackage (use_ccache
                   NAME use_ccache
                   GIT_TAG ${USE_CCACHE_VERSION} # main
                   GITHUB_REPOSITORY seqan/cmake-scripts
                   SOURCE_SUBDIR ccache
                   SYSTEM TRUE
                   EXCLUDE_FROM_ALL TRUE
)
