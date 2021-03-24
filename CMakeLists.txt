cmake_minimum_required (VERSION 3.12)

## CUSTOMISE

# Define the application name and version.
project (jstmap VERSION 0.0.1 LANGUAGES CXX)
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 REQUIRED HINTS lib/seqan3/build_system)

# Dependency: SeqAn.
set (SEQAN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/lib/seqan/include")
find_package (SeqAn REQUIRED HINTS "lib/seqan/util/cmake" NO_DEFAULT_PATH)

## APPLICATION

add_subdirectory (libjst)
add_subdirectory (jstmap-index)
add_subdirectory (jstmap-search)
add_subdirectory (jstmap-simulate)
add_subdirectory (jstmap)

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

add_subdirectory (test/api EXCLUDE_FROM_ALL)
add_subdirectory (test/cli EXCLUDE_FROM_ALL)
