
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/search/*.hpp")

### Base interface target for the index subcommand.
add_library (jstmap_search_base INTERFACE)
target_include_directories (jstmap_search_base INTERFACE ../jstmap-search seqan3::seqan3)
target_link_libraries (jstmap_search_base INTERFACE seqan3::seqan3)
add_library (jstmap::search::base ALIAS jstmap_search_base)

### Create object library for better build times
add_library(jstmap_search_input_queries OBJECT jstmap/search/load_queries.cpp
                                               jstmap/search/load_queries.hpp
                                               jstmap/search/global_types.hpp)
target_include_directories (jstmap_search_input_queries PUBLIC jstmap::search::base)
target_link_libraries (jstmap_search_input_queries PUBLIC jstmap::search::base)

add_library(jstmap_search_input_jst OBJECT jstmap/search/load_jst.cpp
                                           jstmap/search/load_jst.hpp
                                           jstmap/search/global_types.hpp)
target_include_directories (jstmap_search_input_jst PUBLIC jstmap::search::base libjst::libjst)
target_link_libraries (jstmap_search_input_jst PUBLIC jstmap::search::base libjst::libjst)

### Building the journaled sequence tree
add_library(jstmap_search_search OBJECT jstmap/search/search_queries.cpp
                                        jstmap/search/search_queries.hpp
                                        jstmap/search/global_types.hpp)
target_include_directories (jstmap_search_search PUBLIC jstmap::search::base libjst::libjst)
target_link_libraries (jstmap_search_search PUBLIC jstmap::search::base libjst::libjst)

### Building the journaled sequence tree
add_library(jstmap_search_output OBJECT jstmap/search/write_results.cpp
                                        jstmap/search/write_results.hpp
                                        jstmap/search/global_types.hpp)
target_include_directories (jstmap_search_output PUBLIC jstmap::search::base libjst::libjst)
target_link_libraries (jstmap_search_output PUBLIC jstmap::search::base libjst::libjst)

### Create static library for index subcommand
add_library (jstmap_search STATIC jstmap/search/search_main.cpp)
target_include_directories (jstmap_search PUBLIC jstmap_search_input_jst jstmap_search_search jstmap_search_input_queries jstmap_search_output)
target_link_libraries (jstmap_search PUBLIC jstmap_search_input_jst jstmap_search_search jstmap_search_input_queries jstmap_search_output)
add_library (jstmap::search ALIAS jstmap_search)

unset (HEADER_LIST)
