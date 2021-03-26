
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/index/*.hpp")

### Base interface target for the index subcommand.
add_library (jstmap_index_base INTERFACE)
target_include_directories (jstmap_index_base INTERFACE ../jstmap-index seqan3::seqan3)
target_compile_features (jstmap_index_base INTERFACE cxx_std_17)
target_link_libraries (jstmap_index_base INTERFACE seqan3::seqan3)
add_library (jstmap::index::base ALIAS jstmap_index_base)

### Create object library for better build times
add_library(jstmap_index_input OBJECT jstmap/index/load_sequence.cpp
                                      jstmap/index/load_sequence.hpp
                                      jstmap/index/global_types.hpp)
target_include_directories (jstmap_index_input PUBLIC jstmap::index::base)
target_link_libraries (jstmap_index_input PUBLIC jstmap::index::base)

### Building the journaled sequence tree
add_library(jstmap_index_build OBJECT jstmap/index/journaled_sequence_tree_builder.cpp
                                      jstmap/index/journaled_sequence_tree_builder.hpp
                                      jstmap/index/serialise_jst.hpp
                                      jstmap/index/global_types.hpp)
target_include_directories (jstmap_index_build PUBLIC jstmap::index::base libjst::libjst)
target_link_libraries (jstmap_index_build PUBLIC jstmap::index::base libjst::libjst)

### Building the journaled sequence tree
add_library(jstmap_serialise_jst OBJECT jstmap/index/serialise_jst.hpp
                                      jstmap/index/global_types.hpp)
target_include_directories (jstmap_serialise_jst PUBLIC jstmap::index::base libjst::libjst)
target_link_libraries (jstmap_serialise_jst PUBLIC jstmap::index::base libjst::libjst)

### Create static library for index subcommand
add_library (jstmap_index STATIC jstmap/index/index_main.cpp
                                 jstmap/index/index_main.hpp)
target_include_directories (jstmap_index PUBLIC jstmap_index_input jstmap_index_build)
target_link_libraries (jstmap_index PUBLIC jstmap_index_input jstmap_index_build)
add_library (jstmap::index ALIAS jstmap_index)

unset (HEADER_LIST)
