
### Base interface target for the global jstmap library.
add_library (jstmap_index_base INTERFACE)
target_include_directories (jstmap_index_base INTERFACE ../jstmap-index)
target_compile_features (jstmap_index_base INTERFACE cxx_std_20)
target_link_libraries (jstmap_index_base INTERFACE libjst::libjst seqan3::seqan3 seqan::seqan2 jstmap::global)

### Create object library for the index creation.
add_library(jstmap_index_create OBJECT jstmap/index/create_index.cpp jstmap/index/create_index.hpp)
target_link_libraries (jstmap_index_create PUBLIC jstmap_index_base)

### Create object library for saveing the index.
add_library(jstmap_index_save OBJECT jstmap/index/save_index.cpp jstmap/index/save_index.hpp)
target_link_libraries (jstmap_index_save PUBLIC jstmap_index_base)

### Create static library for the index main.
add_library(jstmap_index_main STATIC jstmap/index/index_main.cpp jstmap/index/index_main.hpp)
target_link_libraries (jstmap_index_main PUBLIC jstmap_index_create jstmap_index_save)
add_library (jstmap::index ALIAS jstmap_index_main)
