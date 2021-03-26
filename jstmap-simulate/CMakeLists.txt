
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/simulate/*.hpp")

### Base interface target for the simulate subcommand.
add_library (jstmap_simulate_base INTERFACE)
target_include_directories (jstmap_simulate_base INTERFACE ../jstmap-simulate seqan3::seqan3)
target_compile_features (jstmap_simulate_base INTERFACE cxx_std_17)
target_link_libraries (jstmap_simulate_base INTERFACE seqan3::seqan3)
add_library (jstmap::simulate::base ALIAS jstmap_simulate_base)

### Create object library for better build times
add_library(jstmap_simulate_input OBJECT jstmap/simulate/load_reference.cpp
                                         jstmap/simulate/load_reference.hpp
                                         jstmap/simulate/global_types.hpp)

target_include_directories (jstmap_simulate_input PUBLIC jstmap::simulate::base)
target_link_libraries (jstmap_simulate_input PUBLIC jstmap::simulate::base)

### Building the journaled sequence tree
add_library(jstmap_simulate_build OBJECT jstmap/simulate/simulate_alignment.cpp
                                         jstmap/simulate/simulate_alignment.hpp
                                         jstmap/simulate/global_types.hpp)

target_include_directories (jstmap_simulate_build PUBLIC jstmap::simulate::base)
target_link_libraries (jstmap_simulate_build PUBLIC jstmap::simulate::base)

### Create static library for simulate subcommand
add_library (jstmap_simulate STATIC jstmap/simulate/simulate_main.cpp
                                 jstmap/simulate/simulate_main.hpp)
target_include_directories (jstmap_simulate PUBLIC jstmap_simulate_input jstmap_simulate_build jstmap_index_build libjst::libjst)
target_link_libraries (jstmap_simulate PUBLIC jstmap_simulate_input jstmap_simulate_build jstmap_index_build libjst::libjst)
add_library (jstmap::simulate ALIAS jstmap_simulate)

unset (HEADER_LIST)
