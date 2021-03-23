
file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/simulate/*.hpp")

### Base interface target for the simulate subcommand.
add_library (jstmap_simulate_base INTERFACE)
target_include_directories (jstmap_simulate_base INTERFACE ../jstmap-simulate seqan3::seqan3)
target_compile_features (jstmap_simulate_base INTERFACE cxx_std_17)
target_link_libraries (jstmap_simulate_base INTERFACE seqan3::seqan3)
add_library (jstmap::simulate::base ALIAS jstmap_simulate_base)

### Create object library for better build times
add_library(jstmap_simulate_input OBJECT jstmap/simulate/global_types.hpp)
#add_library(jstmap_simulate_input OBJECT jstmap/simulate/load_sequence.cpp
#                                      jstmap/simulate/load_sequence.hpp
#                                      jstmap/simulate/global_types.hpp)
target_include_directories (jstmap_simulate_input PUBLIC jstmap::simulate::base)
target_link_libraries (jstmap_simulate_input PUBLIC jstmap::simulate::base)

### Building the journaled sequence tree
add_library(jstmap_simulate_build OBJECT jstmap/simulate/global_types.hpp)
# add_library(jstmap_simulate_build OBJECT jstmap/simulate/journaled_sequence_tree_builder.cpp
#                                       jstmap/simulate/journaled_sequence_tree_builder.hpp
#                                       jstmap/simulate/global_types.hpp)
target_include_directories (jstmap_simulate_build PUBLIC jstmap::simulate::base libjst::libjst)
target_link_libraries (jstmap_simulate_build PUBLIC jstmap::simulate::base libjst::libjst)

### Create static library for simulate subcommand
add_library (jstmap_simulate STATIC jstmap/simulate/simulate_main.cpp
                                 jstmap/simulate/simulate_main.hpp)
# target_include_directories (jstmap_simulate PUBLIC jstmap::simulate::base libjst::libjst)
# target_link_libraries (jstmap_simulate PUBLIC jstmap::simulate::base libjst::libjst)
target_include_directories (jstmap_simulate PUBLIC jstmap_simulate_input jstmap_simulate_build)
target_link_libraries (jstmap_simulate PUBLIC jstmap_simulate_input jstmap_simulate_build)
add_library (jstmap::simulate ALIAS jstmap_simulate)

unset (HEADER_LIST)