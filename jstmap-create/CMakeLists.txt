file(GLOB HEADER_LIST CONFIGURE_DEPENDS "jstmap/create/*.hpp")
### Base interface target for the build subcommand.
add_library (jstmap_create_base INTERFACE)
target_include_directories (jstmap_create_base INTERFACE ../jstmap-create
                                                         libjst::libjst
                                                         seqan3::seqan3
                                                         jstmap::global
                                                         ${SEQAN_INCLUDE_DIRS})
target_compile_features (jstmap_create_base INTERFACE cxx_std_20)
target_link_libraries (jstmap_create_base INTERFACE libjst::libjst seqan3::seqan3 seqan::seqan2 jstmap::global)
add_library (jstmap::create::base ALIAS jstmap_create_base)

### Create object library for better build times
add_library(jstmap_create_input OBJECT jstmap/create/load_sequence.cpp
                                       jstmap/create/load_sequence.hpp)
target_link_libraries (jstmap_create_input PUBLIC jstmap::create::base)

### Building the journaled sequence tree
add_library(jstmap_create_builder OBJECT jstmap/create/journaled_sequence_tree_builder.cpp
                                         jstmap/create/journaled_sequence_tree_builder.hpp
                                         jstmap/create/serialise_jst.hpp)
target_link_libraries (jstmap_create_builder PUBLIC jstmap::create::base)

### Building the journaled sequence tree
add_library(jstmap_serialise_jst OBJECT jstmap/create/serialise_jst.hpp)
target_link_libraries (jstmap_serialise_jst PUBLIC jstmap::create::base libjst::libjst)

### Parsing the vcf file
add_library(jstmap_create_vcf_parser OBJECT jstmap/create/vcf_parser2.cpp
                                            jstmap/create/stripped_vcf_record.cpp
                                            jstmap/create/vcf_parser.hpp
                                            jstmap/create/stripped_vcf_record.hpp)
target_link_libraries (jstmap_create_vcf_parser PUBLIC jstmap::create::base)
### Create static library for build subcommand
add_library (jstmap_create STATIC jstmap/create/create_main.cpp
                                  jstmap/create/create_main.hpp)
target_link_libraries (jstmap_create PUBLIC
                                     jstmap_create_input
                                     jstmap_create_builder
                                     jstmap_serialise_jst
                                     jstmap_create_vcf_parser)
add_library (jstmap::create ALIAS jstmap_create)

unset (HEADER_LIST)
