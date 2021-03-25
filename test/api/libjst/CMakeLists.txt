cmake_minimum_required (VERSION 3.8)

macro (add_libjst_test test_filename)
    add_api_test(${test_filename} "libjst::libjst")
endmacro()

add_libjst_test (journaled_sequence_tree_test.cpp)
add_libjst_test (journaled_sequence_tree_traversal_test.cpp)
add_libjst_test (journal_decorator_test.cpp)
add_libjst_test (journal_decorator_iterator_test.cpp)

add_subdirectory (detail)
