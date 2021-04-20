cmake_minimum_required (VERSION 3.14)

set (LIBJST_TEST_DEPENDENCIES "jstmap::simulate" "libjst::libjst")

macro (add_libjst_test test_filename)
    add_api_test(${test_filename} ${LIBJST_TEST_DEPENDENCIES})
endmacro()

add_libjst_test (journaled_sequence_tree_test.cpp)
add_libjst_test (journaled_sequence_tree_traversal_test.cpp)
add_libjst_test (journal_sequence_tree_fuzzy_test.cpp)

add_libjst_test (journal_sequence_tree_search_test.cpp)
target_use_datasources (journal_sequence_tree_search_test FILES sim_refx5.jst)

add_libjst_test (journal_decorator_test.cpp)
add_libjst_test (journal_decorator_iterator_test.cpp)

add_subdirectories ()

unset (LIBJST_TEST_DEPENDENCIES)
