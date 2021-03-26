cmake_minimum_required (VERSION 3.8)

macro (add_jstmap_simulate_test test_filename)
    add_api_test(${test_filename} "jstmap::simulate")
endmacro ()

add_jstmap_simulate_test (load_reference_test.cpp)
target_use_datasources (load_reference_test FILES sim_reads_ref1x10.fa)

add_jstmap_simulate_test (simulate_alignment_test.cpp)
target_use_datasources (simulate_alignment_test FILES sim_reads_ref1x10.fa)
