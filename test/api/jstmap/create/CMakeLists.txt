cmake_minimum_required (VERSION 3.14)

macro (add_jstmap_create_test test_filename)
    add_api_test(${test_filename} "jstmap::create")
endmacro ()

# add_jstmap_create_test (build_journaled_sequence_tree_test.cpp)

# add_jstmap_create_test (load_sequence_test.cpp)
# target_use_datasources (load_sequence_test FILES in.fasta)

# add_jstmap_create_test (vcf_parser_test.cpp)
# target_use_datasources (vcf_parser_test FILES in.fasta
#                                               empty.fa
#                                               sim_ref_10Kb.fasta.gz
#                                               sim_ref_10Kb_SNPs.vcf
#                                               sim_ref_10Kb_SNPs_haplotypes.fasta.gz
#                                               sim_ref_10Kb_SNP_INDELs.vcf
#                                               sim_ref_10Kb_SNP_INDELs_haplotypes.fasta.gz
#                                               sim_ref_10Kb_no_variants.vcf)
