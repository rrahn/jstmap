// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>

#include <jstmap/simulate/load_reference.hpp>
#include <jstmap/simulate/simulate_alignment.hpp>

TEST(jstmap_simulate, simulate_alignment)
{
    std::filesystem::path sequence_file{DATADIR"sim_reads_ref1x10.fa"};
    jstmap::aligned_sequence_t reference = jstmap::load_reference(sequence_file);
    jstmap::alignment_t alignment = jstmap::simulate_alignment(reference, 0.39); // ceil(20*0.39) = 8 => 4 SNPs, 2 Inserts, 2 Deletes
    EXPECT_EQ(alignment.first.size(), 22u);
    EXPECT_EQ(alignment.second.size(), 22u);
    size_t variants = 0;
    for(size_t i = 0; i < alignment.first.size(); ++i){
      variants += (alignment.first[i] != alignment.second[i]);
    }
    EXPECT_EQ(variants, 8u);
    alignment.first.erase(std::remove(alignment.first.begin(), alignment.first.end(), seqan3::gapped<seqan3::dna5>{}.assign_char('-')), alignment.first.end());
    jstmap::aligned_sequence_t stripped_reference = alignment.first;
    EXPECT_EQ(reference, stripped_reference);
}

TEST(jstmap_simulate, simulate_alignment_error_rate_zero)
{
    std::filesystem::path sequence_file{DATADIR"sim_reads_ref1x10.fa"};
    jstmap::aligned_sequence_t reference = jstmap::load_reference(sequence_file);
    jstmap::alignment_t alignment = jstmap::simulate_alignment(reference, 0); // 20*0.4 = 8 => 4 SNPs, 2 Inserts, 2 Deletes
    EXPECT_EQ(alignment.first.size(), 20u);
    EXPECT_EQ(alignment.second.size(), 20u);
    size_t variants = 0;
    for(size_t i = 0; i < alignment.first.size(); ++i){
      variants += (alignment.first[i] != alignment.second[i]);
    }
    EXPECT_EQ(variants, 0u);
}

TEST(jstmap_simulate, simulate_alignment_negative_error_rate)
{
    std::filesystem::path sequence_file{DATADIR"sim_reads_ref1x10.fa"};
    jstmap::aligned_sequence_t reference = jstmap::load_reference(sequence_file);
    EXPECT_THROW(jstmap::alignment_t alignment = jstmap::simulate_alignment(reference, -1), std::runtime_error);
}

TEST(jstmap_simulate, simulate_alignment_too_high_error_rate)
{
    std::filesystem::path sequence_file{DATADIR"sim_reads_ref1x10.fa"};
    jstmap::aligned_sequence_t reference = jstmap::load_reference(sequence_file);
    EXPECT_THROW(jstmap::alignment_t alignment = jstmap::simulate_alignment(reference, 2), std::runtime_error);
}
