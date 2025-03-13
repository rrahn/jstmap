// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <jstmap/global/load_jst.hpp>

TEST(jstmap_index, load_jst)
{
    using spm::operator""_dna4;

    std::filesystem::path jst_file{DATADIR"sim_refx5_p0.jst"};
    auto jst = jstmap::load_jst(jst_file);

    EXPECT_EQ(jst.size(), 5u);
    EXPECT_RANGE_EQ(jst.sequence_at(0),
                    "TATGCACCAGAGTATGGAAGCATAAGCTCTGCATGCAAAGGTACATCAGATCCTGCGGTTGGGTGCCAACCCAAGTGTGTTCACGGGCGC"_dna4);
    EXPECT_RANGE_EQ(jst.sequence_at(1),
                    "TTGACAGACATCGGAGGATGGTGCACACTCACTCGACCAGCGCAAAGCACAGGATCTCACGGGCGGACATCTCTTAGGTCAGTCATCGTGGAGGAATGCT"_dna4);
    EXPECT_RANGE_EQ(jst.sequence_at(2),
                    "TGTACGTTCTTTTGGCTTCCCCTAACACGGCGGGCGTCTCCGGTACGTATCCTGTCGGTACACCCCTTAAGCCCCTAGGCCCGAAGAACATAGCGCATTTCACGCTCTCT"_dna4);
    EXPECT_RANGE_EQ(jst.sequence_at(3),
                    "ACGAATGACCGCAACGATCAAATGGGCGAGAACAACTAATTCCGATTCATGGGGTTTGTGGATTGTGACACAGCGCGCCCGCTAC"_dna4);
    EXPECT_RANGE_EQ(jst.sequence_at(4),
                    "TGCGGGACGTGAGGACGCCCAATTCTGCCAAGGATTATTTAGGGTGTTTCACTAGAGTTATGCGCCGACCCCGGTTGGACCAGCTTGCATTCGAAACTGCGTTA"_dna4);
}

TEST(jstmap_index, load_jst_empty_path)
{
    EXPECT_THROW(jstmap::load_jst(std::filesystem::path{}), std::runtime_error);
}

TEST(jstmap_index, load_jst_unkown_path)
{
    std::filesystem::path unknown_jst_file{DATADIR"unknown.jst"};
    EXPECT_THROW(jstmap::load_jst(unknown_jst_file), std::runtime_error);
}
