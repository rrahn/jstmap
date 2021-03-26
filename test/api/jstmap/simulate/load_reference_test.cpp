// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <filesystem>

#include <jstmap/simulate/load_reference.hpp>

TEST(jstmap_index, load_reference)
{
    std::filesystem::path sequence_file{DATADIR"sim_reads_ref1x10.fa"};
    jstmap::aligned_sequence_t reference = jstmap::load_reference(sequence_file);
    EXPECT_EQ(reference.size(), 20u);
}

TEST(jstmap_index, load_reference_empty_path)
{
    EXPECT_THROW(jstmap::load_reference(std::filesystem::path{}), std::runtime_error);
}

TEST(jstmap_index, load_reference_unkown_path)
{
    std::filesystem::path unknown_sequence_file{DATADIR"unknown.fa"};
    EXPECT_THROW(jstmap::load_reference(unknown_sequence_file), std::runtime_error);
}

TEST(jstmap_index, load_reference_empty_file)
{
    std::filesystem::path empty_file{DATADIR"empty.fa"};
    EXPECT_THROW(jstmap::load_reference(empty_file), std::invalid_argument);
}
