// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

// Including required debug stream types for range_eq test.
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <jstmap/create/load_sequence.hpp>

#include "test_data.hpp"

TEST(jstmap_index, load_sequence)
{
    using spm::operator""_dna4;

    std::filesystem::path input_file{DATADIR"in.fasta"};
    std::vector actual_sequences = jstmap::load_sequences(input_file);

    std::vector expected_sequences{
        "ACGTTTGATTCGCG"_dna4,
        "TCGGGGGATTCGCG"_dna4
    };

    EXPECT_RANGE_EQ(actual_sequences, expected_sequences);
}
