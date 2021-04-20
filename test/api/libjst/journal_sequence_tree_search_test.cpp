// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <set>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/search/naive_search.hpp>
#include <libjst/search/state_manager_stack.hpp>
#include <libjst/journaled_sequence_tree.hpp>

#include "test_utility.hpp" // make_gaped

using namespace std::literals;

using sequence_t = seqan3::dna5_vector;

struct jst_search_fixture
{
    std::filesystem::path jst_file{};
    sequence_t pattern{};

    template <typename char_t>
    friend seqan3::debug_stream_type<char_t> & operator<<(seqan3::debug_stream_type<char_t> & stream,
                                                          jst_search_fixture const & fixture)
    {
        return stream << "path: " << fixture.jst_file << " pattern: " << fixture.pattern;
    }
};

struct search_test : public ::testing::TestWithParam<jst_search_fixture>
{
    using jst_t = libjst::journaled_sequence_tree<sequence_t>;

    jst_t jst{};
    std::set<libjst::context_position> expected_hits{};

    void SetUp() override
    {
        jst = libjst::test::load_jst<jst_t>(GetParam().jst_file);

        generate_hits();
    }

    size_t context_size() const noexcept
    {
        return GetParam().pattern.size();
    }

    bool hit_exist(libjst::context_position hit)
    {
        return expected_hits.contains(hit);
    }

private:

    void generate_hits()
    {
        // Iterate over the sequences:
        for (size_t i = 0; i < jst.size(); ++i)
        {
            auto sequence = jst.sequence_at(i);
            size_t last_sequence_pos = sequence.size() - context_size() + 1;
            for (size_t j = 0; j < last_sequence_pos; ++j)
                if (std::ranges::equal(GetParam().pattern, sequence | seqan3::views::slice(j, j + context_size())))
                    expected_hits.emplace(i, j);
        }
    }
};

TEST_P(search_test, search)
{
    libjst::naive_pattern_searcher searcher{GetParam().pattern, libjst::search_state_manager_stack<size_t>{}};
    auto jst_range_agent = jst.range_agent(context_size(), searcher.state_manager()); // already pushing a branch.

    // TODO: The position should be generate from a const iterator as well.
    searcher(jst_range_agent, [&] (auto & it)
    {
        for (libjst::context_position hit : it.positions())
            EXPECT_TRUE(hit_exist(hit)) << "hit: " << hit;
    });
}

// ----------------------------------------------------------------------------
// Test cases
// ----------------------------------------------------------------------------

using seqan3::operator""_dna5;

INSTANTIATE_TEST_SUITE_P(search_with_hit_at_begin, search_test, testing::Values(jst_search_fixture
{
    .jst_file{DATADIR"sim_refx5.jst"},
    .pattern{"TGCGGGACG"_dna5}
}));

INSTANTIATE_TEST_SUITE_P(search_with_hit_at_end, search_test, testing::Values(jst_search_fixture
{
    .jst_file{DATADIR"sim_refx5.jst"},
    .pattern{"GGAGGAATGCT"_dna5}
}));

INSTANTIATE_TEST_SUITE_P(search_with_one_hit_in_middle, search_test, testing::Values(jst_search_fixture
{
    .jst_file{DATADIR"sim_refx5.jst"},
    .pattern{"GGGCGAGAACAACTAATTCCG"_dna5}
}));

INSTANTIATE_TEST_SUITE_P(search_with_hits_in_some_sequences, search_test, testing::Values(jst_search_fixture
{
    .jst_file{DATADIR"sim_refx5.jst"},
    .pattern{"tat"_dna5}
}));

INSTANTIATE_TEST_SUITE_P(search_with_many_hits_in_all_sequences, search_test, testing::Values(jst_search_fixture
{
    .jst_file{DATADIR"sim_refx5.jst"},
    .pattern{"CA"_dna5}
}));

INSTANTIATE_TEST_SUITE_P(search_with_zero_hits, search_test, testing::Values(jst_search_fixture
{
    .jst_file{DATADIR"sim_refx5.jst"},
    .pattern{"GGGCGAGAACAACTAATTCCA"_dna5}
}));
