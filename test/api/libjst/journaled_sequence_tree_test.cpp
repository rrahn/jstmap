// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <cereal/archives/json.hpp> // output archive for testing
#include <cereal/types/string.hpp> // sereialise std::string

#include <seqan3/alphabet/adaptation/char.hpp> // allow std::string be recognised as seqan3::sequence
#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/journaled_sequence_tree.hpp>

#include "test_utility.hpp" // make_gaped

using namespace std::literals;

struct journaled_sequence_tree_fixture : public ::testing::Test
{
    using sequence_t = std::string;
    using jst_t = libjst::journaled_sequence_tree<sequence_t>;

    using aligned_sequence_t = std::vector<seqan3::gapped<char>>;
    using alignment_t = std::pair<aligned_sequence_t, aligned_sequence_t>;

    sequence_t reference{"aaaabbbbcccc"};

    alignment_t alignment1{libjst::test::make_gapped("aaaabbbbcccc------"sv),
                           libjst::test::make_gapped("------------aabbcc"sv)};
    alignment_t alignment2{libjst::test::make_gapped("aaaabbbbcccc------"sv),
                           libjst::test::make_gapped("------------abcabc"sv)};
    alignment_t alignment3{libjst::test::make_gapped("aaaa--bbbb--cccc--"sv),
                           libjst::test::make_gapped("----cc----aa----bb"sv)};
};

TEST_F(journaled_sequence_tree_fixture, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<jst_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<jst_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<jst_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<jst_t>);
    EXPECT_TRUE((std::is_constructible_v<jst_t, sequence_t>)); // Set explicit reference sequence.
    EXPECT_TRUE((std::is_constructible_v<jst_t, sequence_t &&>)); // Set explicit reference sequence.
    // I might be able to convert a journaled sequence to an alignment.
    // TODO: Construct from multiple sequence alignment: reference is consensus sequence
}

TEST_F(journaled_sequence_tree_fixture, reference)
{
    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    EXPECT_EQ(jst.reference(), reference);
}

TEST_F(journaled_sequence_tree_fixture, size)
{
    jst_t jst{std::move(reference)};

    EXPECT_EQ(jst.size(), 0u);
}

TEST_F(journaled_sequence_tree_fixture, add)
{
    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    jst.add(alignment1); // Yes we can verify that the first sequence is the
    EXPECT_EQ(jst.size(), 1u);

    jst.add(alignment2); // Yes we can verify that the first sequence is the
    EXPECT_EQ(jst.size(), 2u);

    jst.add(alignment3); // Yes we can verify that the first sequence is the
    EXPECT_EQ(jst.size(), 3u);

    alignment_t alignment_wrong_reference{libjst::test::make_gapped("aaaabbbbccc-----x"sv), alignment1.second};
    EXPECT_THROW(jst.add(alignment_wrong_reference), std::invalid_argument);

    alignment_t alignment_wrong_order{alignment1.second, alignment1.first};
    EXPECT_THROW(jst.add(alignment_wrong_order), std::invalid_argument);
}

TEST_F(journaled_sequence_tree_fixture, context_enumerator)
{
    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    jst.add(alignment1);
    jst.add(alignment2);
    jst.add(alignment3);

    auto context_enumerator = jst.context_enumerator(4u);

    auto it = context_enumerator.begin();
    EXPECT_RANGE_EQ(*it, "aabb"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "abbc"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "bbcc"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "abca"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "bcab"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "cabc"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "ccaa"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "caab"sv);
    ++it;
    EXPECT_RANGE_EQ(*it, "aabb"sv);
    ++it;
    EXPECT_TRUE(it == context_enumerator.end());
}

// The test data serialised to disk.
inline constexpr std::string_view expected_output =
R"json({
    "value0": "aaaabbbbcccc",
    "value1": [
        "aabbcc",
        "abcabc",
        "ccaabb"
    ]
})json";

TEST_F(journaled_sequence_tree_fixture, save)
{
    std::stringstream output_stream{};

    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    jst.add(alignment1);
    jst.add(alignment2);
    jst.add(alignment3);

    {
        cereal::JSONOutputArchive output_archive(output_stream);
        jst.save(output_archive);
    }

    EXPECT_EQ(output_stream.str(), expected_output);
}

TEST_F(journaled_sequence_tree_fixture, load)
{
    std::stringstream archive_stream{expected_output.data()};
    jst_t jst{};

    {
        cereal::JSONInputArchive input_archive(archive_stream);
        jst.load(input_archive);
    }

    EXPECT_EQ(jst.size(), 3u);
    EXPECT_EQ(jst.reference(), "aaaabbbbcccc"sv);
}
