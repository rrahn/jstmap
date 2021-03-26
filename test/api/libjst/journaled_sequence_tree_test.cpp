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

    jst.add(alignment1);
    EXPECT_EQ(jst.size(), 1u);

    jst.add(alignment2);
    EXPECT_EQ(jst.size(), 2u);

    jst.add(alignment3);
    EXPECT_EQ(jst.size(), 3u);

    alignment_t alignment_wrong_reference{libjst::test::make_gapped("aaaabbbbccc-----x"sv), alignment1.second};
    EXPECT_THROW(jst.add(alignment_wrong_reference), std::invalid_argument);

    alignment_t alignment_wrong_order{alignment1.second, alignment1.first};
    EXPECT_THROW(jst.add(alignment_wrong_order), std::invalid_argument);
}

TEST_F(journaled_sequence_tree_fixture, sequence_at)
{
    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    jst.add(alignment1);
    jst.add(alignment2);
    jst.add(alignment3);

    auto target_sequence = [] (auto const & alignment)
    {
        std::string sequence;
        std::ranges::for_each(std::get<1>(alignment), [&] <typename alphabet_t> (seqan3::gapped<alphabet_t> const & c)
        {
            sequence.push_back(seqan3::to_char(c));
        });

        auto const tmp = std::ranges::remove_if(sequence, [] (char const c) { return c == '-'; });
        sequence.erase(tmp.begin(), tmp.end());
        return sequence;
    };

    EXPECT_RANGE_EQ(jst.sequence_at(0),  target_sequence(alignment1));
    EXPECT_RANGE_EQ(jst.sequence_at(1),  target_sequence(alignment2));
    EXPECT_RANGE_EQ(jst.sequence_at(2),  target_sequence(alignment3));
    EXPECT_THROW(jst.sequence_at(3), std::out_of_range);
    EXPECT_THROW(jst.sequence_at(-1), std::out_of_range);
}

TEST_F(journaled_sequence_tree_fixture, context_enumerator)
{
    sequence_t tmp_reference{reference};
    jst_t jst{std::move(tmp_reference)};

    jst.add(alignment1);
    jst.add(alignment2);
    jst.add(alignment3);

    auto context_enumerator = jst.context_enumerator(4u);

    auto advance_supported_context = [&] (auto & it)
    {
        while (it != context_enumerator.end() && it.positions().empty())
            ++it;
    };

    auto it = context_enumerator.begin();
    advance_supported_context(it);
    EXPECT_RANGE_EQ(*it, "ccaa"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "caab"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "aabb"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "aabb"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "abbc"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "bbcc"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "abca"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "bcab"sv);
    advance_supported_context(++it);
    EXPECT_RANGE_EQ(*it, "cabc"sv);
    advance_supported_context(++it);
    EXPECT_TRUE(it == context_enumerator.end());
}

// The test data serialised to disk.
inline constexpr std::string_view expected_output =
R"json({
    "value0": "aaaabbbbcccc",
    "value1": [
        {
            "value0": {
                "value0": 0,
                "value1": {
                    "index": 2,
                    "data": {
                        "value0": {
                            "value0": 12
                        }
                    }
                }
            },
            "value1": [
                true,
                true,
                false
            ]
        },
        {
            "value0": {
                "value0": 12,
                "value1": {
                    "index": 0,
                    "data": {
                        "value0": {
                            "value0": [
                                97,
                                97,
                                98,
                                98,
                                99,
                                99
                            ]
                        }
                    }
                }
            },
            "value1": [
                true,
                false,
                false
            ]
        },
        {
            "value0": {
                "value0": 12,
                "value1": {
                    "index": 0,
                    "data": {
                        "value0": {
                            "value0": [
                                97,
                                98,
                                99,
                                97,
                                98,
                                99
                            ]
                        }
                    }
                }
            },
            "value1": [
                false,
                true,
                false
            ]
        },
        {
            "value0": {
                "value0": 0,
                "value1": {
                    "index": 2,
                    "data": {
                        "value0": {
                            "value0": 4
                        }
                    }
                }
            },
            "value1": [
                false,
                false,
                true
            ]
        },
        {
            "value0": {
                "value0": 4,
                "value1": {
                    "index": 0,
                    "data": {
                        "value0": {
                            "value0": [
                                99,
                                99
                            ]
                        }
                    }
                }
            },
            "value1": [
                false,
                false,
                true
            ]
        },
        {
            "value0": {
                "value0": 4,
                "value1": {
                    "index": 2,
                    "data": {
                        "value0": {
                            "value0": 4
                        }
                    }
                }
            },
            "value1": [
                false,
                false,
                true
            ]
        },
        {
            "value0": {
                "value0": 8,
                "value1": {
                    "index": 0,
                    "data": {
                        "value0": {
                            "value0": [
                                97,
                                97
                            ]
                        }
                    }
                }
            },
            "value1": [
                false,
                false,
                true
            ]
        },
        {
            "value0": {
                "value0": 8,
                "value1": {
                    "index": 2,
                    "data": {
                        "value0": {
                            "value0": 4
                        }
                    }
                }
            },
            "value1": [
                false,
                false,
                true
            ]
        },
        {
            "value0": {
                "value0": 12,
                "value1": {
                    "index": 0,
                    "data": {
                        "value0": {
                            "value0": [
                                98,
                                98
                            ]
                        }
                    }
                }
            },
            "value1": [
                false,
                false,
                true
            ]
        }
    ],
    "value2": 3
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
