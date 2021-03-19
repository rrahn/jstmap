// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>

#include <libjst/detail/transform_to_delta_events.hpp>

struct transform_to_delta_events_test : public ::testing::Test
{
    using delta_event_t = libjst::detail::delta_event<char>;

    std::string base_sequence{"aaaaccccggggtttt"};

    auto make_alignment(std::string sequence1, std::string sequence2)
    {
        EXPECT_EQ(sequence1.size(), sequence2.size());

        using gapped_sequence_t = std::vector<seqan3::gapped<char>>;

        gapped_sequence_t aligned_reference{};
        gapped_sequence_t aligned_target{};

        for (auto it1 = sequence1.begin(), it2 = sequence2.begin(); it1 < sequence1.end(); ++it1, ++it2)
        {
            if (*it1 == '-')
            {
                EXPECT_TRUE(*it2 != '-');
                aligned_reference.emplace_back(seqan3::gap{});
                aligned_target.emplace_back(*it2);
            }
            else if (*it2 == '-')
            {
                EXPECT_TRUE(*it1 != '-');
                aligned_reference.emplace_back(*it1);
                aligned_target.emplace_back(seqan3::gap{});
            }
            else
            {
                aligned_reference.emplace_back(*it1);
                aligned_target.emplace_back(*it2);
            }
        }

        EXPECT_EQ(aligned_reference.size(), aligned_target.size());

        return std::pair{aligned_reference, aligned_target};
    }
};

TEST_F(transform_to_delta_events_test, substitutions)
{
    using namespace std::literals;

    { // at begin
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence, target.replace(0, 4, "gggg"s));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{0, libjst::detail::delta_kind_substitution{"gggg"s}}));
    }

    { // in middle
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence, target.replace(4, 4, "gggg"s));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{4, libjst::detail::delta_kind_substitution{"gggg"s}}));
    }

    { // at end
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence, target.replace(12, 4, "gggg"s));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{12, libjst::detail::delta_kind_substitution{"gggg"s}}));
    }

    { // multiple
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence,
                                        target.replace(0, 2, "gg"s).replace(3, 2, "tt"s).replace(6, 2, "aa"));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 3u);
        EXPECT_EQ(delta_events[0], (delta_event_t{0, libjst::detail::delta_kind_substitution{"gg"s}}));
        EXPECT_EQ(delta_events[1], (delta_event_t{3, libjst::detail::delta_kind_substitution{"tt"s}}));
        EXPECT_EQ(delta_events[2], (delta_event_t{6, libjst::detail::delta_kind_substitution{"aa"s}}));
    }
}

TEST_F(transform_to_delta_events_test, deletions)
{
    using namespace std::literals;

    { // at begin
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence, target.replace(0, 4, "----"s));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{0, libjst::detail::delta_kind_deletion{4}}));
    }

    { // in middle
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence, target.replace(5, 4, "----"s));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{5, libjst::detail::delta_kind_deletion{4}}));
    }

    { // at end
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence, target.replace(11, 5, "-----"s));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{11, libjst::detail::delta_kind_deletion{5}}));
    }

    { // multiple
        std::string target{base_sequence};
        auto alignment = make_alignment(base_sequence,
                                        target.replace(0, 2, "--"s).replace(3, 2, "--"s).replace(6, 3, "---"));
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 3u);
        EXPECT_EQ(delta_events[0], (delta_event_t{0, libjst::detail::delta_kind_deletion{2}}));
        EXPECT_EQ(delta_events[1], (delta_event_t{3, libjst::detail::delta_kind_deletion{2}}));
        EXPECT_EQ(delta_events[2], (delta_event_t{6, libjst::detail::delta_kind_deletion{3}}));
    }
}

TEST_F(transform_to_delta_events_test, insertions)
{
    using namespace std::literals;

    { // at begin
        std::string reference{base_sequence};
        auto alignment = make_alignment(reference.replace(0, 3, "---"s), base_sequence);
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{0, libjst::detail::delta_kind_insertion{"aaa"s}}));
    }

    { // in middle
        std::string reference{base_sequence};
        auto alignment = make_alignment(reference.replace(4, 3, "---"s), base_sequence);
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{4, libjst::detail::delta_kind_insertion{"ccc"s}}));
    }

    { // at end
        std::string reference{base_sequence};
        auto alignment = make_alignment(reference.replace(13, 3, "---"s), base_sequence);
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 1u);
        EXPECT_EQ(delta_events[0], (delta_event_t{13, libjst::detail::delta_kind_insertion{"ttt"s}}));
    }

    { // multiple
        std::string reference{base_sequence};
        auto alignment = make_alignment(reference.replace(0, 3, "---"s).replace(4, 1, "-"s).replace(11, 5, "-----"s),
                                        base_sequence);
        auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

        EXPECT_EQ(delta_events.size(), 3u);
        EXPECT_EQ(delta_events[0], (delta_event_t{0, libjst::detail::delta_kind_insertion{"aaa"s}}));
        EXPECT_EQ(delta_events[1], (delta_event_t{1, libjst::detail::delta_kind_insertion{"c"s}}));
        EXPECT_EQ(delta_events[2], (delta_event_t{7, libjst::detail::delta_kind_insertion{"gtttt"s}}));
    }
}

TEST_F(transform_to_delta_events_test, mixed)
{
    using namespace std::literals;

    // -- aa g - cc ggg - t tt t
    // aa aa c c -- ggg g a -- t

    std::string reference{base_sequence};
    std::string target{base_sequence};
    reference.replace(0, 2, "--"s).replace(4, 1, "g"s).replace(5, 1, "-"s).replace(11, 1, "-"s);
    target.replace(6, 2, "--").replace(6, 2, "--"s).replace(12, 1, "a"s).replace(13, 2, "--"s);

    EXPECT_EQ(reference, "--aag-ccggg-tttt");
    EXPECT_EQ(target,    "aaaacc--gggga--t");
    auto alignment = make_alignment(reference, target);
    auto delta_events = libjst::detail::transform_to_delta_events<char>(alignment);

    EXPECT_EQ(delta_events.size(), 7u);
    EXPECT_EQ(delta_events[0], (delta_event_t{0, libjst::detail::delta_kind_insertion{"aa"s}}));
    EXPECT_EQ(delta_events[1], (delta_event_t{2, libjst::detail::delta_kind_substitution{"c"s}}));
    EXPECT_EQ(delta_events[2], (delta_event_t{3, libjst::detail::delta_kind_insertion{"c"s}}));
    EXPECT_EQ(delta_events[3], (delta_event_t{3, libjst::detail::delta_kind_deletion{2}}));
    EXPECT_EQ(delta_events[4], (delta_event_t{8, libjst::detail::delta_kind_insertion{"g"s}}));
    EXPECT_EQ(delta_events[5], (delta_event_t{8, libjst::detail::delta_kind_substitution{"a"s}}));
    EXPECT_EQ(delta_events[6], (delta_event_t{9, libjst::detail::delta_kind_deletion{2}}));
}
