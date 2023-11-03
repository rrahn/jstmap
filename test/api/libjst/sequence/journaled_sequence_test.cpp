// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>

#include <iostream>
#include <string>
#include <vector>

#include <libjst/sequence/journaled_sequence.hpp>

namespace test
{
    template <std::ranges::forward_range t>
    struct debug_sequence
    {
        t value;

        template <typename char_t, typename char_traits_t>
        friend std::basic_ostream<char_t, char_traits_t> & operator<<(std::basic_ostream<char_t, char_traits_t> & stream, debug_sequence const & seq)
        {
            for (auto && c : seq.value)
                stream << c;
            return stream;
        }
    };
}

SCENARIO("A journaled sequence can be initialized", "[sequence][journaled_sequence]")
{
    GIVEN("A default initialized journaled sequence")
    {
        libjst::journaled_sequence<std::vector<char>> journaled_sequence{};
        THEN("The journaled sequence is empty")
        {
            CHECK(journaled_sequence.empty());
            REQUIRE(journaled_sequence.size() == 0u);
        }
        AND_THEN("The journaled sequence models std::ranges::random_access_range")
        {
            REQUIRE(std::ranges::random_access_range<decltype(journaled_sequence)>);
        }
        WHEN("insert() is called")
        {
            std::vector sequence{'A', 'C', 'G', 'T'};
            auto it = journaled_sequence.insert(journaled_sequence.begin(), sequence);
            THEN("The journaled sequence is not empty and its range spells the inserted sequence")
            {
                CHECK_FALSE(journaled_sequence.empty());
                REQUIRE(journaled_sequence.size() == sequence.size());
                INFO("actual " << test::debug_sequence{journaled_sequence});
                INFO("expected " << test::debug_sequence{sequence});
                REQUIRE(std::ranges::equal(journaled_sequence, sequence));
            }
            AND_THEN("The returned iterator is the begin of the inserted sequence")
            {
                REQUIRE(it == journaled_sequence.begin());
            }
        }
        AND_WHEN("erase() is called")
        {
            auto it = journaled_sequence.erase(journaled_sequence.begin(), journaled_sequence.end());
            THEN("The journaled sequence is not changed")
            {
                CHECK(journaled_sequence.empty());
                REQUIRE(journaled_sequence.size() == 0u);
            }
            AND_THEN("The returned iterator is the end of the journaled sequence")
            {
                REQUIRE(it == journaled_sequence.end());
            }
        }
    }
}

SCENARIO("Modifying a journaled sequence", "[sequence][journaled_sequence]")
{
    using namespace std::literals;

    GIVEN("A non-empty journaled sequence initialized with some sequence")
    {
        std::vector sequence{'A', 'C', 'G', 'T'};
        libjst::journaled_sequence journaled_sequence{sequence};
        THEN("The journaled sequence is not empty and its range spells the original sequence")
        {
            CHECK_FALSE(journaled_sequence.empty());
            REQUIRE(journaled_sequence.size() == sequence.size());
            REQUIRE(std::ranges::equal(journaled_sequence, sequence));
        }
        GIVEN("A second sequence to insert")
        {
            std::vector insert_sequence{'T', 'G', 'C', 'A'};
            WHEN("This sequence is inserted in the middle of the journaled sequence")
            {
                auto it = journaled_sequence.insert(journaled_sequence.begin() + 2, insert_sequence);
                THEN("The journaled sequence represents the modified sequence")
                {
                    std::vector expected_sequence{'A', 'C', 'T', 'G', 'C', 'A', 'G', 'T'};
                    REQUIRE(std::ranges::equal(journaled_sequence, expected_sequence));
                    REQUIRE(journaled_sequence.size() == expected_sequence.size());
                }
                AND_THEN("The returned iterator is the begin of the inserted sequence")
                {
                    REQUIRE(it == journaled_sequence.begin() + 2);
                }
                AND_THEN("The original sequence is not changed")
                {
                    REQUIRE(std::ranges::equal(sequence, "ACGT"s));
                }
            }
            WHEN("This sequence is inserted at the beginning of the journaled sequence")
            {
                auto it = journaled_sequence.insert(journaled_sequence.begin(), insert_sequence);
                THEN("The journaled sequence represents the modified sequence")
                {
                    std::vector expected_sequence{'T', 'G', 'C', 'A', 'A', 'C', 'G', 'T'};
                    REQUIRE(std::ranges::equal(journaled_sequence, expected_sequence));
                    REQUIRE(journaled_sequence.size() == expected_sequence.size());
                }
                AND_THEN("The returned iterator is the begin of the inserted sequence")
                {
                    REQUIRE(it == journaled_sequence.begin());
                }
                AND_THEN("The original sequence is not changed")
                {
                    REQUIRE(std::ranges::equal(sequence, "ACGT"s));
                }
            }
            WHEN("This sequence is inserted at the end of the journaled sequence")
            {
                auto it = journaled_sequence.insert(journaled_sequence.end(), insert_sequence);
                THEN("The journaled sequence represents the modified sequence")
                {
                    std::vector expected_sequence{'A', 'C', 'G', 'T', 'T', 'G', 'C', 'A'};
                    REQUIRE(std::ranges::equal(journaled_sequence, expected_sequence));
                    REQUIRE(journaled_sequence.size() == expected_sequence.size());
                }
                AND_THEN("The returned iterator is the begin of the inserted sequence")
                {
                    REQUIRE(it == journaled_sequence.end() - insert_sequence.size());
                }
                AND_THEN("The original sequence is not changed")
                {
                    REQUIRE(std::ranges::equal(sequence, "ACGT"s));
                }
            }
        }
        WHEN("A single element is erased in the middle")
        {
            auto it = journaled_sequence.erase(journaled_sequence.begin() + 2);
            THEN("The journaled sequence represents the modified sequence")
            {
                std::vector expected_sequence{'A', 'C', 'T'};
                REQUIRE(std::ranges::equal(journaled_sequence, expected_sequence));
                REQUIRE(journaled_sequence.size() == expected_sequence.size());
            }
            AND_THEN("The returned iterator points to first element after the deletion") {
                REQUIRE(it == journaled_sequence.begin() + 2);
            }
            AND_THEN("The original sequence is not changed") {
                REQUIRE(std::ranges::equal(sequence, "ACGT"s));
            }
        }
        WHEN("A range is erased")
        {
            auto it = journaled_sequence.erase(journaled_sequence.begin() + 1, journaled_sequence.begin() + 3);
            THEN("The journaled sequence represents the modified sequence")
            {
                std::vector expected_sequence{'A', 'T'};
                REQUIRE(std::ranges::equal(journaled_sequence, expected_sequence));
                REQUIRE(journaled_sequence.size() == expected_sequence.size());
            }
            AND_THEN("The returned iterator points to first element after the deletion") {
                REQUIRE(*it == *(journaled_sequence.begin() + 1));
            }
            AND_THEN("The original sequence is not changed") {
                REQUIRE(std::ranges::equal(sequence, "ACGT"s));
            }
        }
        GIVEN("A second sequence")
        {
            std::vector replace_sequence{'T', 'G', 'C', 'A'};
            WHEN("Replacing the entire journaled sequence with the second sequence")
            {
                journaled_sequence.replace(journaled_sequence.begin(), journaled_sequence.end(), replace_sequence);
                THEN("The journaled sequence represents the second sequence")
                {
                    REQUIRE(std::ranges::equal(journaled_sequence, replace_sequence));
                    REQUIRE(journaled_sequence.size() == replace_sequence.size());
                }
                AND_THEN("The original sequence is not changed") {
                    REQUIRE(std::ranges::equal(sequence, "ACGT"s));
                }
            }
            WHEN("Replacing a segment of the journaled sequence with the second sequence")
            {
                auto it = journaled_sequence.replace(journaled_sequence.begin() + 1, journaled_sequence.begin() + 3, replace_sequence);
                THEN("The journaled sequence represents the modified sequence")
                {
                    std::vector expected_sequence{'A', 'T', 'G', 'C', 'A', 'T'};
                    REQUIRE(std::ranges::equal(journaled_sequence, expected_sequence));
                    REQUIRE(journaled_sequence.size() == expected_sequence.size());
                }
                AND_THEN("The returned iterator points to first element after the replacement") {
                    REQUIRE(it == journaled_sequence.begin() + 1);
                }
                AND_THEN("The original sequence is not changed") {
                    REQUIRE(std::ranges::equal(sequence, "ACGT"s));
                }
            }
        }
    }
}
