// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>

#include <span>
#include <string>
#include <vector>

#include <libjst/journal/any_sequence.hpp>

SCENARIO("any_sequence can be constructed for types that are convertible to std::span<char const>", "[any_sequence]")  {
    GIVEN("an any_sequence over a std::span<char const>") {
        using any_sequence_t = libjst::any_sequence<std::span<char const>>;
        WHEN("default initialized") {

            any_sequence_t seq{};

            THEN("seq does not contain a value") {
                REQUIRE(seq.has_value() == false);
                REQUIRE(static_cast<bool>(seq) == false);
            } AND_THEN("accessing seq's value throws libjst::bad_sequence_access") {
                REQUIRE_THROWS_AS(seq.value(), libjst::bad_sequence_access);
            }
        }

        WHEN("initialized with a vector of chars") {
            std::vector<char> source{'A', 'C', 'G', 'T'};
            any_sequence_t seq{source};

            THEN("seq does contain a value") {
                REQUIRE(seq.has_value() == true);
                REQUIRE(static_cast<bool>(seq) == true);
            } AND_THEN("accessing seq's value returns the vector.") {
                REQUIRE(std::ranges::equal(seq.value(), source));
                REQUIRE(std::ranges::equal(*seq, source));
            }
        }

        WHEN("initialized with a std::string") {
            std::string source{"ACGT"};
            any_sequence_t seq{source};

            THEN("seq does contain a value") {
                REQUIRE(seq.has_value() == true);
                REQUIRE(static_cast<bool>(seq) == true);
            } AND_THEN("accessing seq's value returns the string.") {
                REQUIRE(std::ranges::equal(seq.value(), source));
                REQUIRE(std::ranges::equal(*seq, source));
            }
        }
    }
}
