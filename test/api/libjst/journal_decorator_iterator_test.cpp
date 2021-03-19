// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>
#include <seqan3/test/../../../unit/range/iterator_test_template.hpp>

#include <libjst/journal_decorator.hpp>

using segment_t = std::span<char>;
using journal_decorator_t = libjst::journal_decorator<segment_t>;
using journal_decorator_iterator = std::ranges::iterator_t<journal_decorator_t>;

template <>
struct iterator_fixture<journal_decorator_iterator> : public ::testing::Test
{
    std::string reference{"aaaaaaaa"};
    std::string ins_segment{"ccccgggggggg"};
    std::string repl_segment{"tttt"};

    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    journal_decorator_t test_range{};
    std::string expected_range{"aaaaccccggggtttt"};

    void SetUp() override
    {
        test_range = journal_decorator_t{std::span{reference}}; // aaaa
        test_range.record_insertion(4, std::span{ins_segment}); // aaaaccccggggggggaaaa
        test_range.record_substitution(16, std::span{repl_segment}); // aaaaccccggggggggtttt
        test_range.record_deletion(9, 13); // aaaaccccggggtttt

        EXPECT_RANGE_EQ(test_range, expected_range);
    }
};

INSTANTIATE_TYPED_TEST_SUITE_P(journal_decorator_iterator_test,
                               iterator_fixture,
                               ::testing::Types<journal_decorator_iterator>, );
