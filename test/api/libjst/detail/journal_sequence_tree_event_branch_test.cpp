// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "journal_sequence_tree_event_test_template.hpp"

INSTANTIATE_TEST_SUITE_P(deletion_event, jst_event_test, testing::Values(jst_event_fixture
{
    .event = {10u, libjst::detail::test::deletion_t{4}, libjst::detail::test::coverage_t{0, 1, 1, 1, 0}},
    .expected_position = 10,
    .category = libjst::detail::test::event_category::branch
}));

INSTANTIATE_TEST_SUITE_P(substitution_event, jst_event_test, testing::Values(jst_event_fixture
{
    .event = {10u, libjst::detail::test::substitution_t{"aaaa"s}, libjst::detail::test::coverage_t{0, 1, 1, 1, 0}},
    .expected_position = 10,
    .category = libjst::detail::test::event_category::branch
}));

INSTANTIATE_TEST_SUITE_P(insertion_event, jst_event_test, testing::Values(jst_event_fixture
{
    .event = {10u, libjst::detail::test::insertion_t{"aaaa"s}, libjst::detail::test::coverage_t{0, 1, 1, 1, 0}},
    .expected_position = 10,
    .category = libjst::detail::test::event_category::branch
}));
