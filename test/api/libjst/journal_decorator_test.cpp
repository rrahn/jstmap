// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/journal_decorator.hpp>

struct journal_decorator_test : public ::testing::Test
{
    std::string sequence{"aaaaccccggggtttt"};
};

// TODO: check iterator concept.

TEST_F(journal_decorator_test, construction)
{
    libjst::journal_decorator<char> jd{};
    EXPECT_EQ(jd.size(), 0u);
    EXPECT_TRUE(jd.empty());

    jd = libjst::journal_decorator{std::span{sequence}};
    EXPECT_EQ(jd.size(), sequence.size());
    EXPECT_FALSE(jd.empty());
}

TEST_F(journal_decorator_test, iterator)
{
    libjst::journal_decorator jd{std::span{sequence}};
    EXPECT_RANGE_EQ(jd, sequence);

    // TODO: Check const iterator
}

TEST_F(journal_decorator_test, record_insertion)
{
    using namespace std::literals;

    std::string segment{"uu"s};

    { // insert in middel
        libjst::journal_decorator jd{std::span{sequence}};

        bool insert_result = jd.record_insertion(8, std::span{segment});

        EXPECT_TRUE(insert_result);
        EXPECT_EQ(jd.size(), sequence.size() + segment.size());

        std::string expected{sequence};
        EXPECT_RANGE_EQ(jd, expected.insert(8, segment));
    }

    // insert in empty
    {
        libjst::journal_decorator<char> jd{};
        bool insert_result = jd.record_insertion(8, std::span{segment});

        EXPECT_FALSE(insert_result);
        EXPECT_EQ(jd.size(), 0u);
        // EXPECT_RANGE_EQ(jd, std::string{}); TODO: check me!

        insert_result = jd.record_insertion(0, std::span{segment});

        EXPECT_TRUE(insert_result);
        EXPECT_EQ(jd.size(), segment.size());
        EXPECT_RANGE_EQ(jd, segment);
    }

    // insert at end
    {
        libjst::journal_decorator jd{std::span{sequence}};
        bool insert_result = jd.record_insertion(sequence.size() + 1, std::span{segment});

        EXPECT_FALSE(insert_result);
        EXPECT_EQ(jd.size(), sequence.size());
        EXPECT_RANGE_EQ(jd, sequence);

        insert_result = jd.record_insertion(sequence.size(), std::span{segment});

        EXPECT_TRUE(insert_result);
        EXPECT_EQ(jd.size(), sequence.size() + segment.size());
        EXPECT_RANGE_EQ(jd, sequence + segment);
    }

    // insert at beginning
    {
        libjst::journal_decorator jd{std::span{sequence}};
        bool insert_result = jd.record_insertion(0, std::span{segment});

        EXPECT_TRUE(insert_result);
        EXPECT_EQ(jd.size(), sequence.size() + segment.size());
        EXPECT_RANGE_EQ(jd, segment + sequence);
    }

    // insert on same position twice.
    {
        libjst::journal_decorator jd{std::span{sequence}};

        bool insert_result = jd.record_insertion(8, std::span{segment});
        EXPECT_TRUE(insert_result);
        EXPECT_EQ(jd.size(), sequence.size() + segment.size());

        insert_result = jd.record_insertion(8, std::span{segment});
        EXPECT_TRUE(insert_result);
        EXPECT_EQ(jd.size(), sequence.size() + (2 * segment.size()));
        std::string expected{sequence};
        EXPECT_RANGE_EQ(jd, expected.insert(8, segment).insert(8, segment));
    }
}

TEST_F(journal_decorator_test, record_insertion_in_empty_journal_sequence)
{
    using namespace std::literals;

    std::string empty_ref{""s};
    std::string single_insertion{"i"s};

    libjst::journal_decorator jd{segment_t{empty_ref}};

    EXPECT_EQ(jd.size(), 0u);
    EXPECT_TRUE(jd.empty());

    bool insert_result = jd.record_insertion(0, segment_t{single_insertion});

    EXPECT_TRUE(insert_result);
    EXPECT_EQ(jd.size(), 1u);
    EXPECT_FALSE(jd.empty());
    EXPECT_RANGE_EQ(jd, "i"s);
}

TEST_F(journal_decorator_test, record_deletion)
{
    using namespace std::literals;

    //-------------------------------------------------------------------------
    // invalid erase

    { // erase from empty journal decorator
        libjst::journal_decorator<char> jd{};

        bool erase_result = jd.record_deletion(0, 10);

        EXPECT_FALSE(erase_result);
        EXPECT_EQ(jd.size(), 0u);
        EXPECT_TRUE(jd.empty());
    }

    { // erase inverted range
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(5, 4);

        EXPECT_FALSE(erase_result);
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());
    }

    { // erase region that goes over the end of the journal decorator
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(5, 17);

        EXPECT_FALSE(erase_result);
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());
    }

    { // erase region that goes over the end of the journal decorator
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(16, 17);

        EXPECT_FALSE(erase_result);
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());
    }

    { // erase empty segment
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(5, 5);

        EXPECT_FALSE(erase_result);
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());
    }

    //-------------------------------------------------------------------------
    // erase from journal decorator with single entry

    { // erase some internal segment
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(4, 8);

        EXPECT_TRUE(erase_result);
        EXPECT_EQ(jd.size(), 12u);
        EXPECT_FALSE(jd.empty());

        std::string expected{sequence};
        EXPECT_RANGE_EQ(jd, expected.erase(4, 4));
    }

    { // erase single element within entry
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(7, 8);

        EXPECT_TRUE(erase_result);
        EXPECT_EQ(jd.size(), 15u);
        EXPECT_FALSE(jd.empty());

        std::string expected{sequence};
        EXPECT_RANGE_EQ(jd, expected.erase(7, 1));
    }

    { // erase entire entry
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(0, 16);

        EXPECT_TRUE(erase_result);
        EXPECT_EQ(jd.size(), 0u);
        EXPECT_TRUE(jd.empty());

        EXPECT_RANGE_EQ(jd, ""s);
    }

    { // erase suffix of entry
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(5, 16);

        EXPECT_TRUE(erase_result);
        EXPECT_EQ(jd.size(), 5u);
        EXPECT_FALSE(jd.empty());

        EXPECT_RANGE_EQ(jd, sequence.substr(0, 5));
    }

    // erase prefix of entry
    {
        libjst::journal_decorator jd{std::span{sequence}};

        bool erase_result = jd.record_deletion(0, 5);

        EXPECT_TRUE(erase_result);
        EXPECT_EQ(jd.size(), 11u);
        EXPECT_FALSE(jd.empty());

        EXPECT_RANGE_EQ(jd, sequence.substr(5, 16));
    }

    //-------------------------------------------------------------------------
    // erase from journal decorator with multiple entries
    // erase between two adjacent entries

    libjst::journal_decorator jd_base{std::span{sequence}};

    EXPECT_TRUE(jd_base.record_deletion(12, 13));
    EXPECT_TRUE(jd_base.record_deletion(8, 9));
    EXPECT_TRUE(jd_base.record_deletion(4, 5));
    EXPECT_TRUE(jd_base.record_deletion(0, 1));

    std::string expected_base{sequence};
    EXPECT_RANGE_EQ(jd_base, expected_base.erase(12, 1).erase(8, 1).erase(4, 1).erase(0, 1));

    { // erase begin from first to last
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(3, 9));
        EXPECT_EQ(jd.size(), 6u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(3, 6));
    }

    { // from last of first to first of last // two elements are delted
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(5, 7));
        EXPECT_EQ(jd.size(), 10u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(5, 2));
    }

    { // from middle to middle - some range is deleted
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(4, 8));
        EXPECT_EQ(jd.size(), 8u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(4, 4));
    }

    { // from middle to end of second - entire second entry is deleted
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(4, 9));
        EXPECT_EQ(jd.size(), 7u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(4, 5));
    }

    { // from begin to middle of second - entire first entry is deleted
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(3, 8));
        EXPECT_EQ(jd.size(), 7u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(3, 5));
    }

    //-------------------------------------------------------------------------
    // erase from journal decorator with multiple entries
    // erase between two distant entries

    { // erase begin from first to last
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(0, 9));
        EXPECT_EQ(jd.size(), 3u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(0, 9));
    }

    { // erase everything
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(0, 12));
        EXPECT_EQ(jd.size(), 0u);
        EXPECT_TRUE(jd.empty());

        EXPECT_RANGE_EQ(jd, ""s);
    }

    { // from last of first to first of last - two elements are delted and all entries in between
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(2, 10));
        EXPECT_EQ(jd.size(), 4u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(2, 8));
    }

    { // from last of first to first of last - two elements are delted and all entries in between
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(1, 11));
        EXPECT_EQ(jd.size(), 2u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(1, 10));
    }

    { // from middle to middle - some elements are delted and all entries in between
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(1, 11));
        EXPECT_EQ(jd.size(), 2u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(1, 10));
    }

    { // from middle to end - some elements are delted and all entries in between
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(1, 12));
        EXPECT_EQ(jd.size(), 1u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(1, 11));
    }

    { // from begin to middle - some elements are delted and all entries in between
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_deletion(0, 11));
        EXPECT_EQ(jd.size(), 1u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.erase(0, 11));
    }
}

TEST_F(journal_decorator_test, record_substitution)
{
    using namespace std::literals;

    std::string segment{"uu"};

    //-------------------------------------------------------------------------
    // replace invalid

    { // replace in empty journal decorator
        libjst::journal_decorator<char> jd{};

        EXPECT_FALSE(jd.record_substitution(0, segment));
        EXPECT_EQ(jd.size(), 0u);
        EXPECT_TRUE(jd.empty());

        EXPECT_FALSE(jd.record_substitution(10, segment));
        EXPECT_EQ(jd.size(), 0u);
        EXPECT_TRUE(jd.empty());
    }

    { // position is after end
        libjst::journal_decorator jd{std::span{sequence}};

        EXPECT_FALSE(jd.record_substitution(17, segment));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());
    }

    { // last position is after end
        libjst::journal_decorator jd{std::span{sequence}};

        EXPECT_FALSE(jd.record_substitution(15, segment));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());
    }

    { // segment is empty
        libjst::journal_decorator jd{std::span{sequence}};

        EXPECT_FALSE(jd.record_substitution(10, std::span<char>{}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());
    }

    //-------------------------------------------------------------------------
    // replace within single entry

    { // replace middle of entry
        libjst::journal_decorator jd{std::span{sequence}};

        EXPECT_TRUE(jd.record_substitution(4, segment));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{sequence};
        EXPECT_RANGE_EQ(jd, expected.replace(4, segment.size(), segment));
    }

    { // replace entire entry
        libjst::journal_decorator jd{std::span{sequence}};

        std::string replace_all(sequence.size(), 'u');
        EXPECT_TRUE(jd.record_substitution(0, std::span{replace_all}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        EXPECT_RANGE_EQ(jd, replace_all);
    }

    { // replace prefix of entry
        libjst::journal_decorator jd{std::span{sequence}};

        EXPECT_TRUE(jd.record_substitution(0, std::span{segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{sequence};
        EXPECT_RANGE_EQ(jd, expected.replace(0, segment.size(), segment));
    }

    { // replace suffix of entry
        libjst::journal_decorator jd{std::span{sequence}};

        EXPECT_TRUE(jd.record_substitution(14, std::span{segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{sequence};
        EXPECT_RANGE_EQ(jd, expected.replace(14, segment.size(), segment));
    }

    //-------------------------------------------------------------------------
    // replace with multiple entries
    // replace between adjacent entries

    libjst::journal_decorator jd_base{std::span{sequence}};
    // create: aa|uu|cc|uu|gg|uu|tt|uu
    jd_base.record_substitution(2, segment);
    jd_base.record_substitution(6, segment);
    jd_base.record_substitution(10, segment);
    jd_base.record_substitution(14, segment);

    std::string expected_base{sequence};
    expected_base.replace(2, segment.size(), segment).
                  replace(6, segment.size(), segment).
                  replace(10, segment.size(), segment).
                  replace(14, segment.size(), segment);

    {
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(5, std::span{segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(5, segment.size(), segment));
    }

    {
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(8, std::span{segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(8, segment.size(), segment));
    }

    { // replace from begin of node to end of next node.
        std::string new_segment{"xxxx"};
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(4, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(4, new_segment.size(), new_segment));
    }

    {
        std::string new_segment{"xxx"};
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(4, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(4, new_segment.size(), new_segment));
    }

    {
        std::string new_segment{"xxx"};
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(5, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(5, new_segment.size(), new_segment));
    }

    {
        std::string new_segment{"x"};
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(5, std::span{new_segment}));
        EXPECT_TRUE(jd.record_substitution(4, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(5, new_segment.size(), new_segment).
                                     replace(4, new_segment.size(), new_segment));
    }

    //-------------------------------------------------------------------------
    // replace with multiple entries
    // replace between distant entries

    { // entire range
        std::string new_segment(sequence.size(), 'y');
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(0, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        EXPECT_RANGE_EQ(jd, new_segment);
    }

    std::string new_segment(sequence.size() - 5, 'y');
    { // replace prefix
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(0, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(0, new_segment.size(), new_segment));
    }

    { // replace suffix
        libjst::journal_decorator jd{jd_base};
        size_t const replace_position = jd.size() - new_segment.size();
        EXPECT_TRUE(jd.record_substitution(replace_position, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(replace_position, new_segment.size(), new_segment));
    }

    { // replace middle from begin
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(3, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(3, new_segment.size(), new_segment));
    }

    {
        libjst::journal_decorator jd{jd_base};

        EXPECT_TRUE(jd.record_substitution(2, std::span{new_segment}));
        EXPECT_EQ(jd.size(), 16u);
        EXPECT_FALSE(jd.empty());

        std::string expected{expected_base};
        EXPECT_RANGE_EQ(jd, expected.replace(2, new_segment.size(), new_segment));
    }
}
