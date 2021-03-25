// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <cereal/archives/json.hpp> // output archive for testing

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include <libjst/detail/delta_event.hpp>

struct delta_operation_fixture : public ::testing::Test
{
    using delta_event_t = libjst::detail::delta_event<char>;

    static constexpr std::string_view expected_substitution_archive =
R"json({
    "value0": 23,
    "value1": {
        "index": 1,
        "data": {
            "value0": {
                "value0": [
                    97,
                    98,
                    99,
                    100
                ]
            }
        }
    }
})json";

    static constexpr std::string_view expected_insertion_archive =
R"json({
    "value0": 5,
    "value1": {
        "index": 0,
        "data": {
            "value0": {
                "value0": [
                    105,
                    106,
                    107,
                    108,
                    109
                ]
            }
        }
    }
})json";

    static constexpr std::string_view expected_deletion_archive =
R"json({
    "value0": 100,
    "value1": {
        "index": 2,
        "data": {
            "value0": {
                "value0": 10
            }
        }
    }
})json";
};

TEST_F(delta_operation_fixture, basic_construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<delta_event_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<delta_event_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<delta_event_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<delta_event_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<delta_event_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<delta_event_t>);
}

TEST_F(delta_operation_fixture, construct_from_substitution)
{
    using namespace std::literals;

    delta_event_t op{10u, libjst::detail::delta_kind_substitution{"abc"s}};

    EXPECT_EQ(op.position(), 10u);
    EXPECT_TRUE(op.is_substitution());
}

TEST_F(delta_operation_fixture, construct_from_insertion)
{
    using namespace std::literals;

    delta_event_t op{10u, libjst::detail::delta_kind_insertion{"abc"s}};

    EXPECT_EQ(op.position(), 10u);
    EXPECT_TRUE(op.is_insertion());
}

TEST_F(delta_operation_fixture, construct_from_deletion)
{
    delta_event_t op{10u, libjst::detail::delta_kind_deletion{3u}};

    EXPECT_EQ(op.position(), 10u);
    EXPECT_TRUE(op.is_deletion());
}

TEST_F(delta_operation_fixture, deletion_size)
{
    using namespace std::literals;

    { // substitution
        delta_event_t op{10u, libjst::detail::delta_kind_substitution{"abc"s}};

        EXPECT_EQ(op.deletion_size(), 3u);
    }

    { // insertion
        delta_event_t op{10u, libjst::detail::delta_kind_insertion{"abc"s}};

        EXPECT_EQ(op.deletion_size(), 0u);
    }

    { // deletion
        delta_event_t op{10u, libjst::detail::delta_kind_deletion{3u}};

        EXPECT_EQ(op.deletion_size(), 3u);
    }
}

TEST_F(delta_operation_fixture, insertion_size)
{
    using namespace std::literals;

    { // substitution
        delta_event_t op{10u, libjst::detail::delta_kind_substitution{"abc"s}};

        EXPECT_EQ(op.insertion_size(), 3u);
    }

    { // insertion
        delta_event_t op{10u, libjst::detail::delta_kind_insertion{"abc"s}};

        EXPECT_EQ(op.insertion_size(), 3u);
    }

    { // deletion
        delta_event_t op{10u, libjst::detail::delta_kind_deletion{3u}};

        EXPECT_EQ(op.insertion_size(), 0u);
    }
}

TEST_F(delta_operation_fixture, sequence)
{
    using namespace std::literals;

    { // substitution
        delta_event_t op{10u, libjst::detail::delta_kind_substitution{"abc"s}};

        EXPECT_RANGE_EQ(op.sequence(), "abc"s);
    }

    { // insertion
        delta_event_t op{10u, libjst::detail::delta_kind_insertion{"abc"s}};

        EXPECT_RANGE_EQ(op.sequence(), "abc"s);
    }

    { // deletion
        delta_event_t op{10u, libjst::detail::delta_kind_deletion{3u}};

        EXPECT_RANGE_EQ(op.sequence(), ""s);
    }
}

TEST_F(delta_operation_fixture, equality)
{
    using namespace std::literals;

    delta_event_t op1{10u, libjst::detail::delta_kind_deletion{3u}};
    delta_event_t op2{10u, libjst::detail::delta_kind_deletion{2u}};
    delta_event_t op3{9u, libjst::detail::delta_kind_deletion{3u}};
    delta_event_t op4{9u, libjst::detail::delta_kind_substitution{"3u"s}};
    delta_event_t op5{9u, libjst::detail::delta_kind_insertion{"3u"s}};

    EXPECT_EQ(op1, op1);
    EXPECT_NE(op1, op2);
    EXPECT_NE(op1, op3);
    EXPECT_NE(op1, op4);
    EXPECT_NE(op1, op5);

    EXPECT_NE(op2, op1);
    EXPECT_EQ(op2, op2);
    EXPECT_NE(op2, op3);
    EXPECT_NE(op2, op4);
    EXPECT_NE(op2, op5);

    EXPECT_NE(op3, op1);
    EXPECT_NE(op3, op2);
    EXPECT_EQ(op3, op3);
    EXPECT_NE(op3, op4);
    EXPECT_NE(op3, op5);

    EXPECT_NE(op4, op1);
    EXPECT_NE(op4, op2);
    EXPECT_NE(op4, op3);
    EXPECT_EQ(op4, op4);
    EXPECT_NE(op4, op5);

    EXPECT_NE(op5, op1);
    EXPECT_NE(op5, op2);
    EXPECT_NE(op5, op3);
    EXPECT_NE(op5, op4);
    EXPECT_EQ(op5, op5);
}

TEST_F(delta_operation_fixture, stream)
{
    using namespace std::literals;

    {
        delta_event_t op{10u, libjst::detail::delta_kind_substitution{"abc"s}};
        std::stringstream sstream;
        sstream << op;
        EXPECT_EQ(sstream.str(), "(10, sub: abc)"s);
    }

    {
        delta_event_t op{10u, libjst::detail::delta_kind_insertion{"abc"s}};
        std::stringstream sstream;
        sstream << op;
        EXPECT_EQ(sstream.str(), "(10, ins: abc)"s);
    }

    {
        delta_event_t op{10u, libjst::detail::delta_kind_deletion{3u}};
        std::stringstream sstream;
        sstream << op;
        EXPECT_EQ(sstream.str(), "(10, del: 3)"s);
    }
}

TEST_F(delta_operation_fixture, save_substitution)
{
    using namespace std::literals;
    std::stringstream archive_stream{};

    delta_event_t substitution_event{23u, libjst::detail::delta_kind_substitution{"abcd"s}};

    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        substitution_event.save(output_archive);
    }

    EXPECT_EQ(archive_stream.str(), expected_substitution_archive);
}

TEST_F(delta_operation_fixture, save_insertion)
{
    using namespace std::literals;
    std::stringstream archive_stream{};

    delta_event_t insertion_event{5u, libjst::detail::delta_kind_insertion{"ijklm"s}};

    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        insertion_event.save(output_archive);
    }

    EXPECT_EQ(archive_stream.str(), expected_insertion_archive);
}

TEST_F(delta_operation_fixture, save_deletion)
{
    std::stringstream archive_stream{};

    delta_event_t deletion_event{100u, libjst::detail::delta_kind_deletion{10}};

    {
        cereal::JSONOutputArchive output_archive(archive_stream);
        deletion_event.save(output_archive);
    }

    EXPECT_EQ(archive_stream.str(), expected_deletion_archive);
}

TEST_F(delta_operation_fixture, load_substitution)
{
    using namespace std::literals;

    std::stringstream archive_stream{expected_substitution_archive.data()};

    delta_event_t substitution_event{};
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        substitution_event.load(input_archive);
    }

    EXPECT_EQ(substitution_event, (delta_event_t{23u, libjst::detail::delta_kind_substitution{"abcd"s}}));
}

TEST_F(delta_operation_fixture, load_insertion)
{
    using namespace std::literals;

    std::stringstream archive_stream{expected_insertion_archive.data()};

    delta_event_t insertion_event{};
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        insertion_event.load(input_archive);
    }

    EXPECT_EQ(insertion_event, (delta_event_t{5u, libjst::detail::delta_kind_insertion{"ijklm"s}}));
}

TEST_F(delta_operation_fixture, load_deletion)
{
    std::stringstream archive_stream{expected_deletion_archive.data()};

    delta_event_t deletion_event{};
    {
        cereal::JSONInputArchive input_archive(archive_stream);
        deletion_event.load(input_archive);
    }

    EXPECT_EQ(deletion_event, (delta_event_t{100u, libjst::detail::delta_kind_deletion{10}}));
}
