// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/adaptation/char.hpp>

#include <libjst/detail/delta_event_shared.hpp>

struct delta_event_shared_fixture : public ::testing::Test
{
    using delta_event_shared_t = libjst::detail::delta_event_shared<char>;
    using alphabet_t = typename delta_event_shared_t::alphabet_type;
    using insertion_t = typename delta_event_shared_t::insertion_type;
    using substitution_t = typename delta_event_shared_t::substitution_type;
    using deletion_t = typename delta_event_shared_t::deletion_type;
    using delta_event_t = typename delta_event_shared_t::delta_event_type;
    using coverage_t = typename delta_event_shared_t::coverage_type;
};

TEST_F(delta_event_shared_fixture, construction)
{
    EXPECT_TRUE(std::is_default_constructible_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<delta_event_shared_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<delta_event_shared_t>);

    // Special constructor.
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, delta_event_t, coverage_t>));
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, size_t, insertion_t, coverage_t>));
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, size_t, deletion_t, coverage_t>));
    EXPECT_TRUE((std::is_constructible_v<delta_event_shared_t, size_t, substitution_t, coverage_t>));
}

TEST_F(delta_event_shared_fixture, construct_from_substitution)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{10u, substitution_t{"abc"s}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node1.position(), 10u);
    EXPECT_TRUE(node1.is_substitution());

    delta_event_shared_t node2{10u, substitution_t{"abc"s}, coverage_t{0, 1, 1, 0, 0, 1, 1}};
    EXPECT_EQ(node2.position(), 10u);
    EXPECT_TRUE(node2.is_substitution());
}

TEST_F(delta_event_shared_fixture, construct_from_insertion)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{10u, insertion_t{"abc"s}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node1.position(), 10u);
    EXPECT_TRUE(node1.is_insertion());

    delta_event_shared_t node2{10u, insertion_t{"abc"s}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node2.position(), 10u);
    EXPECT_TRUE(node2.is_insertion());
}

TEST_F(delta_event_shared_fixture, construct_from_deletion)
{
    delta_event_shared_t node1{delta_event_t{10u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node1.position(), 10u);
    EXPECT_TRUE(node1.is_deletion());

    delta_event_shared_t node2{10u, deletion_t{1}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node2.position(), 10u);
    EXPECT_TRUE(node2.is_deletion());
}

TEST_F(delta_event_shared_fixture, position)
{
    delta_event_shared_t node{delta_event_t{10u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node.position(), 10u);
    EXPECT_EQ(std::as_const(node).position(), 10u);
}

TEST_F(delta_event_shared_fixture, coverage)
{
    delta_event_shared_t node{delta_event_t{10u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_EQ(node.coverage(), (coverage_t{0, 1, 1, 0, 0, 1, 1}));
    EXPECT_EQ(std::as_const(node).coverage(), (coverage_t{0, 1, 1, 0, 0, 1, 1}));
}

TEST_F(delta_event_shared_fixture, equality)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{10u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node2{delta_event_t{10u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node3{delta_event_t{9u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node4{delta_event_t{9u, substitution_t{"a"s}}, coverage_t{0, 1, 1, 0, 0, 1, 0}};

    EXPECT_EQ(std::as_const(node1), std::as_const(node1));
    EXPECT_NE(std::as_const(node1), std::as_const(node2));
    EXPECT_NE(std::as_const(node1), std::as_const(node3));
    EXPECT_NE(std::as_const(node1), std::as_const(node4));

    EXPECT_NE(std::as_const(node2), std::as_const(node1));
    EXPECT_EQ(std::as_const(node2), std::as_const(node2));
    EXPECT_NE(std::as_const(node2), std::as_const(node3));
    EXPECT_NE(std::as_const(node2), std::as_const(node4));

    EXPECT_NE(std::as_const(node3), std::as_const(node1));
    EXPECT_NE(std::as_const(node3), std::as_const(node2));
    EXPECT_EQ(std::as_const(node3), std::as_const(node3));
    EXPECT_NE(std::as_const(node3), std::as_const(node4));

    EXPECT_NE(std::as_const(node4), std::as_const(node1));
    EXPECT_NE(std::as_const(node4), std::as_const(node2));
    EXPECT_NE(std::as_const(node4), std::as_const(node3));
    EXPECT_EQ(std::as_const(node4), std::as_const(node4));
}

TEST_F(delta_event_shared_fixture, less)
{
    using namespace std::literals;

    delta_event_shared_t node1{delta_event_t{0u, deletion_t{2}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node2{delta_event_t{0u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node3{delta_event_t{0u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node4{delta_event_t{0u, substitution_t{"abc"s}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node5{delta_event_t{10u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};
    delta_event_shared_t node6{delta_event_t{10u, deletion_t{1}}, coverage_t{0, 1, 1, 0, 0, 1, 0}};
    delta_event_shared_t node7{delta_event_t{10u, substitution_t{"abc"s}}, coverage_t{0, 1, 1, 0, 0, 1, 1}};

    EXPECT_FALSE(std::as_const(node1) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node1) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node1) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node1) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node1) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node1) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node1) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node2) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node2) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node2) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node2) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node2) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node2) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node2) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node3) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node3) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node3) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node3) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node3) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node3) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node3) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node4) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node4) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node4) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node4) < std::as_const(node4));
    EXPECT_TRUE (std::as_const(node4) < std::as_const(node5));
    EXPECT_TRUE (std::as_const(node4) < std::as_const(node6));
    EXPECT_TRUE (std::as_const(node4) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node5) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node4));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node5));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node6));
    EXPECT_FALSE(std::as_const(node5) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node6) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node4));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node5));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node6));
    EXPECT_FALSE(std::as_const(node6) < std::as_const(node7));

    EXPECT_FALSE(std::as_const(node7) < std::as_const(node2));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node3));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node1));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node4));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node5));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node6));
    EXPECT_FALSE(std::as_const(node7) < std::as_const(node7));
}

TEST_F(delta_event_shared_fixture, fill_vector)
{
    std::vector nodes
    {
        delta_event_shared_t{delta_event_t{10u, deletion_t{1}}, coverage_t{1, 0}},
        delta_event_shared_t{delta_event_t{11u, deletion_t{1}}, coverage_t{0, 1}},
    };

    nodes.emplace_back(12u, insertion_t{"xxx"}, coverage_t{1, 1});

    EXPECT_EQ(nodes.size(), 3u);
}

TEST_F(delta_event_shared_fixture, formatted_output)
{
    using namespace std::literals;

    delta_event_shared_t del{10u, deletion_t{1}, coverage_t{1, 0, 1, 0}};
    delta_event_shared_t ins{11u, insertion_t{"ii"s}, coverage_t{1, 1, 0, 0}};
    delta_event_shared_t sub{12u, substitution_t{"sss"s}, coverage_t{0, 0, 1, 1}};

    std::stringstream sstream;

    sstream << del << '\n' << ins << '\n' << sub;
    EXPECT_EQ(sstream.str(), ("(10, del: 1) ~ <1010>\n(11, ins: ii) ~ <1100>\n(12, sub: sss) ~ <0011>"s));
}
