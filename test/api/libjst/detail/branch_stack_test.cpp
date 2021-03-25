// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <concepts>
#include <type_traits>

#include <libjst/detail/branch_stack.hpp>

namespace libjst::test
{

class branch
{
    int _value{};
public:

    constexpr branch() = default;
    explicit constexpr branch(int value) noexcept : _value{value}
    {}

    constexpr operator int() const noexcept
    {
        return _value;
    }
};

using branch_stack_t = libjst::detail::branch_stack<branch>;

} // namespace libjst::test

TEST(branch_stack_test, construction)
{
    EXPECT_TRUE(std::is_nothrow_default_constructible_v<libjst::test::branch_stack_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<libjst::test::branch_stack_t>);
    EXPECT_TRUE(std::is_nothrow_move_constructible_v<libjst::test::branch_stack_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<libjst::test::branch_stack_t>);
    EXPECT_TRUE(std::is_nothrow_move_assignable_v<libjst::test::branch_stack_t>);
    EXPECT_TRUE(std::is_nothrow_destructible_v<libjst::test::branch_stack_t>);
}

TEST(branch_stack_test, associated_types)
{
    EXPECT_TRUE((std::same_as<libjst::test::branch_stack_t::container_type, std::vector<libjst::test::branch>>));
    EXPECT_TRUE((std::same_as<libjst::test::branch_stack_t::value_type, libjst::test::branch>));
    EXPECT_TRUE((std::same_as<libjst::test::branch_stack_t::reference, libjst::test::branch &>));
    EXPECT_TRUE((std::same_as<libjst::test::branch_stack_t::const_reference, libjst::test::branch const &>));
    EXPECT_TRUE((std::same_as<libjst::test::branch_stack_t::size_type, size_t>));
}

TEST(branch_stack_test, push)
{
    libjst::test::branch_stack_t stack{};

    EXPECT_TRUE(stack.empty());

    stack.push(libjst::test::branch{0});
    stack.push(libjst::test::branch{1});
    stack.push(libjst::test::branch{2});

    EXPECT_EQ(stack.size(), 3u);
    EXPECT_EQ(stack.top(), 2);
    EXPECT_EQ(stack.branch_at(0), 0);
    EXPECT_EQ(stack.branch_at(1), 1);
    EXPECT_EQ(stack.branch_at(2), 2);
}

TEST(branch_stack_test, emplace)
{
    libjst::test::branch_stack_t stack{};

    EXPECT_TRUE(stack.empty());

    stack.emplace(0);
    stack.emplace(1);
    stack.emplace(2);

    EXPECT_EQ(stack.size(), 3u);
    EXPECT_EQ(stack.top(), 2);
    EXPECT_EQ(stack.branch_at(0), 0);
    EXPECT_EQ(stack.branch_at(1), 1);
    EXPECT_EQ(stack.branch_at(2), 2);
}

TEST(branch_stack_test, pop)
{
    libjst::test::branch_stack_t stack{};
    stack.emplace(0);
    stack.emplace(1);
    stack.emplace(2);

    EXPECT_EQ(stack.size(), 3u);
    stack.pop();
    EXPECT_EQ(stack.size(), 2u);
    stack.pop();
    EXPECT_EQ(stack.size(), 1u);
    stack.pop();
    EXPECT_EQ(stack.size(), 0u);
}

TEST(branch_stack_test, top)
{
    libjst::test::branch_stack_t stack{};

    EXPECT_TRUE(stack.empty());

    stack.emplace(0);
    EXPECT_EQ(stack.top(), 0);
    stack.emplace(1);
    EXPECT_EQ(stack.top(), 1);
    stack.emplace(2);
    EXPECT_EQ(stack.top(), 2);

    stack.pop();
    EXPECT_EQ(stack.top(), 1);
    stack.pop();
    EXPECT_EQ(stack.top(), 0);
}

TEST(branch_stack_test, empty)
{
    libjst::test::branch_stack_t stack{};

    EXPECT_TRUE(stack.empty());

    stack.emplace(0);
    stack.emplace(1);
    stack.emplace(2);

    EXPECT_FALSE(stack.empty());
}

TEST(branch_stack_test, size)
{
    libjst::test::branch_stack_t stack{};

    EXPECT_EQ(stack.size(), 0u);

    stack.emplace(0);
    stack.emplace(1);
    stack.emplace(2);

    EXPECT_EQ(stack.size(), 3u);
}

TEST(branch_stack_test, branch_at)
{
    libjst::test::branch_stack_t stack{};

    EXPECT_TRUE(stack.empty());

    stack.emplace(0);
    stack.emplace(1);
    stack.emplace(2);

    EXPECT_EQ(stack.branch_at(0), 0);
    EXPECT_EQ(stack.branch_at(1), 1);
    EXPECT_EQ(stack.branch_at(2), 2);
}

TEST(branch_stack_test, base_branch)
{
    libjst::test::branch_stack_t stack{};

    EXPECT_TRUE(stack.empty());

    stack.emplace(0);
    EXPECT_EQ(stack.base_branch(), 0);
    stack.emplace(1);
    EXPECT_EQ(stack.base_branch(), 0);
    stack.emplace(2);
    EXPECT_EQ(stack.base_branch(), 0);
}
