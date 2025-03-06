// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>

#include <libjst/reference_sequence/sequence_breakpoint_concept.hpp>

namespace libjst::test
{
    using tag_invoke_breakend = uint32_t;

    class breakpoint_with_tag_invoke
    {
        constexpr friend tag_invoke_breakend tag_invoke(libjst::tag_t<libjst::low_breakend>, breakpoint_with_tag_invoke const &)
        {
            return 0u;
        }

        constexpr friend tag_invoke_breakend tag_invoke(libjst::tag_t<libjst::high_breakend>, breakpoint_with_tag_invoke const &)
        {
            return 1u;
        }

        constexpr friend std::ptrdiff_t tag_invoke(libjst::tag_t<libjst::breakend_span>, breakpoint_with_tag_invoke const &)
        {
            return 1;
        }

    private:

        constexpr friend bool operator==(breakpoint_with_tag_invoke const &, breakpoint_with_tag_invoke const &) = default;
        constexpr friend auto operator<=>(breakpoint_with_tag_invoke const &, breakpoint_with_tag_invoke const &) = default;
    };

    using member_breakend = uint32_t;

    class breakpoint_with_member
    {

    public:
        constexpr auto low_breakend() const noexcept
        {
            return member_breakend{0u};
        }

        constexpr auto high_breakend() const noexcept
        {
            return member_breakend{1u};
        }

        constexpr auto breakend_span() const noexcept
        {
            return 1;
        }

    private:

        constexpr friend bool operator==(breakpoint_with_member const &, breakpoint_with_member const &) = default;
        constexpr friend auto operator<=>(breakpoint_with_member const &, breakpoint_with_member const &) = default;
    };

    struct subtractable_breakend
    {

        int32_t value{};

        private:

        constexpr int32_t friend operator-(subtractable_breakend const & lhs, subtractable_breakend const & rhs)
        {
            return lhs.value - rhs.value;
        }

        constexpr bool friend operator==(subtractable_breakend const &, subtractable_breakend const &) = default;
        constexpr auto friend operator<=>(subtractable_breakend const &, subtractable_breakend const &) = default;

    };

    class breakpoint_no_breakend_span_member
    {
        public:

        constexpr auto low_breakend() const noexcept
        {
            return subtractable_breakend{0u};
        }

        constexpr auto high_breakend() const noexcept
        {
            return subtractable_breakend{1u};
        }

        constexpr auto breakend_span() = delete;

    private:

        constexpr friend bool operator==(breakpoint_no_breakend_span_member const &, breakpoint_no_breakend_span_member const &) = default;
        constexpr friend auto operator<=>(breakpoint_no_breakend_span_member const &, breakpoint_no_breakend_span_member const &) = default;
    };


} // namespace libjst::test

SCENARIO("Getting the low breakend of an breakpoint", "[reference_sequence][sequence_breakpoint][low_breakend]")
{
    GIVEN("A type that implements low_breakend as member function")
    {
        libjst::test::breakpoint_with_member obj{};
        THEN("the low breakend can be retrieved")
        {
            REQUIRE(std::same_as<libjst::low_breakend_t<libjst::test::breakpoint_with_member>, libjst::test::member_breakend>);
            REQUIRE(libjst::low_breakend(obj) == libjst::test::member_breakend{0u});
        }
    }
    GIVEN("A type that implements low_breakend as tag_invoke function")
    {
        libjst::test::breakpoint_with_tag_invoke obj{};
        THEN("the low breakend can be retrieved")
        {
            REQUIRE(std::same_as<libjst::low_breakend_t<libjst::test::breakpoint_with_tag_invoke>, libjst::test::tag_invoke_breakend>);
            REQUIRE(libjst::low_breakend(obj) == libjst::test::tag_invoke_breakend{0u});
        }
    }
}

SCENARIO("Getting the high breakend of a breakpoint", "[reference_sequence][sequence_breakpoint][high_breakend]")
{
    GIVEN("A type that implements high_breakend as member function")
    {
        libjst::test::breakpoint_with_member obj{};
        THEN("the high breakend can be retrieved")
        {
            REQUIRE(std::same_as<libjst::high_breakend_t<libjst::test::breakpoint_with_member>, libjst::test::member_breakend>);
            REQUIRE(libjst::high_breakend(obj) == libjst::test::member_breakend{1u});
        }
    }
    GIVEN("A type that implements high_breakend as tag_invoke function")
    {
        libjst::test::breakpoint_with_tag_invoke obj{};
        THEN("the high breakend can be retrieved")
        {
            REQUIRE(std::same_as<libjst::high_breakend_t<libjst::test::breakpoint_with_tag_invoke>, libjst::test::tag_invoke_breakend>);
            REQUIRE(libjst::high_breakend(obj) == libjst::test::tag_invoke_breakend{1u});
        }
    }
}

SCENARIO("Getting the span of a breakpoint", "[reference_sequence][sequence_breakpoint][breakend_span]")
{
    GIVEN("A type that implements breakend_span as tag_invoke function")
    {
        libjst::test::breakpoint_with_tag_invoke breakpoint{};
        THEN("the span of the breakends can be retrieved")
        {
            REQUIRE(std::integral<libjst::breakend_span_t<libjst::test::breakpoint_with_tag_invoke>>);
            REQUIRE(libjst::breakend_span(breakpoint) == 1);
        }
    }
    GIVEN("A type that implements breakend_span as member function")
    {
        libjst::test::breakpoint_with_member breakpoint{};
        THEN("the span of the breakends can be retrieved")
        {
            REQUIRE(std::integral<libjst::breakend_span_t<libjst::test::breakpoint_with_member>>);
            REQUIRE(libjst::breakend_span(breakpoint) == 1);
        }
    }
    GIVEN("A type whose breakends are subtractable and no overload for breakend_span exists")
    {
        libjst::test::breakpoint_no_breakend_span_member breakpoint{};
        THEN("the span of the breakpoints can be retrieved")
        {
            REQUIRE(std::integral<libjst::breakend_span_t<libjst::test::breakpoint_no_breakend_span_member>>);
            REQUIRE(libjst::breakend_span(breakpoint) == 1);
        }
    }
}

SCENARIO("Testing the sequence_breakpoint concept of different types", "[reference_sequence][sequence_breakpoint][concept]")
{
    GIVEN("libjst::test::breakpoint_with_member")
    {
        THEN("the sequence_breakpoint concept evaluates to true")
        {
            REQUIRE(libjst::sequence_breakpoint<libjst::test::breakpoint_with_member>);
        }
    }
    GIVEN("libjst::test::breakpoint_with_tag_invoke")
    {
        THEN("the sequence_breakpoint concept evaluates to true")
        {
            REQUIRE(libjst::sequence_breakpoint<libjst::test::breakpoint_with_tag_invoke>);
        }
    }
    GIVEN("libjst::test::breakpoint_no_breakend_span_member")
    {
        THEN("the sequence_breakpoint concept evaluates to true")
        {
            REQUIRE(libjst::sequence_breakpoint<libjst::test::breakpoint_no_breakend_span_member>);
        }
    }
    GIVEN("non-breakpoint types")
    {
        THEN("the sequence_breakpoint concept evaluates to false")
        {
            REQUIRE(!libjst::sequence_breakpoint<int>);
            REQUIRE(!libjst::sequence_breakpoint<libjst::test::subtractable_breakend>);
        }
    }
}
