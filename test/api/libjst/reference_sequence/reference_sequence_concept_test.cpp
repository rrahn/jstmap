// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <catch2/catch_test_macros.hpp>

#include <ranges>
#include <string>
#include <string_view>

#include <libjst/reference_sequence/sequence_breakpoint_simple.hpp>
#include <libjst/reference_sequence/reference_sequence_concept.hpp>

namespace libjst::test
{

    class type_with_tag_invoke_overload
    {
    private:
        std::string _data{"AAAACCCCGGGGTTTT"};

        using iterator_t = std::ranges::iterator_t<std::string>;
        using breakpoint_t = libjst::sequence_breakpoint_simple<std::iter_difference_t<iterator_t>>;
    public:

        constexpr auto begin() noexcept { return _data.begin(); }
        constexpr auto begin() const noexcept { return _data.begin(); }
        constexpr auto end() noexcept { return _data.end(); }
        constexpr auto end() const noexcept { return _data.end(); }
    private:

        constexpr friend auto tag_invoke(libjst::tag_t<libjst::to_breakpoint>,
                                         type_with_tag_invoke_overload const & self,
                                         iterator_t low_it,
                                         iterator_t high_it) noexcept
        {
            return breakpoint_t{.low = std::ranges::distance(std::ranges::begin(self._data), low_it),
                                .high = std::ranges::distance(std::ranges::begin(self._data), high_it)};
        }

        constexpr friend auto tag_invoke(libjst::tag_t<libjst::breakpoint_slice>,
                                         type_with_tag_invoke_overload const & self,
                                         breakpoint_t const & breakpoint) noexcept
        {
            return std::string_view{std::ranges::next(std::ranges::begin(self._data), breakpoint.low),
                                    std::ranges::next(std::ranges::begin(self._data), breakpoint.high)};
        }
    };

    class type_with_member_overload
    {
    private:
        std::string _data{"AAAACCCCGGGGTTTT"};

        using iterator_t = std::ranges::iterator_t<std::string>;
        using breakpoint_t = libjst::sequence_breakpoint_simple<std::iter_difference_t<iterator_t>>;
    public:

            constexpr auto begin() noexcept { return _data.begin(); }
            constexpr auto begin() const noexcept { return _data.begin(); }
            constexpr auto end() noexcept { return _data.end(); }
            constexpr auto end() const noexcept { return _data.end(); }

        constexpr auto to_breakpoint(iterator_t low_it, iterator_t high_it) const noexcept
        {
            return breakpoint_t{.low = std::ranges::distance(std::ranges::begin(_data), low_it),
                                .high = std::ranges::distance(std::ranges::begin(_data), high_it)};
        }

        constexpr auto breakpoint_slice(breakpoint_t const & breakpoint) const noexcept
        {
            return std::string_view{std::ranges::next(std::ranges::begin(_data), breakpoint.low),
                                    std::ranges::next(std::ranges::begin(_data), breakpoint.high)};
        }
    };

    class other_breakpoint
    {

    public:
        constexpr auto low_breakend() const noexcept
        {
            return std::integral_constant<size_t, 0u>{};
        }

        constexpr auto high_breakend() const noexcept
        {
            return std::integral_constant<size_t, 1u>{};
        }

        constexpr auto breakend_span() const noexcept
        {
            return 1;
        }

    private:

        constexpr friend bool operator==(other_breakpoint const &, other_breakpoint const &) = default;
        constexpr friend auto operator<=>(other_breakpoint const &, other_breakpoint const &) = default;
    };


} // namespace libjst::test

SCENARIO("Converting a pair of iterators to a breakpoint", "[reference_sequence][to_breakpoint]")
{
    using simple_breakpoint_t = libjst::sequence_breakpoint_simple<std::iter_difference_t<std::ranges::iterator_t<std::string>>>;
    GIVEN("A type that implements to_breakpoint as member function")
    {
        libjst::test::type_with_member_overload obj{};
        THEN("the breakpoint can be retrieved")
        {
            REQUIRE(std::same_as<libjst::sequence_breakpoint_t<libjst::test::type_with_member_overload>,
                                 simple_breakpoint_t>);
            REQUIRE(libjst::to_breakpoint(obj, std::ranges::begin(obj), std::ranges::end(obj)) == simple_breakpoint_t{0, 16});
        }
    }
    GIVEN("A type that implements to_breakpoint as tag_invoke function")
    {
        libjst::test::type_with_tag_invoke_overload obj{};
        THEN("the breakpoint can be retrieved")
        {
            REQUIRE(std::same_as<libjst::sequence_breakpoint_t<libjst::test::type_with_tag_invoke_overload>,
                                 simple_breakpoint_t>);
            REQUIRE(libjst::to_breakpoint(obj, std::ranges::begin(obj), std::ranges::end(obj)) == simple_breakpoint_t{0, 16});
        }
    }
    GIVEN("A pure std::string")
    {
        std::string str{"AAAACCCCGGGGTTTT"};
        THEN("the breakpoint can be retrieved via the default implementation")
        {
            REQUIRE(std::same_as<libjst::sequence_breakpoint_t<std::string>, simple_breakpoint_t>);
            REQUIRE(libjst::to_breakpoint(str, std::ranges::begin(str), std::ranges::end(str)) == simple_breakpoint_t{0, 16});
        }
    }
}

SCENARIO("Getting the slice of a reference sequence below a given breakpoint", "[reference_sequence][breakpoint_slice]")
{
    std::string_view expected_slice{"AAAACCCCGGGGTTTT"};
    GIVEN("A type that implements breakpoint_slice as member function")
    {
        libjst::test::type_with_member_overload obj{};
        using breakpoint_t = libjst::sequence_breakpoint_t<libjst::test::type_with_member_overload>;
        THEN("the breakpoint slice can be retrieved")
        {
            REQUIRE(std::same_as<libjst::breakpoint_slice_t<libjst::test::type_with_member_overload>, std::string_view>);
            REQUIRE(libjst::breakpoint_slice(obj, breakpoint_t{.low = 0, .high = 16}) == expected_slice);
        }
    }
    GIVEN("A type that implements breakpoint_slice as tag_invoke function")
    {
        libjst::test::type_with_tag_invoke_overload obj{};
        using breakpoint_t = libjst::sequence_breakpoint_t<libjst::test::type_with_member_overload>;
        THEN("the breakpoint slice can be retrieved")
        {
            REQUIRE(std::same_as<libjst::breakpoint_slice_t<libjst::test::type_with_tag_invoke_overload>, std::string_view>);
            REQUIRE(libjst::breakpoint_slice(obj, breakpoint_t{.low = 0, .high = 16}) == expected_slice);
        }
    }
    GIVEN("A pure std::string")
    {
        std::string obj{"AAAACCCCGGGGTTTT"};
        using breakpoint_t = libjst::sequence_breakpoint_t<std::string>;
        THEN("the breakpoint slice can be retrieved via the default implementation")
        {
            REQUIRE(std::same_as<libjst::breakpoint_slice_t<std::string>, std::string_view>);
            auto breakpoint_slice = libjst::breakpoint_slice(obj, breakpoint_t{.low = 0, .high = 16});
            REQUIRE(std::string_view{breakpoint_slice.begin(), breakpoint_slice.end()} == expected_slice);
        }
    }
}

SCENARIO("Testing the reference_sequence concept of different types", "[reference_sequence][concept]")
{
    GIVEN("libjst::test::type_with_tag_invoke_overload")
    {
        libjst::test::type_with_tag_invoke_overload obj{};
        THEN("it models reference_sequence")
        {
            REQUIRE(libjst::reference_sequence<libjst::test::type_with_tag_invoke_overload>);
        }
    }
    GIVEN("libjst::test::type_with_member_overload")
    {
        libjst::test::type_with_member_overload obj{};
        THEN("it models reference_sequence")
        {
            REQUIRE(libjst::reference_sequence<libjst::test::type_with_member_overload>);
        }
    }
    GIVEN("std::string")
    {
        std::string obj{"AAAACCCCGGGGTTTT"};
        THEN("it models reference_sequence")
        {
            REQUIRE(libjst::reference_sequence<std::string>);
        }
    }
}

SCENARIO("Testing the sequence_breakpoint_for concept for different types", "[reference_sequence_for][concept]")
{
    GIVEN("A type of a reference sequence")
    {
        using refseq_t = libjst::test::type_with_tag_invoke_overload;
        WHEN("checking the sequence_breakpoint_for concept with its own breakpoint type")
        {
            using breakpoint_t = libjst::sequence_breakpoint_t<refseq_t>;
            THEN("it models the concept")
            {
                REQUIRE(libjst::sequence_breakpoint<breakpoint_t>);
                REQUIRE(libjst::sequence_breakpoint_for<breakpoint_t, refseq_t>);
            }
        }
        WHEN("checking the sequence_breakpoint_for concept with a convertible breakpoint type")
        {
            using breakpoint_t = libjst::sequence_breakpoint_simple<uint16_t>;
            THEN("it models the concept")
            {
                REQUIRE(libjst::sequence_breakpoint<breakpoint_t>);
                REQUIRE(libjst::sequence_breakpoint_for<breakpoint_t, refseq_t>);
            }
        }
        WHEN("checking the sequence_breakpoint_for concept with an incompatible breakpoint type")
        {
            using breakpoint_t = libjst::test::other_breakpoint;
            THEN ("breakpoint_t is not convertible to the sequence's breakpoint type")
            {
                REQUIRE(!std::convertible_to<breakpoint_t, libjst::sequence_breakpoint_t<refseq_t>>);
            }
            THEN ("it does not model the concept")
            {
                REQUIRE(libjst::sequence_breakpoint<breakpoint_t>);
                REQUIRE(!libjst::sequence_breakpoint_for<breakpoint_t, refseq_t>);
            }
        }
    }
}
