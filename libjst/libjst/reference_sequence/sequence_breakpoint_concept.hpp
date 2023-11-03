// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides concept for sequence breakpoint.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <utility>

#include <libjst/utility/tag_invoke.hpp>

namespace libjst
{

    namespace exposition_only
    {
        template <typename t>
        concept pair_like = requires (t && obj)
        {
            typename std::tuple_element_t<0, std::remove_cvref_t<t>>;
            typename std::tuple_element_t<1, std::remove_cvref_t<t>>;

            requires std::tuple_size_v<std::remove_cvref_t<t>> == 2;

            { std::get<0>(std::forward<t>(obj)) } -> std::convertible_to<std::tuple_element_t<0, std::remove_cvref_t<t>>>;
            { std::get<1>(std::forward<t>(obj)) } -> std::convertible_to<std::tuple_element_t<1, std::remove_cvref_t<t>>>;
        };
    }

    namespace _low_breakend
    {
        inline constexpr struct _tag
        {

        public:
            template <typename breakpoint_t>
                requires (libjst::tag_or_member_invocable<_tag, breakpoint_t const &>)
            constexpr auto operator()(breakpoint_t const &breakpoint) const
                noexcept(libjst::nothrow_tag_or_member_invocable<_tag, breakpoint_t const &>)
                    -> libjst::tag_or_member_invoke_result_t<_tag, breakpoint_t const &>
            {
                return libjst::tag_or_member_invoke(_tag{}, breakpoint);
            }

            template <typename breakpoint_t>
                requires (!libjst::tag_or_member_invocable<_tag, breakpoint_t const &>)
            constexpr auto operator()(breakpoint_t const &breakpoint) const
            {
                if constexpr (exposition_only::pair_like<breakpoint_t>)
                    return std::get<0>(breakpoint);
            }

        private:

            template <typename breakpoint_t>
                requires requires(breakpoint_t const &breakpoint) { { breakpoint.low_breakend() }; }
            constexpr friend auto tag_invoke(tag_t<libjst::get_tag_member>, _tag const &, breakpoint_t const &) noexcept
            {
                return std::mem_fn(&breakpoint_t::low_breakend);
            }

        } low_breakend;
    }
    using _low_breakend::low_breakend;

    template <typename breakpoint_t>
        requires std::invocable<_low_breakend::_tag, breakpoint_t>
    using low_breakend_t = std::invoke_result_t<_low_breakend::_tag, breakpoint_t>;

    namespace _high_breakend
    {
        inline constexpr struct _tag
        {

        public:
            template <typename breakpoint_t>
                requires (libjst::tag_or_member_invocable<_tag, breakpoint_t const &>)
            constexpr auto operator()(breakpoint_t const &breakpoint) const
                noexcept(libjst::nothrow_tag_or_member_invocable<_tag, breakpoint_t const &>)
                    -> libjst::tag_or_member_invoke_result_t<_tag, breakpoint_t const &>
            {
                return libjst::tag_or_member_invoke(_tag{}, breakpoint);
            }

            template <typename breakpoint_t>
                requires (!libjst::tag_or_member_invocable<_tag, breakpoint_t const &>)
            constexpr auto operator()(breakpoint_t const &breakpoint) const
            {
                if constexpr (exposition_only::pair_like<breakpoint_t>)
                    return std::get<1>(breakpoint);
            }

        private:

            template <typename breakpoint_t>
                requires requires(breakpoint_t const &breakpoint) { { breakpoint.high_breakend() }; }
            constexpr friend auto tag_invoke(tag_t<libjst::get_tag_member>, _tag const &, breakpoint_t const &) noexcept
            {
                return std::mem_fn(&breakpoint_t::high_breakend);
            }
        } high_breakend;
    }
    using _high_breakend::high_breakend;

    template <typename breakpoint_t>
        requires std::invocable<_high_breakend::_tag, breakpoint_t>
    using high_breakend_t = std::invoke_result_t<_high_breakend::_tag, breakpoint_t>;

    namespace _breakend_span
    {
        template <typename t>
        concept breakends_subtractable =
            requires (libjst::high_breakend_t<t> const & lhs, libjst::low_breakend_t<t> const & rhs) {
            { lhs - rhs } -> std::integral;
        };

        inline constexpr struct _tag
        {

        public:
            template <typename breakpoint_t>
                requires libjst::tag_or_member_invocable<_tag, breakpoint_t const &>
            constexpr auto operator()(breakpoint_t const &breakpoint) const
                noexcept(libjst::nothrow_tag_or_member_invocable<_tag, breakpoint_t const &>)
                    -> libjst::tag_or_member_invoke_result_t<_tag, breakpoint_t const &>
            {
                return libjst::tag_or_member_invoke(_tag{}, breakpoint);
            }

            template <typename breakpoint_t>
                requires (!libjst::tag_or_member_invocable<_tag, breakpoint_t const &>)
            constexpr decltype(auto) operator()(breakpoint_t const &breakpoint) const
            {
                if constexpr (breakends_subtractable<breakpoint_t>) {
                    return libjst::high_breakend(breakpoint) - libjst::low_breakend(breakpoint);
                }
            }

        private:

            template <typename breakpoint_t>
                requires requires(breakpoint_t const &breakpoint) { { breakpoint.breakend_span() }; }
            constexpr friend auto tag_invoke(tag_t<libjst::get_tag_member>, _tag const &, breakpoint_t const &) noexcept
            {
                return std::mem_fn(&breakpoint_t::breakend_span);
            }
        } breakend_span;
    }
    using _breakend_span::breakend_span;

    template <typename breakpoint_t>
        requires std::invocable<_breakend_span::_tag, breakpoint_t>
    using breakend_span_t = std::invoke_result_t<_breakend_span::_tag, breakpoint_t>;

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    /**
     * @concept
     * @brief Concept describing the constraints of a sequence breakpoint.
     * @tparam object_t The object type to check.
     *
     * A sequence breakpoint is a type that models `std::totally_ordered`.
     * Given an instance `object_t const & object` the following interfaces are valid:
     * @rst
     *   * ``libjst::low_breakend(object)``
     *   * ``libjst::high_breakend(object)``
     *   * ``libjst::breakend_span(object)``
     * @endrst
     * and where the breakend span returns a std::integral value.
     */
    template <typename object_t>
    concept sequence_breakpoint =
        std::totally_ordered<object_t> &&
        requires(std::remove_reference_t<object_t> const &object) {
            { libjst::low_breakend(object) };
            { libjst::high_breakend(object) };
            { libjst::breakend_span(object) } -> std::integral;
        };

} // namespace libjst
