// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides concept for sequence variants.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <libjst/reference_sequence/sequence_breakpoint_concept.hpp>
#include <libjst/reference_sequence/sequence_breakpoint_simple.hpp>
#include <libjst/reference_sequence/sequence_concept.hpp>
#include <libjst/utility/tag_invoke.hpp>

namespace libjst
{

    namespace _to_breakpoint
    {
        inline constexpr struct _tag
        {

        public:
            template <libjst::sequence sequence_t>
                requires (tag_or_member_invocable<_tag, sequence_t>)
            constexpr auto operator()(sequence_t &&sequence, std::ranges::iterator_t<sequence_t> low, std::ranges::iterator_t<sequence_t> high) const
                noexcept(libjst::nothrow_tag_or_member_invocable<_tag, sequence_t>)
                    -> tag_or_member_invoke_result_t<_tag, sequence_t>
            {
                return libjst::tag_or_member_invoke(_tag{}, (sequence_t &&)sequence, std::move(low), std::move(high));
            }

            // Possible default for standard flat sequences would be to return simple breakpoints with positions
            template <libjst::sequence sequence_t>
                requires (!tag_or_member_invocable<_tag, sequence_t>)
            constexpr auto operator()(sequence_t &&sequence, std::ranges::iterator_t<sequence_t> low, std::ranges::iterator_t<sequence_t> high) const
            {
                if constexpr (std::ranges::random_access_range<sequence_t>) {
                    auto min_high = std::max(low, high);
                    return libjst::sequence_breakpoint_simple{
                        .low = std::ranges::distance(std::ranges::begin(sequence), std::move(low)),
                        .high = std::ranges::distance(std::ranges::begin(sequence), std::move(min_high))
                    };
                }
            }
        private:

            template <typename sequence_t>
                requires requires (sequence_t &&sequence, std::ranges::iterator_t<sequence_t> l, std::ranges::iterator_t<sequence_t> h) {
                    { std::forward<sequence_t>(sequence).to_breakpoint(l, h) } -> libjst::sequence_breakpoint;
                }
            constexpr friend auto tag_invoke(tag_t<libjst::member_invoke>, _tag const &, sequence_t && sequence, std::ranges::iterator_t<sequence_t> low, std::ranges::iterator_t<sequence_t> last)
                noexcept(noexcept(std::forward<sequence_t>(sequence).to_breakpoint(low, last)))
                -> decltype(std::forward<sequence_t>(sequence).to_breakpoint(low, last))
            {
                return std::forward<sequence_t>(sequence).to_breakpoint(low, last);
            }
        } to_breakpoint;
    }
    using _to_breakpoint::to_breakpoint;

    template <typename sequence_t>
    using sequence_breakpoint_t = std::invoke_result_t<_to_breakpoint::_tag, sequence_t, std::ranges::iterator_t<sequence_t>, std::ranges::iterator_t<sequence_t>>;

    namespace _breakpoint_slice
    {
        inline constexpr struct _tag
        {
        private:

            template <typename breakpoint_t>
            inline static constexpr bool has_integral_breakends =
                std::integral<std::remove_reference_t<libjst::low_breakend_t<breakpoint_t>>> &&
                std::integral<std::remove_reference_t<libjst::high_breakend_t<breakpoint_t>>>;

        public:
            template <libjst::sequence sequence_t, libjst::sequence_breakpoint breakpoint_t>
                requires (tag_or_member_invocable<_tag, sequence_t, breakpoint_t>)
            constexpr auto operator()(sequence_t &&sequence, breakpoint_t &&breakpoint) const
                noexcept(libjst::nothrow_tag_or_member_invocable<_tag, sequence_t, breakpoint_t>)
                    -> tag_or_member_invoke_result_t<_tag, sequence_t, breakpoint_t>
            {
                return libjst::tag_or_member_invoke(_tag{}, (sequence_t &&)sequence, (breakpoint_t &&)breakpoint);
            }

            template <libjst::sequence sequence_t, libjst::sequence_breakpoint breakpoint_t>
                requires (!tag_or_member_invocable<_tag, sequence_t, breakpoint_t>)
            constexpr auto operator()(sequence_t &&sequence, breakpoint_t &&breakpoint) const
            {
                if constexpr (has_integral_breakends<breakpoint_t>)
                {
                    auto segment_begin = std::ranges::begin(sequence) + libjst::low_breakend((breakpoint_t &&)breakpoint);
                    auto segment_end = std::ranges::begin(sequence) + libjst::high_breakend((breakpoint_t &&)breakpoint);
                    if constexpr (std::ranges::contiguous_range<sequence_t>) {
                        return std::span{std::move(segment_begin), std::move(segment_end)};
                    } else {
                        return std::ranges::subrange{std::move(segment_begin), std::move(segment_end)};
                    }
                }
            }

        private:

            template <typename sequence_t, typename breakpoint_t>
                requires requires (sequence_t &&sequence, breakpoint_t &&breakpoint) {
                    { std::forward<sequence_t>(sequence).breakpoint_slice(std::forward<breakpoint_t>(breakpoint)) };// -> libjst::sequence;
                }
            constexpr friend auto tag_invoke(tag_t<libjst::member_invoke>, _tag const &, sequence_t && sequence, breakpoint_t && breakpoint)
                noexcept(noexcept(std::forward<sequence_t>(sequence).breakpoint_slice(std::forward<breakpoint_t>(breakpoint))))
                -> decltype(std::forward<sequence_t>(sequence).breakpoint_slice(std::forward<breakpoint_t>(breakpoint)))
            {
                return std::forward<sequence_t>(sequence).breakpoint_slice(std::forward<breakpoint_t>(breakpoint));
            }
        } breakpoint_slice;
    }
    using _breakpoint_slice::breakpoint_slice;

    template <typename sequence_t>
    using breakpoint_slice_t = std::invoke_result_t<_breakpoint_slice::_tag, sequence_t, libjst::sequence_breakpoint_t<sequence_t>>;

    // ----------------------------------------------------------------------------
    // Concept defintions
    // ----------------------------------------------------------------------------

    template <typename object_t>
    concept reference_sequence = libjst::sequence<object_t> &&
        requires(std::remove_reference_t<object_t> const & obj)
    {
        typename libjst::sequence_breakpoint_t<object_t>;
        typename libjst::breakpoint_slice_t<object_t>;

        { libjst::to_breakpoint(obj, std::ranges::begin(obj), std::ranges::end(obj)) } -> libjst::sequence_breakpoint;
        { libjst::breakpoint_slice(obj, libjst::to_breakpoint(obj, std::ranges::begin(obj), std::ranges::end(obj))) } -> libjst::sequence;
    };

    template <typename breakpoint_t, typename sequence_t>
    concept sequence_breakpoint_for =
        libjst::sequence_breakpoint<breakpoint_t> &&
        libjst::reference_sequence<sequence_t> &&
        requires (breakpoint_t && breakpoint, sequence_t && sequence)
    {
        { libjst::breakpoint_slice(std::forward<sequence_t>(sequence), std::forward<breakpoint_t>(breakpoint)) } -> libjst::sequence;
    };

    template <typename object_t>
    concept preserving_reference_sequence =
        libjst::reference_sequence<object_t> &&
        libjst::reference_sequence<libjst::breakpoint_slice_t<object_t>> &&
        std::convertible_to<libjst::breakpoint_slice_t<libjst::breakpoint_slice_t<object_t>>, libjst::breakpoint_slice_t<object_t>> &&
        std::convertible_to<libjst::sequence_breakpoint_t<libjst::breakpoint_slice_t<object_t>>, libjst::sequence_breakpoint_t<object_t>>;

} // namespace libjst
