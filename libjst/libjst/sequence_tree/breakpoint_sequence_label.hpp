// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides basic breakpoint sequence label implementation.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <libjst/reference_sequence/sequence_breakpoint_concept.hpp>
#include <libjst/reference_sequence/sequence_concept.hpp>
#include <libjst/utility/member_type_trait.hpp>
#include <libjst/utility/tag_invoke.hpp>

namespace libjst
{

template <sequence sequence_t, sequence_breakpoint breakpoint_t> class breakpoint_sequence_label
{
    /// @name Member types
    /// @{
  private:
    using low_breakend_type = libjst::low_breakend_t<breakpoint_t>;
    using high_breakend_type = libjst::high_breakend_t<breakpoint_t>;
    /// @}

    /// @name Member variables
    /// @{
  private:
    sequence_t _sequence{};
    breakpoint_t _breakpoint{};
    /// @}

    /// @name Member functions
    /// @{
  public:
    constexpr breakpoint_sequence_label() = default; //!< Default.

    explicit constexpr breakpoint_sequence_label(sequence_t sequence, breakpoint_t breakpoint) noexcept(
        std::is_nothrow_move_constructible_v<breakpoint_t> &&std::is_nothrow_move_constructible_v<sequence_t>)
        : _sequence{std::move(sequence)}, _breakpoint{std::move(breakpoint)}
    {
    }

    explicit constexpr breakpoint_sequence_label(
        sequence_t sequence, low_breakend_type low_breakend,
        high_breakend_type high_breakend) noexcept(std::is_nothrow_move_constructible_v<low_breakend_type>
                                                       &&std::is_nothrow_move_constructible_v<high_breakend_type>
                                                           &&std::is_nothrow_move_constructible_v<sequence_t>)
        : breakpoint_sequence_label{std::move(sequence),
                                    breakpoint_t{std::move(low_breakend), std::move(high_breakend)}}
    {
    }

    constexpr sequence_t const &sequence() const & noexcept
    {
        return _sequence;
    }

    constexpr sequence_t &&sequence() && noexcept
    {
        return std::move(_sequence);
    }
    /// @}

    /// @name Non-member functions
    /// @{
  public:
    template <typename tag_t, typename self_t>
        requires std::same_as<std::remove_cvref_t<self_t>, breakpoint_sequence_label> &&
                 std::invocable<tag_t, member_type_t<self_t, breakpoint_t>>
    constexpr friend auto tag_invoke(tag_t const &tag, self_t &&self) noexcept(
        std::is_nothrow_invocable_v<tag_t, member_type_t<self_t, breakpoint_t>>)
        -> std::invoke_result_t<tag_t, member_type_t<self_t, breakpoint_t>>
    {
        return tag((member_type_t<self_t, breakpoint_t> &&)self._breakpoint);
    }
    /// @}
};
} // namespace libjst
