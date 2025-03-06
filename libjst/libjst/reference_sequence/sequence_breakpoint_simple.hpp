// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides sequence breakend span.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

namespace libjst
{
    /**
     * @brief A simple sequence breakpoint using integral types to define the breakends.
     *
     * @tparam breakend_t The type of the breakend to use. Must model `std::integral`.
     */
    template <std::integral breakend_t>
    struct sequence_breakpoint_simple
    {
        ///@name Member variables
        ///@{
        breakend_t low{};
        breakend_t high{};
        ///@}

        ///@name Member functions
        ///@{
        constexpr breakend_t & low_breakend() & noexcept { return low; }
        constexpr breakend_t const & low_breakend() const & noexcept { return low; }
        constexpr breakend_t && low_breakend() && noexcept { return std::move(low); }
        constexpr breakend_t const && low_breakend() const && noexcept { return std::move(low); }

        constexpr breakend_t & high_breakend() & noexcept { return high; }
        constexpr breakend_t const & high_breakend() const & noexcept { return high; }
        constexpr breakend_t && high_breakend() && noexcept { return std::move(high); }
        constexpr breakend_t const && high_breakend() const && noexcept { return std::move(high); }
        ///@}

        ///@name Conversion operators
        ///@{
        template <typename convertible_breakend_t>
            requires (!std::same_as<convertible_breakend_t, breakend_t>) &&
                     (std::convertible_to<convertible_breakend_t, breakend_t>)
        constexpr operator sequence_breakpoint_simple<convertible_breakend_t>() const noexcept
        {
            return {low, high};
        }
        ///@}

        ///@name Comparison operators
        ///@{
    public:
        friend constexpr bool operator==(sequence_breakpoint_simple const &, sequence_breakpoint_simple const &) noexcept = default;
        friend constexpr auto operator<=>(sequence_breakpoint_simple const &, sequence_breakpoint_simple const &) noexcept = default;
        ///@}
    };
} // namespace libjst
