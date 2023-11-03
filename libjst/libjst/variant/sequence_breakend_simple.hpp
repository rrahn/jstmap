// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides breakend.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst
{
    class sequence_breakend_simple
    {

        ///@name Member types
        ///@{
    public:
        using value_type = uint32_t;
        ///@}

        ///@name Constructors, destructor and assignment
        ///@{
    public:
        constexpr sequence_breakend_simple() noexcept = default;
        constexpr sequence_breakend_simple(value_type value) noexcept : _value{value}
        {
        }
        ///@}

        ///@name Member variables
        ///@{
    private:
        value_type _value{};
        ///@}

        ///@name Member functions
        ///@{
    public:
        constexpr value_type value() const noexcept
        {
            return _value;
        }
        ///@}

        ///@name Comparison operators
        ///@{
    public:
        friend constexpr bool operator==(sequence_breakend_simple const &, sequence_breakend_simple const &) noexcept = default;
        friend constexpr auto operator<=>(sequence_breakend_simple const &, sequence_breakend const &) noexcept = default;
        ///@}

        ///@name Difference operator
        ///@{
    public:
        friend constexpr std::ptrdiff_t operator-(sequence_breakend_simple const & lhs, sequence_breakend_simple const & rhs) noexcept
        {
            return lhs.value() - rhs.value();
        }
        ///@}
    };
} // namespace libjst
