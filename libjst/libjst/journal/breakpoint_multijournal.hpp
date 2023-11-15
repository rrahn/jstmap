// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md===================
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides implementation for the multisequence journal.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <set>

#include <libjst/journal/any_sequence.hpp>
#include <libjst/reference_sequence/reference_sequence_concept.hpp>
#include <libjst/utility/member_type_trait.hpp>

namespace libjst
{
    // TODO: add concept for coverage
    template <libjst::reference_sequence source_t>
    class breakpoint_multijournal
    {
        /// @name Member types
        /// @{
    private:
        class record_impl;

        using breakpoint_map_type = std::multiset<record_impl, std::less<void>>;
    public:
        using source_type = source_t;
        using sequence_type = libjst::breakpoint_slice_t<source_type const>;
        using breakpoint_type = libjst::sequence_breakpoint_t<source_type const>;
        using iterator = std::ranges::iterator_t<breakpoint_map_type const>;
        ///@}

        /// @name Member variables
        /// @{
    private:
        source_type _source{};
        breakpoint_map_type _breakpoint_map{};
        /// @}

        /// @name Member functions
        /// @{
    public:

        constexpr breakpoint_multijournal()
            requires std::default_initializable<source_type> &&
                     std::default_initializable<breakpoint_map_type> = default;

        constexpr breakpoint_multijournal(source_type source)
            noexcept(std::is_nothrow_move_constructible_v<source_type>)
            : _source{std::move(source)}
        {
        }

        constexpr source_type const & source() const & noexcept
        {
            return _source;
        }

        constexpr source_type source() && noexcept
        {
            return std::move(_source);
        }
        /// @}

        /// @name Iterators
        /// @{
    public:

        constexpr iterator begin() const noexcept
        {
            return _breakpoint_map.begin();
        }

        constexpr iterator end() const noexcept
        {
            return _breakpoint_map.end();
        }
        /// @}

        /// @name Modifier
        /// @{
    public:

        template <typename concrete_sequence_t>
            requires std::convertible_to<concrete_sequence_t, sequence_type>
        constexpr iterator record(breakpoint_type breakpoint, concrete_sequence_t && sequence)
        {
            return _breakpoint_map.emplace(std::move(breakpoint), std::forward<concrete_sequence_t>(sequence));
        }
        /// @}

        /// @name Capacities
        /// @{
    public:
        constexpr size_t size() const noexcept
        {
            return _breakpoint_map.size();
        }

        constexpr size_t max_size() const noexcept
        {
            return _breakpoint_map.max_size();
        }

        constexpr bool empty() const noexcept
        {
            return _breakpoint_map.empty();
        }
        /// @}

        /// @name Lookup
        /// @{
        public:

            template <typename breakend_t>
                requires std::convertible_to<breakend_t, std::remove_cvref_t<libjst::low_breakend_t<breakpoint_type>> const &>
            constexpr iterator lower_bound(breakend_t const & breakend) const & noexcept
            {
                using low_breakend_type = std::remove_cvref_t<libjst::low_breakend_t<breakpoint_type>>;
                return _breakpoint_map.lower_bound(static_cast<low_breakend_type const &>(breakend));
            }

            template <typename breakend_t>
                requires std::convertible_to<breakend_t, std::remove_cvref_t<libjst::low_breakend_t<breakpoint_type>> const &>
            constexpr iterator upper_bound(breakend_t const & breakend) const & noexcept
            {
                using low_breakend_type = std::remove_cvref_t<libjst::low_breakend_t<breakpoint_type>>;
                return _breakpoint_map.upper_bound(static_cast<low_breakend_type const &>(breakend));
            }
        /// @}
    };

    template <libjst::reference_sequence source_t>
    class breakpoint_multijournal<source_t>::record_impl
    {
        /// @name Member variables
        /// @{
    private:
        breakpoint_type _breakpoint{};
        any_sequence<sequence_type> _sequence{};
        /// @}

        /// @name Member functions
        /// @{
    public:

        constexpr record_impl() = default;

        template <typename concrete_sequence_t>
        constexpr record_impl(breakpoint_type breakpoint, concrete_sequence_t && sequence)
            noexcept(std::is_nothrow_move_constructible_v<breakpoint_type> &&
                     std::is_nothrow_constructible_v<any_sequence<sequence_type>, concrete_sequence_t>)
            : _breakpoint{std::move(breakpoint)}, _sequence{std::forward<concrete_sequence_t>(sequence)}
        {}

        constexpr sequence_type sequence() const noexcept
        {
            return _sequence.value();
        }
        /// @}

        /// @name Non-member functions
        /// @{
    public:

        template <typename tag_t, typename self_t>
            requires std::same_as<std::remove_cvref_t<self_t>, record_impl> &&
                     std::invocable<tag_t, member_type_t<self_t, breakpoint_type>>
        constexpr friend auto tag_invoke(tag_t const & tag, self_t && self)
            noexcept(std::is_nothrow_invocable_v<tag_t, member_type_t<self_t, breakpoint_type>>)
            -> std::invoke_result_t<tag_t, member_type_t<self_t, breakpoint_type>>
        {
            return tag(std::forward<member_type_t<self_t, breakpoint_type>>(self._breakpoint));
        }

        constexpr friend bool operator==(record_impl const & lhs, record_impl const & rhs) noexcept
        {
            return lhs._breakpoint == rhs._breakpoint && lhs.equivalence_rank() == rhs.equivalence_rank();
        }

        constexpr friend std::weak_ordering operator<=>(record_impl const & lhs, record_impl const & rhs) noexcept
        {
            if (auto cmp = lhs._breakpoint <=> rhs._breakpoint; cmp != 0)
                return cmp;

            return lhs.equivalence_rank() <=> rhs.equivalence_rank();
        }

        constexpr friend std::weak_ordering operator<=>(record_impl const & lhs,
                                                        std::remove_cvref_t<libjst::low_breakend_t<breakpoint_type>> const & breakend) noexcept
        {
            return libjst::low_breakend(lhs._breakpoint) <=> breakend;
        }
        /// @}

        /// @name Utilities
        /// @{
    private:
        constexpr std::ptrdiff_t equivalence_rank() const noexcept
        {
            return libjst::breakend_span(_breakpoint) - std::ranges::size(sequence());
        }
        /// @}
    };
} // namespace libjst
