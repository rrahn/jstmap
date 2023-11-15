// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the
// 3-clause BSD-License shipped with this file and also available at:
// https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/**
 * @file
 * @brief Provides breakpoint tree node.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <concepts>
#include <iterator>
#include <optional>
#include <ranges>

#include <libjst/reference_sequence/reference_sequence_concept.hpp>
#include <libjst/reference_sequence/sequence_breakpoint_concept.hpp>
#include <libjst/sequence_tree/breakpoint_sequence_label.hpp>
#include <libjst/sequence_tree/breakpoint_sequence_tree_sentinel.hpp>

namespace libjst
{

template <typename journal_t> class breakpoint_sequence_tree_node
{
    /// @name Member types
    /// @{
  private:
    using journal_iterator_type = std::ranges::iterator_t<journal_t const &>;
    using record_type = std::iter_reference_t<journal_iterator_type>;

    using low_breakend_type = std::remove_cvref_t<low_breakend_t<record_type>>;
    using high_breakend_type = std::remove_cvref_t<high_breakend_t<record_type>>;
    using breakpoint_type = std::pair<low_breakend_type, high_breakend_type>;
    static_assert(std::same_as<low_breakend_type, high_breakend_type>,
                  "At the moment it is expected that low and high breakends have "
                  "the same type.");

    using sequence_type =
        typename journal_t::sequence_type; // TODO: maybe change to
                                           // alt_sequence_t<std::iter_reference_t<journal_iterator_type>>
    using label_type = breakpoint_sequence_label<sequence_type, breakpoint_type>;
    /// @}

    /// @name Member variables
    /// @{
  private:
    journal_t const *_journal{};
    journal_iterator_type _prev_breakpoint{};
    journal_iterator_type _next_breakpoint{};
    std::optional<label_type> _label{};
    /// @}

    /// @name Member functions
    /// @{
  public:
    constexpr breakpoint_sequence_tree_node() = default; //!< Default.

    explicit constexpr breakpoint_sequence_tree_node(journal_t const &journal) noexcept : _journal{&journal}
    {
        // now initialize the first sequence.
        _next_breakpoint = _journal->begin();
        auto [src_low, src_high] = libjst::to_breakpoint(_journal->source(), std::ranges::begin(_journal->source()),
                                                         std::ranges::end(_journal->source()));
        auto breakpoint = (_next_breakpoint != _journal->end())
                              ? breakpoint_type{std::move(src_low), libjst::low_breakend(*_next_breakpoint)}
                              : breakpoint_type{std::move(src_low), std::move(src_high)};
        auto ref_slice = libjst::breakpoint_slice(_journal->source(), breakpoint);

        _label = label_type{std::move(ref_slice), libjst::low_breakend(std::move(breakpoint)),
                            libjst::high_breakend(std::move(breakpoint))};
    }

  private:
    explicit constexpr breakpoint_sequence_tree_node(
        journal_t const &journal, journal_iterator_type prev_breakpoint, journal_iterator_type next_breakpoint,
        std::optional<label_type> label) noexcept(std::is_nothrow_move_constructible_v<label_type>
                                                      &&std::is_nothrow_move_constructible_v<journal_iterator_type>)
        : _journal{&journal}, _prev_breakpoint{std::move(prev_breakpoint)},
          _next_breakpoint{std::move(next_breakpoint)}, _label{std::move(label)}
    {
    }
    /// @}

    /// @name Observers
    /// @{
  public:
    constexpr label_type const &value() const &
    {
        return _label.value();
    }

    constexpr label_type &&value() &&
    {
        return std::move(_label).value();
    }

    constexpr bool is_nil() const noexcept
    {
        return !_label.has_value();
    }
    /// @}

    /// @name Node operations
    /// @{
    constexpr breakpoint_sequence_tree_node next_ref() const noexcept
    {
        if (_next_breakpoint == _journal->end())
            return breakpoint_sequence_tree_node{*_journal, _next_breakpoint, _next_breakpoint, std::nullopt};

        journal_iterator_type child_prev_breakpoint = _next_breakpoint;
        journal_iterator_type child_next_breakpoint = std::ranges::next(child_prev_breakpoint);

        if (is_overlapping(child_prev_breakpoint, child_next_breakpoint))
        {
            child_next_breakpoint = _journal->lower_bound(libjst::high_breakend(*child_prev_breakpoint));
        }

        assert(child_prev_breakpoint != _journal->end());
        auto child_low_breakend = libjst::high_breakend(value());

        auto child_high_breakend =
            (child_next_breakpoint != _journal->end()) ? libjst::low_breakend(*child_next_breakpoint) : max_breakend();
        assert(child_high_breakend - child_low_breakend >= 0);

        auto ref_slice =
            libjst::breakpoint_slice(_journal->source(), std::tie(child_low_breakend, child_high_breakend));

        label_type child_label{std::move(ref_slice), std::move(child_low_breakend), std::move(child_high_breakend)};

        return breakpoint_sequence_tree_node{*_journal, std::move(child_prev_breakpoint),
                                             std::move(child_next_breakpoint), std::move(child_label)};
    }

    constexpr breakpoint_sequence_tree_node next_alt() const noexcept
    {
        if (is_alt_node() || _next_breakpoint == _journal->end())
            return breakpoint_sequence_tree_node{*_journal, _prev_breakpoint, _next_breakpoint, std::nullopt};

        label_type child_label{(*_next_breakpoint).sequence(), libjst::low_breakend(*_next_breakpoint),
                               libjst::high_breakend(*_next_breakpoint)};
        return breakpoint_sequence_tree_node{*_journal, _next_breakpoint, _next_breakpoint, std::move(child_label)};
    }
    /// @}

    /// @name Utilities
    /// @{
  private:
    constexpr bool is_alt_node() const noexcept
    {
        return _prev_breakpoint == _next_breakpoint;
    }

    constexpr bool is_overlapping(journal_iterator_type last_breakpoint,
                                  journal_iterator_type next_breakpoint) const noexcept
    {
        return next_breakpoint != _journal->end() && is_alt_node() &&
               libjst::low_breakend(*next_breakpoint) < libjst::high_breakend(*last_breakpoint);
    }

    constexpr auto max_breakend() const noexcept
    {
        return libjst::high_breakend(libjst::to_breakpoint(_journal->source(), std::ranges::begin(_journal->source()),
                                                           std::ranges::end(_journal->source())));
    }
    /// @}

    /// @name Non-member functions
    /// @{
  public:
    constexpr friend bool operator==(breakpoint_sequence_tree_node const &lhs,
                                     breakpoint_sequence_tree_sentinel const &) noexcept
    {
        return lhs._prev_breakpoint == lhs._journal->end();
    }
    /// @}
};
} // namespace libjst
