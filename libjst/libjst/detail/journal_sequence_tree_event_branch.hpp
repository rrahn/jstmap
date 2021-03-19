// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::journal_sequence_tree_event_branch.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <libjst/detail/journal_sequence_tree_event_base.hpp>

namespace libjst::detail
{

/*!\brief Represents a wrapped delta event as a branch event within the libjst::detail::journaled_sequence_tree.
 * \tparam delta_event_t The type of the wrapped delta event.
 *
 * \details
 * \copydetails libjst::detail::journal_sequence_tree_event_base
 */
template <typename delta_event_t>
class journal_sequence_tree_event_branch :
    public journal_sequence_tree_event_base<journal_sequence_tree_event_branch<delta_event_t>, delta_event_t>
{
private:
    //!\brief The base class.
    using base_t = journal_sequence_tree_event_base<journal_sequence_tree_event_branch<delta_event_t>, delta_event_t>;

    friend base_t;
public:
    /*!\name Associated types
     * \{
     */
    using typename base_t::coverage_type; //!\brief The coverage type.
    using typename base_t::size_type; //!\brief The size type.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_event_branch() = default; //!< Default.
    constexpr journal_sequence_tree_event_branch(journal_sequence_tree_event_branch const &) = default; //!< Default.
    constexpr journal_sequence_tree_event_branch(journal_sequence_tree_event_branch &&) = default; //!< Default.
    constexpr journal_sequence_tree_event_branch & operator=(journal_sequence_tree_event_branch const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_event_branch & operator=(journal_sequence_tree_event_branch &&)
        = default; //!< Default.
    ~journal_sequence_tree_event_branch() = default; //!< Default.

    /*!\brief Constructs a new event wrapper for the given delta event.
     * \param[in] delta_event The delta event to wrap.
     */
    explicit constexpr journal_sequence_tree_event_branch(delta_event_t * delta_event) noexcept
        : base_t{delta_event}
    {}
    //!\}

private:
    //!\brief Returns the position of the wrapped delta event.
    size_type position_impl() const noexcept
    {
        assert(base_t::_delta_event != nullptr);

        return base_t::_delta_event->position();
    }
};

}  // namespace libjst::detail
