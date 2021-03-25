// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::journal_sequence_tree_event_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

namespace libjst::detail
{

/*!\brief An abstract crtp base class for the jst events.
 * \tparam derived_t The type of the derived class.
 * \tparam delta_event_t The type of the wrapped delta event.
 *
 * \details
 *
 * Wraps a delta event and overrides the position interface based on the given derived class.
 * This wrapper class is used by the libjst::journaled_sequence_tree to keep track of the branch events and join events
 * in the tree, without storing the same delta event multiple times. The differentiation between join and branch events
 * is required by the traversal over the journal sequence tree.
 */
template <typename derived_t, typename delta_event_t>
class journal_sequence_tree_event_base
{
private:
    using coverage_type = typename delta_event_t::coverage_type; //!< The coverage type.
    using size_type = typename delta_event_t::size_type; //!< The size type.

    friend derived_t;

    delta_event_t * _delta_event{}; //!< The wrapped delta event.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_event_base() = default; //!< Default.
    constexpr journal_sequence_tree_event_base(journal_sequence_tree_event_base const &) = default; //!< Default.
    constexpr journal_sequence_tree_event_base(journal_sequence_tree_event_base &&) = default; //!< Default.
    constexpr journal_sequence_tree_event_base & operator=(journal_sequence_tree_event_base const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_event_base & operator=(journal_sequence_tree_event_base &&)
        = default; //!< Default.
    ~journal_sequence_tree_event_base() = default; //!< Default.

    /*!\brief Constructs a new event wrapper for the given delta event.
     * \param[in] delta_event The delta event to wrap.
     */
    explicit constexpr journal_sequence_tree_event_base(delta_event_t * delta_event) noexcept :
        _delta_event{delta_event}
    {}
    //!\}

public:
    /*!\name Element access
     * \{
     */
    //!\brief Returns the coverage of the wrapped delta event.
    constexpr coverage_type const & coverage() const noexcept
    {
        assert(_delta_event != nullptr);
        return _delta_event->coverage();
    }

    /*!\brief Returns the position of the wrapped delta event.
     *
     * \returns The position of the event.
     *
     * \details
     *
     * The branch event position corresponds to the position of the wrapped delta event, while the join event position
     * corresponds to the position of the wrapped delta event plus the deletion size of this event.
     */
    constexpr size_type position() const noexcept
    {
        return as_derived().position_impl();
    }

    //!\brief Returns a pointer to the wrapped event.
    constexpr delta_event_t * event_handle() noexcept
    {
        return _delta_event;
    }

    //!\overload
    constexpr delta_event_t const * event_handle() const noexcept
    {
        return _delta_event;
    }
    //!\}

    /*!\name Comparison
     * \{
     */
    //!\brief Compares the wrapped events for equality.
    constexpr bool operator==(journal_sequence_tree_event_base const & rhs) const noexcept
    {
        return event_handle() == rhs.event_handle();
    }

    /*!\brief Compares the wrapped events by their position.
     *
     * \param[in] rhs The right hand side of the comparison.
     *
     * \returns std::weak_ordering
     *
     * \details
     *
     * Compares the position of `this` with the position of `rhs`. If both positions are equivalent, the order depends
     * on the delta kind. In the branch event the order is `insertion < substitution < deletion` and in the join event
     * the order is `deletion < substitution < insertion`.
     */
    constexpr std::weak_ordering operator<=>(journal_sequence_tree_event_base const & rhs) const noexcept
    {
        if (std::weak_ordering ordering = position() <=> rhs.position(); ordering == std::weak_ordering::equivalent)
            return delta_index() <=> rhs.delta_index();
        else
            return ordering;
    }

    //!\brief Compares the wrapped event with another position.
    constexpr std::weak_ordering operator<=>(size_type const & rhs) const noexcept
    {
        return position() <=> rhs;
    }
    //!\}

private:
    /*!\brief Returns the possibly modified index of the delta variant of the wrapped delta event.
     *
     * \details
     *
     * The table shows the index returned by this method depending on the delta kind and the event type:
     *
     * |   delta kind  | branch event | join event |
     * |---------------|--------------|------------|
     * | insertion     |      0       |      2     |
     * | substitution  |      1       |      1     |
     * | deletion      |      2       |      0     |
     */
    constexpr size_t delta_index() const noexcept
    {
        return as_derived().delta_index_impl();
    }

    //!\brief Returns this casted to the derived type.
    derived_t & as_derived()
    {
        return static_cast<derived_t &>(*this);
    }

    //!\brief overload
    derived_t const & as_derived() const
    {
        return static_cast<derived_t const &>(*this);
    }
};

}  // namespace libjst::detail
