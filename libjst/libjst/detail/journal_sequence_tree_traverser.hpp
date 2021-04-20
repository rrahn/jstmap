// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::detail::journal_sequence_tree_traverser.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <numeric>
#include <ranges>

#include <seqan3/range/views/slice.hpp>
#include <seqan3/range/views/zip.hpp>

#include <libjst/context_position.hpp>
#include <libjst/journal_decorator.hpp>
#include <libjst/detail/branch_stack.hpp>
#include <libjst/utility/logger.hpp>

namespace libjst::detail
{

/*!\brief The basic class to traverse a libjst::journaled_sequence_tree.
 *
 * \tparam derived_t The derived type of this traverser.
 * \tparam jst_t The type of the libjst::journaled_sequence_tree.
 *
 * \details
 *
 * This crtp-base class implements the actual traversal over the libjst::journaled_sequence_tree. It expands each
 * subtree depending on the given context size and the variants contained in the tree. Unsupported branches are not
 * traversed.
 */
template <typename derived_t, typename jst_t>
class journal_sequence_tree_traverser
{
private:

    //!\brief Grant access to the derived class.
    friend derived_t;

    /*!\name Event types
     * \{
     */
    //!\brief The type of the branch event queue.
    using branch_event_queue_type = decltype(std::declval<jst_t>()._branch_event_queue);
    //!\brief The type of the branch event queue iterator.
    using branch_event_queue_iterator = std::ranges::iterator_t<branch_event_queue_type const>;
    //!\brief The type of the branch event.
    using branch_event_type = std::iter_value_t<branch_event_queue_iterator>;

    //!\brief The type of the join event queue.
    using join_event_queue_type = decltype(std::declval<jst_t>()._join_event_queue);
    //!\brief The type of the join event queue iterator.
    using join_event_queue_iterator = std::ranges::iterator_t<join_event_queue_type const>;
    //!\brief The type of the join event.
    using join_event_type = std::iter_value_t<join_event_queue_iterator>;
    //!\}

    /*!\name Delta event types
     * \{
     */
    using delta_event_shared_type = typename jst_t::delta_event_shared_type; //!< The shared delta event type.
    using segment_type = typename delta_event_shared_type::segment_type; //!< The segment type.
    using coverage_type = typename delta_event_shared_type::coverage_type; //!< The coverage type.
    using size_type = typename delta_event_shared_type::size_type; //!< The size type.
    using journal_decorator_type = journal_decorator<segment_type>; //!< The journal decorator type.
    //!\}

    //!\brief The type of the sequence context.
    using sequence_context_type = decltype(std::declval<journal_decorator_type const &>() | seqan3::views::slice(0, 1));

    //!\brief The type to store the positions of a particular sequence context.
    using context_positions_type = std::vector<context_position>;

    struct branch;

    //!\brief Represents the status of the branch creation.
    enum branch_creation_status : uint8_t
    {
        success, //!< A new branch could be created.
        no_support, //!< No new branch was created as it was not supported by any of the sequences.
        success_with_deletion //!< A new branch could be created covering a deletion.
    };

    jst_t const * _jst_host; //!< The referenced jst
    size_t _context_size{}; //!< The size of the underlying context.
    std::vector<int32_t> _sequence_offsets{}; //!< The context position offsets.
    join_event_queue_iterator _join_event_it{}; //!< The global join event iterator to update the context positions.

    branch_stack<branch> _branch_stack{}; //!< The internal stack of branches.

    using context_position_type = context_position; //!< The type representing a single context position.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_traverser() = default; //!< Default.
    constexpr journal_sequence_tree_traverser(journal_sequence_tree_traverser const &) = default; //!< Default.
    constexpr journal_sequence_tree_traverser(journal_sequence_tree_traverser &&) = default; //!< Default.
    constexpr journal_sequence_tree_traverser & operator=(journal_sequence_tree_traverser const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_traverser & operator=(journal_sequence_tree_traverser &&) = default; //!< Default.
    ~journal_sequence_tree_traverser() = default; //!< Default.

    /*!\brief Constructs the journaled sequence tree cursor for a given libjst::journaled_sequence_tree and a context
     *        size.
     *
     * \param[in] jst A pointer to a const journaled sequence tree.
     * \param[in] context_size The context size to use for the cursor.
     *
     * \details
     *
     * The cursor is initialised with the given context size and already represents the first context.
     */
    journal_sequence_tree_traverser(jst_t const * jst, size_t const context_size) noexcept :
        _jst_host{jst},
        _context_size{context_size}
    {
        assert(_jst_host != nullptr);
        assert(_context_size > 0);

        // Initialise the context position offsets and the first join event for updating the context positions.
        _sequence_offsets.resize(_jst_host->size(), 0);
        _join_event_it = std::ranges::begin(join_event_queue());

        // Initialise the base branch covering the reference sequence.
        _branch_stack.emplace(
            size_type{0},  // current context position.
            std::ranges::size(_jst_host->reference()), // current branch end position.
            0, // branch offset.
            nullptr, // pointer to the delta event causing the branch.
            std::ranges::begin(branch_event_queue()), // next branch event to consider for this branch.
            _join_event_it, // next join event to consider for this branch.
            journal_decorator_type{std::span{_jst_host->reference()}}, // the journal decorator of the current branch
            coverage_type(_jst_host->size(), true) // the current branch coverage.
        );

        // Initialise the first branch if any exists at the first position.
        while (on_branch_event())
        {
            if (branch_creation_status state = create_branch(); state == branch_creation_status::success)
            { // The branch could be created, i.e. there is a valid branch at position 0 with at least one sequence
              // covering this branch event.
                assert(!is_base_branch());
                assert(active_branch().coverage.any());
                break; // terminate loop since we found a candidate.
            }
            else if (state == branch_creation_status::success_with_deletion)
            { // The branch could be created but is immediately dropped because the deletion is at the beginning.
              // There is no context spanning over the deletion here. The alternative branch has been updated
              // and points to the next branch event.
                assert(!is_base_branch()); // Cannot be in the base branch.
                drop_branch(); // remove branch and continue to see if there is another branch.
            }
        }
    }
    //!\}

    /*!\name Host member access
     * \{
     */
    //!\brief Returns the event queue containing the branch events from the underlying host.
    branch_event_queue_type const & branch_event_queue() const noexcept
    {
        assert(_jst_host != nullptr);

        return _jst_host->_branch_event_queue;
    }

    //!\brief Returns the event queue containing the join events from the underlying host.
    join_event_queue_type const & join_event_queue() const noexcept
    {
        assert(_jst_host != nullptr);

        return _jst_host->_join_event_queue;
    }

    //!\brief Returns the event queue from the underlying host.
    auto const & reference() const noexcept
    {
        assert(_jst_host != nullptr);

        return _jst_host->reference();
    }
    //!\}

    /*!\name Branch stack operations
     * \{
     */
    //!\brief Returns the active branch.
    branch & active_branch() noexcept
    {
        return _branch_stack.top();
    }

    //!\overload
    branch const & active_branch() const noexcept
    {
        return _branch_stack.top();
    }

    //!\brief Checks if the active branch is the base branch.
    bool is_base_branch() const noexcept
    {
        return _branch_stack.size() == 1u;
    }

    /*!\brief Tests wether the enumerator is at the end.
     *
     * \returns `true`, if no more contexts are available, `false` otherwise.
     */
    bool at_end() const noexcept
    {
        return _branch_stack.empty();
    }

    //!\brief Pushes a new branch on the branch stack.
    void push_branch(branch && new_branch) noexcept
    {
        _branch_stack.push(std::move(new_branch));
        as_derived()->notify_push();
    }

    //!\brief Removes the current branch from the branch stack.
    void drop_branch() noexcept
    {
        assert(!_branch_stack.empty());

        _branch_stack.pop();
        as_derived()->notify_pop();
    }

    //!\brief Makes a new branch at the current position and switches to this branch.
    branch_creation_status create_branch()
    {
        assert(!at_end());
        assert(active_branch().branch_event_it != std::ranges::end(branch_event_queue()));
        assert(on_branch_event());

        // Create a branch from the current one.
        branch new_branch{active_branch()};

        // Update the delta event and the coverage.
        new_branch.delta_event = active_branch().branch_event_it->event_handle();
        active_branch().branch_event_it = next_branch_event(active_branch().branch_event_it);
        update_coverage(new_branch);

        // Terminate early if this branch is not supported by any sequence.
        if (new_branch.coverage.none())
            return branch_creation_status::no_support;

        // Apply the delta event to update the current journal decorator.
        record_delta_event(new_branch);

        new_branch.offset += event_offset(new_branch.delta_event);

        // The max end position will be determined and possibly cropped if the natural end of the
        // journal sequence tree is reached before (e.g. the end of the reference sequence).
        size_type const max_end_position = (is_base_branch())
                                         ? new_branch.delta_event->position() + _context_size +
                                           new_branch.delta_event->insertion_size() - 1
                                         : branch_max_end_position();

        new_branch.branch_end_position = std::min(_branch_stack.base_branch().branch_end_position + new_branch.offset,
                                                  max_end_position);

        new_branch.branch_event_it = find_next_relative_branch_event(new_branch);

        // Only return with no support if this branch is a deletion and is at end of the branch because the deletion
        // is at the end of the reference sequence. But try again if there are more events left, because there could
        // be an insertion at the end of the branch such that it is extended with a valid context.
        if (bool is_deletion = new_branch.delta_event->is_deletion(); is_deletion &&
            (new_branch.at_end() && !has_more_branch_events(new_branch)))
        {
            return branch_creation_status::no_support;
        }
        else
        { // Push the new branch onto the stack, which automatically switches the active branch to the new one.
            push_branch(std::move(new_branch));
            return is_deletion ? branch_creation_status::success_with_deletion : branch_creation_status::success;
        }
    }

    /*!\brief Terminates all branches that have been fully visited.
     *
     * \details
     *
     * After advancing the context the end of a branch might have been reached and the consumed branch will be deleted.
     * The next alternative branch might also have been consumed already, e.g. if the previous branch was extended by an
     * insertion close to the end of the reference sequence. Or the alternative branch has no support by any of the
     * sequences, i.e. the coverage is 0 for each sequence.
     */
    void terminate_consumed_branches() noexcept
    {
        assert(!at_end());

        // Only terminate a branch if it reached the end and either the branch end position is also the branch max
        // end position or there are no more branch events left.
        auto reached_branch_end = [&] () -> bool
        {
            return active_branch().at_end() &&
                   (active_branch().branch_end_position == branch_max_end_position() ||
                   !has_more_branch_events(active_branch()));
        };

        // Drop all visited or unsupported branches except the base branch.
        // Check for 0 coverage only in the very last.
        while (!is_base_branch() && (reached_branch_end() || active_branch().coverage.none()))
            drop_branch();
    }
    //!\}

    /*!\name Branch operations
     * \{
     */

    //!\brief Returns the iterator to the next branch event not the end of the branch queue.
    branch_event_queue_iterator next_branch_event(branch_event_queue_iterator event_it) const noexcept
    {
        return std::ranges::next(event_it, 1, std::ranges::end(branch_event_queue()));
    }

    /*!\brief Finds the next relative valid branch event.
     *
     * \returns The branch iterator to the first branch event that comes relatively after the deleted/replaced region.
     *
     * \details
     *
     * To find the next branch event two cases need to be considered. In the first case the current branch iterator
     * points to a substitution or an deletion. The next iterator will then point to the first branch event that comes
     * after the deleted segment in the reference.
     * In the second case, the branch iterator points to an insertion. In this case, all insertions at the same position
     * are first skipped and then the next valid branch is searched. The rational here is that insertions are joined
     * at the same position where they are spawned, because they "consume" no character from the reference sequence.
     * But one sequence cannot cover two insertions at the same position as this would result in an ambiguous sequence
     * context (depending on which insertion is traversed first). Thus, all insertions at the current position are
     * skipped and then the regular search is conducted.
     */
    branch_event_queue_iterator find_next_relative_branch_event(branch const & new_branch) const noexcept
    {
        auto is_local_insertion = [&] (auto const & event) constexpr
        {
            return event.event_handle()->is_insertion() && event.position() == new_branch.delta_event->position();
        };

        auto it = std::ranges::find_if_not(next_branch_event(new_branch.branch_event_it),
                                           std::ranges::end(branch_event_queue()),
                                           is_local_insertion);

        auto event_lies_behind_deletion = [&] (auto const & event) constexpr
        {
            return event.position() >= new_branch.delta_event->position() + new_branch.delta_event->deletion_size();
        };

        return std::ranges::find_if(it, std::ranges::end(branch_event_queue()), event_lies_behind_deletion);
    }

    //!\brief Tests whether the context position of the active branch lies on a branch event.
    bool on_branch_event() const noexcept
    {
        // Either there are no branches left or the context postion is equal to the postion of the branch event
        // relative to the position of the active branch.
        return has_more_branch_events(active_branch()) &&
               active_branch().context_position == active_branch().relative_branch_event_position();
    }

    /*!\brief Returns the maximal end position of a branch.
     *
     * \returns the maximal end position of the current branch.
     *
     * \details
     *
     * Only one branch will be processed at a time. Thus, the origin of the current branch will be at the second
     * position within the branch stack from the end. The maximal branch position can be computed with the stored
     * branch position as well as the context size and the respective insertion size of the original branch.
     */
    size_type branch_max_end_position() const noexcept
    {
        assert(_branch_stack.size() > 1);

        branch const & branch_origin = _branch_stack.branch_at(1);
        return branch_position() + _context_size + branch_origin.delta_event->insertion_size() - 1;
    }

    size_type branch_position() const noexcept
    {
        assert(_branch_stack.size() > 1);

        return _branch_stack.branch_at(1).delta_event->position();
    }

    /*!\brief Tests wether the full context is available in the active branch.
     *
     * \returns `true` if the full context can be dereferenced, `false` otherwise.
     */
    bool has_full_context_in_branch() const noexcept
    {
        return context_end_position() >= _context_size;
    }

    /*!\brief Tests wether the active branch has more branch events left.
     *
     * \param[in] branch The branch to check if it has more branch events.
     *
     * \returns `true` if the branch iterator is not at the end of the branch queue, `false` otherwise.
     */
    bool has_more_branch_events(branch const & branch) const noexcept
    {
        return branch.branch_event_it != std::ranges::end(branch_event_queue());
    }

    //!\brief Returns the begin position of the context in the active branch.
    size_type context_begin_position() const noexcept
    {
        return active_branch().context_position + 1 - _context_size;
    }

    //!\brief Returns the end position of the context in the active branch.
    size_type context_end_position() const noexcept
    {
        return active_branch().context_position + 1;
    }

    /*!\brief Returns the relative delta event offset.
     *
     * returns An integer value representing the relative offset generated by this delta event.
     *
     * \details
     *
     * Depending on the kind of the delta event the returned value can be negative, 0, or positive.
     * A negative value represents the number of deleted characters with respect to the reference sequence, a zero means
     * that with respect to the reference sequence no character was deleted or added and a positive value accounts for
     * the number of characters inserted relative to the reference sequence.
     */
    auto event_offset(delta_event_shared_type const * delta_event) const noexcept
    {
        using signed_size_t = std::make_signed_t<size_type>;

        return static_cast<signed_size_t>(delta_event->insertion_size()) -
               static_cast<signed_size_t>(delta_event->deletion_size());
    }

    //!\brief Records the represented delta event in the journal decorator of the new branch.
    void record_delta_event(branch & new_branch)
    {
        using substitution_t = delta_event_shared_type::substitution_type;
        using insertion_t = delta_event_shared_type::insertion_type;
        using deletion_t = delta_event_shared_type::deletion_type;

        size_type position = new_branch.delta_event->position() + new_branch.offset;
        journal_decorator_type & jd = new_branch.journal_decorator;
        std::visit([&] (auto & delta_kind)
        {
            seqan3::detail::multi_invocable
            {
                [&] (substitution_t const & e) { jd.record_substitution(position, std::span{e.value()}); },
                [&] (insertion_t const & e) { jd.record_insertion(position, std::span{e.value()}); },
                [&] (deletion_t const & e) { jd.record_deletion(position, position + e.value()); }
            } (delta_kind);
        }, new_branch.delta_event->delta_variant());
    }

    /*!\brief Updates the coverages of the new branch and the parent branch.
     *
     * \details
     *
     * If the new branch comes from the base branch it will be set to the original event coverage, otherwise to the
     * conjunction of parant branch coverage and the branch event coverage. Hence, the new branch is only represented
     * by all sequences that share the current branch and the same branch history.
     * The parent branch coverage is updated accordingly, to represent all sequences without the new branch event.
     */
    void update_coverage(branch & new_branch) noexcept
    {
        coverage_type const & event_coverage = new_branch.delta_event->coverage();
        new_branch.coverage = (is_base_branch()) ? event_coverage : event_coverage & active_branch().coverage;
        active_branch().coverage &= ~(event_coverage);
    }

    /*!\brief Updates the relative context position offset for each sequence.
     *
     * \details
     *
     * When the context advances in any branch, the head position of the respective context might lie behind a join
     * event. In this case, the relative position offsets for all sequences are updated with the respective event
     * offset of the pointed-to delta event, if it is shared by the respective sequence id.
     * The first branch whose context head position fulfills the condition issues the update.
     * This update affects all branches globally and is done only once for every join event.
     */
    void update_relative_sequence_offsets() noexcept
    {
        // Little helper lambda to check if there is really a join event before the current head of the context
        auto is_join_event_before_context_head = [&] ()
        {
            return _join_event_it != std::ranges::end(join_event_queue()) &&
                    context_begin_position() >= _join_event_it->position() &&
                    (is_base_branch() || // Always in base branch.
                     (context_begin_position() <= branch_position() && // In branch only if context position is less &
                      _join_event_it->event_handle() != _branch_stack.branch_at(1).delta_event)); // not the same event.
        };

        for (; has_full_context_in_branch() && is_join_event_before_context_head(); ++_join_event_it)
        {
            coverage_type const & join_coverage = _join_event_it->coverage();
            for (unsigned id = 0; id < _sequence_offsets.size(); ++id)
                _sequence_offsets[id] += (join_coverage[id]) ? event_offset(_join_event_it->event_handle()) : 0;
        }
    }

    /*!\brief Keep the base branch coverage up-to-date with the joined branches.
     *
     * \details
     *
     * When the base branch is advanced and the context head position lies on a join event, the coverage of the base
     * branch will be updated with the coverages of all events pointed-to by the join events on this position.
     */
    void update_base_branch_coverage() noexcept
    {
        auto head_on_join_event = [&] ()
        {
            return active_branch().join_event_it != std::ranges::end(join_event_queue()) &&
                   context_begin_position() == active_branch().join_event_it->position();
        };

        if (is_base_branch() && head_on_join_event())
        {

            for (; head_on_join_event(); ++active_branch().join_event_it)
                active_branch().coverage |= active_branch().join_event_it->coverage();

            // In the base branch the coverage needs to be refined. When joining delta events the coverage will
            // be updated accordingly, which might overwrite a particular bit in the coverage that was previously
            // unset by another delta event that lies inside of the current context of the base branch.

            // The branch_event_it points to the first branch event that lies behind the tail of the current context.
            // Go back to find first branch event that does not affect the coverage of the current context.
            auto first_it = std::ranges::find_if(std::reverse_iterator{active_branch().branch_event_it},
                                                 std::reverse_iterator{std::ranges::begin(branch_event_queue())},
                                                 [&] (auto const & event)
            {
                return (event.event_handle()->is_insertion() && event.position() == context_begin_position()) ||
                       (event.position() < context_begin_position());
            }).base();

            assert(first_it == std::ranges::end(branch_event_queue()) ||
                   first_it->position() >= context_begin_position());

            // Now determine the supported coverage within this context.
            std::ranges::for_each(first_it, active_branch().branch_event_it, [&] (auto const & branch_event)
            {
                active_branch().coverage &= ~branch_event.coverage();
            });
        }
    }

    /*!\brief Advances to the next position.
     *
     * \returns `true` if not at end, `false` otherwise.
     *
     * \details
     *
     * When called the context of the current branch is advanced by one position.
     * After advancing the context the following sequence of operations is performed if applicable:
     * * terminate all consumed branches
     * * update relative context position of each sequence
     * * update base branch coverage with joined branches
     * * create a new branch
     * * terminate base branch if fully consumed
     */
    bool advance()
    {
        assert(!at_end());

        ++active_branch().context_position;
        terminate_consumed_branches(); // terminate all branches
        update_relative_sequence_offsets();
        update_base_branch_coverage();

        // If on a branch event, try to create a new branch.
        while (on_branch_event())
        {
            if (branch_creation_status status = create_branch(); status == branch_creation_status::success)
                break;
            else if (status == branch_creation_status::no_support)
                terminate_consumed_branches();
        }

        // Either on the base branch or there must be at least one supported sequence.
        assert(is_base_branch() || active_branch().coverage.any());

        // Drop the base branch if fully consumed as well.
        if (is_base_branch() && active_branch().at_end() && !has_more_branch_events(active_branch()))
            drop_branch();

        return at_end();
    }

    //!\brief Advances the context by one position.
    //!\returns `true` if there is full context available, otherwise `false`.
    bool next_context()
    {
        return advance() || has_full_context_in_branch();
    }

    /*!\brief Returns the current context.
     *
     * \returns A range representing the current context.
     *
     * \details
     *
     * The returned range is a slice over the generated journal decorator of the current branch.
     */
    sequence_context_type current_context() const noexcept
    {
        return active_branch().journal_decorator
             | seqan3::views::slice(context_begin_position(), context_end_position());
    }

    //!\brief Returns the current value pointed to by the traverser.
    //!\returns The current pointed-to value.
    std::ranges::range_value_t<journal_decorator_type> current_value() const noexcept
    {
        assert(active_branch().context_position < active_branch().journal_decorator.size());
        return *(active_branch().journal_decorator.begin() + active_branch().context_position);
    }

   /*!\brief Compute the branch coverage that is valid for the current branch.
     *
     * \returns The valid branch coverage for the current context.
     *
     * \details
     *
     * The traversal does not keep a valid coverage all the time. When asking for the positions the iterator
     * first refines the current coverage to produce a supported coverage, effectively deferring these operations until
     * they are requested by the caller.
     */
    coverage_type determine_supported_context_coverage() noexcept
    {
        // Just use coverage in these cases.
        if (branch_event_queue().empty() || is_base_branch() || context_begin_position() >= branch_position())
            return active_branch().coverage;

        // On non-base branches the coverage will not be updated when joining delta events.
        // Instead, this update is done here only when requested.
        // To do this, the interval of join events which lie between the head of the context and the original
        // branch event is found, and then the supported coverage is determined.
        auto join_begin = std::ranges::find_if(active_branch().join_event_it,
                                               std::ranges::end(join_event_queue()),
                                               [&] (auto const & join_event)
        {
            return join_event.position() > context_begin_position();
        });

        // Find the end of the interval.
        auto join_end = std::ranges::find_if(join_begin,
                                             std::ranges::end(join_event_queue()),
                                             [&] (auto const & join_event)
        {
            return (join_event.event_handle() == _branch_stack.branch_at(1).delta_event) ||
                    join_event.position() > branch_position();
        });

        // Update the join event so in case we need to call context position on the same branch again, we can
        // start from a closer event.
        active_branch().join_event_it = join_begin;

        // Determine the supported coverage.
        coverage_type unsupported = std::accumulate(join_begin,
                                                    join_end,
                                                    coverage_type(_sequence_offsets.size(), false),
                                                    [] (coverage_type coverage, auto const & event)
        {
            return coverage | event.coverage();
        });

        return active_branch().coverage & ~unsupported;
    }
    //!\}

    //!\brief Cast this to its derived type.
    derived_t * as_derived() noexcept
    {
        return static_cast<derived_t *>(this);
    }

    //!\overload
    derived_t const * as_derived() const noexcept
    {
        return static_cast<derived_t const *>(this);
    }
};

//!\brief The type of a traversal branch.
template <typename derived_t, typename jst_t>
struct journal_sequence_tree_traverser<derived_t, jst_t>::branch
{
    size_type context_position{}; //<! The current tail position of the moving window.
    size_type branch_end_position{}; //!<! The end position of the branch.
    std::make_signed_t<size_type> offset{}; //<! The offset generated by the current branch.
    delta_event_shared_type const * delta_event{}; //!< The pointer to the current delta event.
    branch_event_queue_iterator branch_event_it{}; //!< The iterator pointing to the next branch event.
    join_event_queue_iterator join_event_it{}; //!< The iterator pointing to the next join event.
    journal_decorator_type journal_decorator{}; //!< The journal decorator representing the current sequence context.
    coverage_type coverage{}; //!< The coverage for this branch.

    //!\brief Checks wether this branch reached its end.
    constexpr bool at_end() const noexcept
    {
        return context_position == branch_end_position;
    }

    //!\brief Returns the position of the delta event.
    constexpr size_type delta_event_position() const noexcept
    {
        return delta_event->position();
    }

    //!\brief Returns the position of the pointed-to branch event.
    constexpr size_type branch_event_position() const noexcept
    {
        return branch_event_it->position();
    }

    //!\brief Returns the position of the pointed-to join event.
    constexpr size_type join_event_position() const noexcept
    {
        return join_event_it->position();
    }

    //!\brief Returns the relative branch event position.
    constexpr size_type relative_branch_event_position() const noexcept
    {
        return branch_event_position() + offset;
    }

};
}  // namespace libjst::detail
