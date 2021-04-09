// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journaled_sequence_tree.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <list>
#include <set>

#include <seqan3/std/algorithm>
#include <seqan3/std/filesystem>
#include <iostream>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

#include <cereal/types/vector.hpp>
#include <cereal/types/list.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/convert.hpp>

#include <libjst/detail/delta_event_shared.hpp>
#include <libjst/detail/journal_sequence_tree_event_branch.hpp>
#include <libjst/detail/journal_sequence_tree_event_join.hpp>
#include <libjst/detail/transform_to_delta_events.hpp>
#include <libjst/journal_decorator.hpp>
#include <libjst/journal_sequence_tree_context_enumerator.hpp>

namespace libjst::no_adl
{
//!\brief Specific class implementation in no_adl namespace to avoid ADL of template arguments.
template <seqan3::sequence sequence_t>
class journaled_sequence_tree
{
public:
    //!\brief The inner type to break ADL lookup-cycle of template parameters.
    class type;
};

//!\brief Implements the actual journaled sequence tree type.
template <seqan3::sequence sequence_t>
class journaled_sequence_tree<sequence_t>::type
{
private:
    using alphabet_t = std::ranges::range_value_t<sequence_t>; //!< The alphabet type of the underlying sequence.

    using delta_event_shared_type = detail::delta_event_shared<alphabet_t>; //!< The shared delta event type.
    //!\brief The type used to mark a delta event type as a branch event.
    using branch_event_type = detail::journal_sequence_tree_event_branch<delta_event_shared_type>;
    //!\brief The type used to mark a delta event type as a join event.
    using join_event_type = detail::journal_sequence_tree_event_join<delta_event_shared_type>;
    //!\brief The coverage type.
    using coverage_type = typename delta_event_shared_type::coverage_type;
    //!\brief The container type that stores all delta events.
    using event_list_type = std::list<delta_event_shared_type>;

    using segment_type = typename delta_event_shared_type::segment_type; //!< The segment type.
    using journal_decorator_type = journal_decorator<segment_type>; //!< The journal decorator type.

    //!\cond
    template <typename>
    friend class detail::journal_sequence_tree_context_enumerator;
    //!\endcond

    sequence_t _reference; //!< The internal reference used for referential compression.
    event_list_type _delta_events{}; //!< The list of stored delta events.
    std::multiset<branch_event_type, std::less<void>> _branch_event_queue{}; //!< The queue of branch events.
    std::multiset<join_event_type, std::less<void>> _join_event_queue{}; //!< The queue of join events.
    size_t _size{}; //!< The sequence count.

public:
    /*!\name Associated types
     * \{
     */
    using sequence_type = sequence_t; //!< The type of the underlying sequence.
    using size_type = typename delta_event_shared_type::size_type; //!< The size type.
    //!\brief The type of the context enumerator.
    using context_enumerator_type = detail::journal_sequence_tree_context_enumerator<type>;
    using event_type = delta_event_shared_type; //!< The internally stored event type.
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr type() = default; //!< Default.
    constexpr type(type const &) = default; //!< Default.
    constexpr type(type &&) = default; //!< Default.
    constexpr type & operator=(type const &) = default; //!< Default.
    constexpr type & operator=(type &&) = default; //!< Default.
    ~type() = default; //!< Default.

    /*!\brief Constructs the journaled sequence tree with a given reference sequence.
     *
     * \param[in] reference The reference sequence to set.
     * \param[in] count The number of represented sequences. Defaults to `0`.
     *
     * \details
     *
     * The journaled sequence tree takes ownership over the passed reference sequence. Accordingly, only temporaries
     * or moved sequences can be used. To access the reference later one can use the
     * libjst::journaled_sequence_tree::reference member function to access the stored reference.
     * Thus, one can use this reference sequence as long as the journaled sequence tree is valid.
     *
     * If count is given, then the constructor initialises the libjst::journaled_sequence_tree to contain `count`
     * sequences. All of these sequences represent the original reference sequence. The sequence variation can be
     * added later by the libjst::no_adl::journaled_sequence_tree::type::insert_event interface.
     */
    type(sequence_type && reference, size_type const count = 0) : _reference{std::move(reference)}, _size(count)
    {}
    //!\}

    //!\brief Returns the stored reference sequence.
    sequence_type const & reference() const
    {
        return _reference;
    }

    /*!\brief Returns the target sequence at the specified index.
     *
     * \param index The index to get the target sequence from.
     *
     * \returns An immutable target sequence at the specified index.
     *
     * \details
     *
     * This function reconstructs the original target sequence. It generates a libjst::journal_decorator from the
     * shared delta events and the reference sequence.
     *
     * ### Complexity
     *
     * Linear in the number of delta events.
     *
     * ### Exception
     *
     * Strong exception guarantee.
     * Throws std::out_of_range exception if the index is out of range.
     */
    journal_decorator_type sequence_at(size_type const index) const
    {
        using namespace std::literals;

        if (index >= size())
            throw std::out_of_range{"The index "s + std::to_string(index) +
                                    " is out of range [0, "s + std::to_string(size()) + ")"s};

        journal_decorator_type target_sequence{reference()};
        int32_t target_position_offset = 0;

        // Go over the branch event list since it is sorted by the reference position.
        std::ranges::for_each(_branch_event_queue, [&] (branch_event_type const & branch_event)
        {
            delta_event_shared_type const * delta_event = branch_event.event_handle();

            if (!delta_event->coverage()[index]) // Continue if event is not covered by target sequence.
                return;

            // Record the event in the journal decorator.
            std::visit([&] (auto const & event_kind)
            {
                using substitution_t = typename delta_event_shared_type::substitution_type;
                using insertion_t = typename delta_event_shared_type::insertion_type;
                using deletion_t = typename delta_event_shared_type::deletion_type;

                size_type target_position = delta_event->position() + target_position_offset;

                seqan3::detail::multi_invocable
                {
                    [&] (substitution_t const & e)
                    {
                        target_sequence.record_substitution(target_position, segment_type{e.value()});
                    },
                    [&] (insertion_t const & e)
                    {
                        target_sequence.record_insertion(target_position, segment_type{e.value()});
                    },
                    [&] (deletion_t const & e)
                    {
                        target_sequence.record_deletion(target_position, target_position + e.value());
                    }
                }(event_kind);
            }, delta_event->delta_variant());

            // Update the target position offset by the insertion and deletion size.
            target_position_offset += static_cast<int32_t>(delta_event->insertion_size() -
                                                           delta_event->deletion_size());
        });

        return target_sequence;
    }

    //!\brief Returns the number of stored sequences.
    size_type size() const
    {
        return _size;
    }

    /*!\brief Inserts a new event to the existing journal sequence tree.
     *
     * \param[in] event The event to insert.
     *
     * \returns `true` if the event could be emplaced, `false` otherwise.
     *
     * \details
     *
     * Only inserts the event
     *  * if the event position is less than `size()` (less than or equal to `size()` in case of an insertion),
     *  * if the join poisition is less than `size()` (less than or equal to `size()` in case of an insertion),
     *  * if no other sequence has an overlapping event,
     *  * the event coverage has no bit set.
     *
     * ### Exception
     *
     * Throws std::length_error if the size of the event coverage is not equal to `size()`.
     * If any exception is thrown, this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in the number of events whose join position is after the branch position of the inserted element.
     */
    bool insert(event_type event)
    {
        if (event.coverage().size() != size())
            throw std::length_error{"The coverage length: " + std::to_string(event.coverage().size()) +
                                    " differs from the actual size: " + std::to_string(size()) + "!"};

        size_type const event_join_position = event.position() + event.deletion_size();
        size_type const max_size = std::ranges::size(reference()) + event.is_insertion();

        if (event.position() >= max_size || event_join_position >= max_size || event.coverage().none())
            return false;

        // Find the first event whose join position is not less than the event position of the insert element.
        auto it = _join_event_queue.lower_bound(event.position());

        for (; it != _join_event_queue.end(); ++it)
        {
            // In general, if the join position of the other event is less than or equal to the begin position of the
            // insert event it can be ignored. The same is true for events whose begin position is greater or equal
            // than the join position of the insert event.
            // The only exception is when both events are insertions at the same position. Then their coverage always
            // needs to be compared. This is checked with the add_one constant.
            event_type const & other_event = *(it->event_handle());
            size_t const add_one = other_event.is_insertion() && event.is_insertion();
            if (((it->position() + add_one) <= event.position()) ||
                (other_event.position() >= (event_join_position + add_one)))
                continue;

            coverage_type shared_coverage = other_event.coverage() & event.coverage();
            if (shared_coverage.any())
                return false;
        }

        add_new_event(std::move(event));
        return true;
    }

    /*!\brief Inserts a new event to the existing journal sequence tree.
     *
     * \tparam args_t A template parameter pack; must model std::constructible_from with
     *                libjst::no_adl::journaled_sequence_tree::type::event_type.
     *
     * \param[in] args The parameter pack to construct the
     *
     * \returns `true` if the event could be inserted, `false` otherwise.
     *
     * \details
     *
     * Only inserts the event
     *  * if the event position is less than `size()` (less than or equal to `size()` in case of an insertion),
     *  * if the join poisition is less than `size()` (less than or equal to `size()` in case of an insertion),
     *  * if no other sequence has an overlapping event,
     *  * the event coverage has no bit set.
     *
     * ### Exception
     *
     * Throws std::length_error if the size of the event coverage is not equal to `size()`.
     * If any exception is thrown, this function has no effect (strong exception guarantee).
     *
     * ### Complexity
     *
     * Linear in the number of events whose end position is after the inserted position.
     */
    template <typename ...args_t>
    //!\cond
        requires std::constructible_from<event_type, args_t...>
    //!\endcond
    bool emplace(args_t && ...args)
    {
        return insert(event_type(std::forward<args_t>(args)...));
    }

    /*!\brief Adds a new sequence to the journaled sequence tree based on the given pairwise alignment.
     *
     * \tparam alignment_t The type of the alignment to add; must model pairwise_alignment.
     *
     * \param[in] alignment The alignment to add to this object.
     *
     * \details
     *
     * This method adds a sequence to the current journaled sequence tree using referential compression. The given
     * alignment is used to determine the transcript to derive the added sequence from the stored reference sequence.
     * Only pairwise alignments are supported.
     * The first sequence of the alignment must be identical to the sequence returned by
     * libjst::journaled_sequence_tree::reference after all gap characters have been removed.
     * The second sequence is added encoded by the given alignment.
     */
    template <typename alignment_t>
    void add(alignment_t && alignment)
    {
        // Step 1: just store full sequences:
        using alphabet_t = std::ranges::range_value_t<sequence_t>;
        using gapped_alphabet_t = seqan3::gapped<alphabet_t>;

        auto & [ref, target] = alignment;

        constexpr auto out_gaps = [] (gapped_alphabet_t const c) -> bool { return c != seqan3::gap{}; };

        if (!std::ranges::equal(ref | std::views::filter(out_gaps), _reference))
        {
            throw std::invalid_argument{"The first aligned sequence must be equal to the reference sequence of this "
                                        "journaled sequence tree without the gaps."};
        }

        // Step 1: Increase coverages by one.
        size_type new_size = size() + 1;
        std::ranges::for_each(_delta_events, [new_size] (delta_event_shared_type & delta_event_shared)
        {
            delta_event_shared.coverage().resize(new_size, false);
        });

        // Step 2: extract the delta and update the event queue and delta events.
        auto delta_events = detail::transform_to_delta_events<alphabet_t>(alignment);
        for (auto event : delta_events)
        {
            auto [first_it, last_it] = _branch_event_queue.equal_range(event.position()); // how to compare?

            if (first_it == last_it) // empty range -> no key found.
            {
                // handle the event.
                add_new_event(std::move(event));
            }
            else // non-empty range
            {
                auto event_it = std::ranges::find_if(first_it, last_it, [&] (branch_event_type const & branch_event)
                {
                    return event == *branch_event.event_handle();
                });

                if (event_it == last_it) // event type is not present
                    add_new_event(std::move(event));
                else // event exists.
                    update_event(event_it);
            }
        }
        ++_size;
        assert(validate_added_sequence_with(size() - 1, target)); // check if the added sequence is valid.
    }

    //!\brief Returns a new context enumerator over the current journaled sequence tree.
    context_enumerator_type context_enumerator(size_t const context_size) const noexcept
    {
        // TODO: Throw if not valid context.
        return detail::journal_sequence_tree_context_enumerator<type>{this, context_size};
    }

    /*!\brief Saves this journaled sequence tree to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in, out] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(_reference, _delta_events, _size);
    }

    /*!\brief Loads this journaled sequence tree from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model seqan3::cereal_input_archive.
     *
     * \param[in, out] archive The archive to serialise this object from.
     */
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(_reference, _delta_events, _size);

        // Refill the queues which just store pointers.
        std::ranges::for_each(_delta_events, [&] (auto & event)
        {
            _branch_event_queue.emplace(std::addressof(event));
            _join_event_queue.emplace(std::addressof(event));
        });
    }

    //!\cond
    // [REMOVE_ME]
    void print_event_queue() const
    {
        auto is_branch_event = [] <typename event_t> (event_t) -> bool
        {
            return std::same_as<event_t, detail::journal_sequence_tree_event_branch<delta_event_shared_type>>;
        };

        auto print_event = [&] (auto const & event)
        {
            std::cout << "["
                      << (is_branch_event(event) ? "b" : "j")
                      << ": "
                      << (*event.event_handle())
                      << "]\n";
        };

        std::ranges::for_each(_branch_event_queue, print_event);
        std::ranges::for_each(_join_event_queue, print_event);
    }
    //!\endcond

private:
    /*!\brief Adds a new event to the journal sequence tree.
     *
     * \param[in] delta_event The delta event to add.
     *
     * \details
     *
     * If a new delta event was detected while converting the aligned target sequence, a new shared delta event will
     * be constructed and added to the collection of delta events. Also the corresponding branch and respectively
     * join event queue will be updated by placing the events into the sorted queue based on their branch,
     * respectively join position.
     *
     * ### Complexity
     *
     * Appending the shared delta event is done in constant time. Inserting the branch and join events depends on
     * complexity of inserting elements in the underlying data structure.
     *
     * ### Exception
     *
     * Basic exception guarantee. If an exception is thrown by any operation, this function has no effect.
     */
    void add_new_event(detail::delta_event<alphabet_t> && delta_event)
    {
        coverage_type new_coverage{};
        new_coverage.resize(size() + 1, false);
        new_coverage.back() = true;

        add_new_event(event_type{std::move(delta_event), std::move(new_coverage)});
    }

    //!\overload
    void add_new_event(event_type && shared_event)
    {
        _delta_events.push_back(std::move(shared_event));
        delta_event_shared_type * event_handle = std::addressof(_delta_events.back());

        auto branch_event_it = _branch_event_queue.emplace(event_handle);
        try
        {
            _join_event_queue.emplace(event_handle);
        }
        catch (...)
        { // Remove the inserted branch event such that event and join queue are in balance again.
          // Note the delta events must not be updated because no event references the inserted element anymore.
          // If 8 bytes cannot be added anymore to the queues, then we need to worry about other things.
            _branch_event_queue.erase(branch_event_it);

            assert(_branch_event_queue.size() == _join_event_queue.size());
            std::rethrow_exception(std::current_exception());
        }
    }

    /*!\brief Updates the coverage of an existing event.
     *
     * \param[in] branch_event_it An iterator to the branch event to be updated.
     *
     * \details
     *
     * If an event already exists, the corresponding coverage will be updated by flipping the last bit.
     */
    template <typename branch_event_queue_iterator_t>
    void update_event(branch_event_queue_iterator_t branch_event_it) noexcept
    {
        assert(branch_event_it->coverage().size() == size() + 1);

        delta_event_shared_type * event_ptr = const_cast<delta_event_shared_type *>(branch_event_it->event_handle());
        event_ptr->coverage().back() = true;
    }

    /*!\brief Validates that the transformation of the given alignment was correct.
     *
     * \param[in] idx The index of the sequence to validate.
     * \param[in] target The aligned target sequence from the input.
     *
     * \returns `true` if the added events can reconstruct the aligned target sequence, otherwise `false`.
     *
     * \details
     *
     * Constructs a journal decorator from the recorded delta events that have a bit at the 'idx' position of the
     * associated coverage and compares this generated sequence with the aligned target sequence without gaps.
     * Note this function is only used during debug builds to validate that the alignment has been correctly
     * transformed into a set of delta events.
     */
    template <typename aligned_sequence_t>
    bool validate_added_sequence_with(size_t const idx, aligned_sequence_t && target) const
    {
        using alphabet_t = std::ranges::range_value_t<sequence_t>;
        using gapped_alphabet_t = seqan3::gapped<alphabet_t>;
        sequence_t pure_target_sequence;
        std::ranges::copy(target | std::views::filter([] (gapped_alphabet_t const c) -> bool { return c != seqan3::gap{}; })
                                 | std::views::transform([](gapped_alphabet_t const c)
                                   {
                                        return c.template convert_to<alphabet_t>();
                                   }),
                          std::cpp20::back_inserter(pure_target_sequence));

        // The target sequence and the journal decorator must be equal.
        return std::ranges::equal(sequence_at(idx), pure_target_sequence);
    }
};
} // namespace libjst::no_adl

namespace libjst
{

/*!\brief A referentially compressed sequence tree over collection of sequences.
 *
 * \tparam sequence_t The type of the sequences to store; must model seqan3::sequence.
 *
 * \details
 *
 * This class stores a collection of sequences in a referentially compressed way to tremendously reduce the
 * memory footprint for storing large collections of sequences with a high similarity. Sequences can be added by
 * adding an alignment between the stored reference sequence and the respective target sequence. This class further
 * supports a special enumerator to enable an efficient, compression parallel traversal over the compressed sequences.
 * This can be used in conjunction with any context based streaming algorithm to speed-up the search against large
 * collection of sequences.
 */
template <seqan3::sequence sequence_t>
using journaled_sequence_tree = typename no_adl::journaled_sequence_tree<sequence_t>::type;
}  // namespace libjst
