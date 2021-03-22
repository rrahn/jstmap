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

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/convert.hpp>

#include <libjst/detail/delta_event_shared.hpp>
#include <libjst/detail/journal_sequence_tree_event_branch.hpp>
#include <libjst/detail/journal_sequence_tree_event_join.hpp>
#include <libjst/detail/transform_to_delta_events.hpp>
#include <libjst/journaled_sequence_tree_cursor.hpp>
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
    using sequence_collection_type = std::vector<sequence_t>; //!< The type to store the sequence collection. [REMOVE_ME]

    using delta_event_shared_type = detail::delta_event_shared<alphabet_t>; //!< The shared delta event type.
    //!\brief The type used to mark a delta event type as a branch event.
    using branch_event_type = detail::journal_sequence_tree_event_branch<delta_event_shared_type>;
    //!\brief The type used to mark a delta event type as a join event.
    using join_event_type = detail::journal_sequence_tree_event_join<delta_event_shared_type>;
    //!\brief The coverage type.
    using coverage_type = typename delta_event_shared_type::coverage_type;
    //!\brief The container type that stores all delta events.
    using event_list_type = std::list<delta_event_shared_type>;

    //!\brief Befriend the cursor type. [REMOVE_ME]
    template <typename>
    friend class journaled_sequence_tree_cursor;

    //!\cond
    template <typename>
    friend class detail::journal_sequence_tree_context_enumerator;
    //!\endcond

    sequence_t _reference; //!< The internal reference used for referential compression.
    sequence_collection_type _sequences; //!< The stored sequences. [REMOVE_ME]

    event_list_type _delta_events{}; //!< The list of stored delta events.
    std::multiset<branch_event_type, std::less<void>> _branch_event_queue{}; //!< The queue of branch events.
    std::multiset<join_event_type, std::less<join_event_type>> _join_event_queue{}; //!< The queue of join events.

public:
    /*!\name Associated types
     * \{
     */
    using sequence_type = sequence_t; //!< The type of the underlying sequence.
    using size_type = typename delta_event_shared_type::size_type; //!< The size type.
    //!\brief The type of the context enumerator.
    using context_enumerator_type = detail::journal_sequence_tree_context_enumerator<type>;
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
     *
     * \details
     *
     * The journaled sequence tree takes ownership over the passed reference sequence. Accordingly, only temporaries
     * or moved sequences can be used. To access the reference later one can use the
     * libjst::journaled_sequence_tree::reference member function to access the stored reference.
     * Thus, one can use this reference sequence as long as the journaled sequence tree is valid.
     */
    type(sequence_type && reference) : _reference{std::move(reference)}
    {}
    //!\}

    //!\brief Returns the stored reference sequence.
    sequence_type const & reference() const
    {
        return _reference;
    }

    //!\brief Returns the number of stored sequences.
    size_type size() const
    {
        return _sequences.size();
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

        // TODO: Remove old code.
        sequence_t pure_target_sequence;
        std::ranges::copy(target | std::views::filter(out_gaps)
                                 | std::views::transform([](gapped_alphabet_t const c)
                                   {
                                        return c.template convert_to<alphabet_t>();
                                   }),
                          std::cpp20::back_inserter(pure_target_sequence));

        _sequences.push_back(std::move(pure_target_sequence));

        assert(validate_added_sequence(size() - 1)); // check if the added sequence is valid.
    }

    //!\brief Returns a new context enumerator over the current journaled sequence tree.
    auto context_enumerator(size_t const context_size) const noexcept
    {
        return libjst::journaled_sequence_tree_cursor<type>{this, context_size};
    }

    //!\brief Returns a new context enumerator over the current journaled sequence tree.
    context_enumerator_type context_enumerator_2(size_t const context_size) const noexcept
    {
        // TODO: Throw if not valid context.
        return detail::journal_sequence_tree_context_enumerator<type>{this, context_size};
    }

    /*!\brief Saves the journaled sequence tree to the given output archive.
     *
     * \tparam output_archive_t The type of the output_archive; must model seqan3::cereal_output_archive.
     *
     * \param[in] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_output_archive output_archive_t>
    void save(output_archive_t & archive) const
    {
        archive(_reference, _sequences);
    }

    /*!\brief Loads the journaled sequence tree from the given input archive.
     *
     * \tparam input_archive_t The type of the input_archive; must model seqan3::cereal_input_archive.
     *
     * \param[in] archive The archive to serialise this object to.
     */
    template <seqan3::cereal_input_archive input_archive_t>
    void load(input_archive_t & archive)
    {
        archive(_reference, _sequences);
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

        _delta_events.emplace_back(std::move(delta_event), std::move(new_coverage));
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
     *
     * \returns `true` if the added events can reconstruct the aligned target sequence, otherwise `false`.
     *
     * \details
     *
     * Constructs a journal decorator from the recorded delta events that have a bit at the 'idx' position of the
     * associated coverage and compares this generated sequence with the aligned target sequence without gaps.
     */
    bool validate_added_sequence(size_t const idx) const
    {
        // Create journal decorator and register all events.
        libjst::journal_decorator jd{std::span{reference()}};
        int32_t target_position_offset = 0;

        // We need to go over the list but it is not sorted.
        for (branch_event_type const & branch_event : _branch_event_queue) // sorted by reference postion
        {
            delta_event_shared_type const * delta_event = branch_event.event_handle();

            if (!delta_event->coverage()[idx]) // Continue if event is not covered by target sequence.
                continue;

            // Record the event in the journal decorator.
            std::visit([&] (auto const & event_kind)
            {
                using substitution_t = typename delta_event_shared_type::substitution_type;
                using insertion_t = typename delta_event_shared_type::insertion_type;
                using deletion_t = typename delta_event_shared_type::deletion_type;

                size_type target_position = delta_event->position() + target_position_offset;

                seqan3::detail::multi_invocable
                {
                    [&] (substitution_t const & e) { jd.record_substitution(target_position, std::span{e.value()}); },
                    [&] (insertion_t const & e) { jd.record_insertion(target_position, std::span{e.value()}); },
                    [&] (deletion_t const & e) { jd.record_deletion(target_position, target_position + e.value()); }
                }(event_kind);
            }, delta_event->delta_variant());

            // Update the target position offset by the insertion and deletion size.
            target_position_offset += delta_event->insertion_size() - static_cast<int32_t>(delta_event->deletion_size());
        }

        // The target sequence and the journal decorator must be equal.
        return std::ranges::equal(jd, _sequences[idx]);
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
 * supports a special cursor to enable an efficient, compression parallel traversal over the stored sequences.
 * This can be used in conjunction with any context based streaming algorithm to speed-up the search against large
 * collection of sequences.
 */
template <seqan3::sequence sequence_t>
using journaled_sequence_tree = typename no_adl::journaled_sequence_tree<sequence_t>::type;
}  // namespace libjst
