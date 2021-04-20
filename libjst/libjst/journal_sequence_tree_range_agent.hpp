// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal_sequence_tree_range_agent.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <numeric>
#include <ranges>

#include <seqan3/range/views/zip.hpp>

#include <libjst/context_position.hpp>
#include <libjst/detail/journal_sequence_tree_traverser.hpp>
#include <libjst/search/stack_notification_registry.hpp>
#include <libjst/utility/logger.hpp>

namespace libjst::detail
{

/*!\brief A range over a libjst::journaled_sequence_tree with an integrated libjst::search_stack_notification_registry.
 * \implements std::ranges::input_range
 *
 * \tparam jst_t The type of the libjst::journaled_sequence_tree.
 *
 * \details
 *
 * This agent provides a range interface to the algorithms. During the traversal a stack is used to keep track of which
 * branch is currently visited inside of the algorithm. To allow external algorithms to keep track of these state
 * changes they can attach a libjst::search_stack_observer during the construction.
 * Subsequently, they will be notified whenever a state changes is applied during the traversal.
 */
template <typename jst_t>
class journal_sequence_tree_range_agent :
    public search_stack_notification_registry,
    protected journal_sequence_tree_traverser<journal_sequence_tree_range_agent<jst_t>, jst_t>
{
private:
    //!\brief The base traversal type.
    using base_t = journal_sequence_tree_traverser<journal_sequence_tree_range_agent<jst_t>, jst_t>;

    //!\brief Grant access to the protected and private data members.
    friend base_t;

    // Imported types.
    using typename base_t::coverage_type;
    using typename base_t::size_type;
    using typename base_t::segment_type;

    // The iterator type.
    class iterator;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_range_agent() = default; //!< Default.
    constexpr journal_sequence_tree_range_agent(journal_sequence_tree_range_agent const &) = default;
        //!< Default.
    constexpr journal_sequence_tree_range_agent(journal_sequence_tree_range_agent &&) = default;
        //!< Default.
    constexpr journal_sequence_tree_range_agent & operator=(journal_sequence_tree_range_agent const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_range_agent & operator=(journal_sequence_tree_range_agent &&)
        = default; //!< Default.
    ~journal_sequence_tree_range_agent() = default; //!< Default.

    /*!\brief Constructs the range agent for a given libjst::journaled_sequence_tree and a context size.
     *
     * \tparam observer_t A template parameter pack over all observers to attach; must model
     *                    libjst::search_stack_observer.
     *
     * \param[in] jst A pointer to a const journaled sequence tree.
     * \param[in] context_size The context size to use for the cursor.
     * \param[in] observer The observers to attach to the stack notification registry.
     *
     * \details
     *
     * The range agent is initialised with the given context size and attaches the observer to the stack registry.
     */
    template <search_stack_observer ...observer_t>
    journal_sequence_tree_range_agent(jst_t const * jst, size_t const context_size, observer_t & ...observer) noexcept :
        search_stack_notification_registry{observer...},
        base_t{jst, context_size}
    {}
    //!\}

    /*!\name Iterator
     * \{
     */
    //!\brief Returns an input iterator to the begin of the range.
    iterator begin()
    {
        return iterator{this};
    }

    //!\brief Const qualified iteration is not allowed.
    iterator begin() const = delete;

    //!\brief Returns the sentinel denoting the end of the range.
    std::default_sentinel_t end() noexcept
    {
        return std::default_sentinel;
    }

    //!\brief Const qualified iteration is not allowed.
    std::default_sentinel_t end() const = delete;
    //!\}
};

/*!\brief The iterator of the range agent.
 * \implements std::input_iterator
 *
 * \details
 *
 * This iterator models an input iterator interface to iterate over each element of the journal sequence tree.
 * This iterator is move only, i.e. to generate multiple iterator from the same journal sequence tree, one needs to
 * get different instances of the range agent. The iterator does not model std::output_iterator, i.e. the referenced
 * value cannot be modified by the caller.
 */
template <typename jst_t>
class journal_sequence_tree_range_agent<jst_t>::iterator
{
public:
    /*!\name Associated types
     * \{
     */
    using context_positions_type = std::vector<libjst::context_position>; //!< The context positions type.
    using value_type = std::ranges::range_value_t<segment_type>; //!< The value type of the iterator.
    using reference = value_type; //!< The reference type which is not assignable.
    using pointer = void; //!< No pointer available.
    using difference_type = std::ptrdiff_t; //!< The difference type.
    using iterator_category = std::input_iterator_tag; //!< The iterator category.
    //!\}

private:
    //!\brief The pointer to the underlying host.
    journal_sequence_tree_range_agent * _host{};
    //!\brief The context positions per context.
    context_positions_type _context_positions{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr iterator() = default; //!< Default.
    constexpr iterator(iterator const &) = delete; //!< Deleted.
    constexpr iterator(iterator &&) = default; //!< Default.
    constexpr iterator & operator=(iterator const &) = delete; //!< Deleted.
    constexpr iterator & operator=(iterator &&) = default; //!< Default.
    ~iterator() = default; //!< Default.

    //!\brief Constructs an instance of this iterator from the given host.
    explicit constexpr iterator(journal_sequence_tree_range_agent * host) : _host{host}
    {
        if (!_host->at_end())
            ++(*this);
    }
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the current context.
    reference operator*() const noexcept(noexcept(_host->current_value()))
    {
        return _host->current_value();
    }

    //!\brief Returns a vector with the positions valid for the current context. Can be empty.
    context_positions_type const & positions() noexcept
    {
        assert(_host != nullptr);

        _context_positions.clear();

        coverage_type branch_coverage = _host->determine_supported_context_coverage();
        size_type context_position = _host->context_begin_position();
        size_type sequence_id{};
        for(auto && [offset, is_covered] : seqan3::views::zip(_host->_sequence_offsets, branch_coverage))
        {
            if (is_covered)
                _context_positions.emplace_back(sequence_id, offset + context_position);

            ++sequence_id;
        };

        return _context_positions;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    //!\brief Advances to the next context.
    iterator & operator++() noexcept
    {
        _host->advance();
        return *this;
    }

    //!\brief Advances to the next context.
    iterator operator++(int) noexcept
    {
        iterator tmp{*this};
        ++(*this);
        return tmp;
    }
    //!\}

    /*!\name Comparison operators
     * \{
     */
    //!\brief Compares with the sentinel.
    bool operator==(std::default_sentinel_t const &) const noexcept
    {
        return _host->at_end();
    }
    //!\}
};

} // namespace libjst::detail
