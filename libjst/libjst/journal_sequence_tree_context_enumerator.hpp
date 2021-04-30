// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides libjst::journal_sequence_tree_context_enumerator.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <iterator>
#include <numeric>
#include <ranges>

#include <seqan3/range/views/zip.hpp>
#include <seqan3/utility/detail/multi_invocable.hpp>

#include <libjst/context_position.hpp>
#include <libjst/detail/journal_sequence_tree_traverser.hpp>
#include <libjst/utility/logger.hpp>

namespace libjst::detail
{

/*!\brief A context enumerator for a libjst::journaled_sequence_tree
 * \implements std::ranges::input_range
 *
 * \tparam jst_t The type of the libjst::journaled_sequence_tree.
 *
 * \details
 *
 * This class provides a range interface over all contexts generated by the underlying journal sequence tree.
 * The enumerator implements a single-pass input range over the all possible contexts, given a context size.
 * Multiple context enumerators can be created for the same journal sequence tree and can be processed in parallel.
 * Note, that the interfaces of a single enumerator are not thread-safe unless they are marked const.
 */
template <typename jst_t>
class journal_sequence_tree_context_enumerator :
    protected journal_sequence_tree_traverser<journal_sequence_tree_context_enumerator<jst_t>, jst_t>
{
private:
    //!\brief The base traversal type.
    using base_t = journal_sequence_tree_traverser<journal_sequence_tree_context_enumerator<jst_t>, jst_t>;
    //!\brief The model type.
    using model_t = typename base_t::model_t;

    //!\brief Grant access to the underlying notification methods.
    friend base_t;

    // Imported types.
    using typename base_t::coverage_type;
    using typename base_t::size_type;
    using typename base_t::sequence_context_type;

    // The iterator type.
    class iterator;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr journal_sequence_tree_context_enumerator() = default; //!< Default.
    constexpr journal_sequence_tree_context_enumerator(journal_sequence_tree_context_enumerator const &) = default;
        //!< Default.
    constexpr journal_sequence_tree_context_enumerator(journal_sequence_tree_context_enumerator &&) = default;
        //!< Default.
    constexpr journal_sequence_tree_context_enumerator & operator=(journal_sequence_tree_context_enumerator const &)
        = default; //!< Default.
    constexpr journal_sequence_tree_context_enumerator & operator=(journal_sequence_tree_context_enumerator &&)
        = default; //!< Default.
    ~journal_sequence_tree_context_enumerator() = default; //!< Default.

    /*!\brief Constructs the context enumerator for a given libjst::journaled_sequence_tree and a context
     *        size.
     *
     * \param[in] jst A pointer to a const journaled sequence tree.
     * \param[in] context_size The context size to use for the cursor.
     * \param[in] begin_pos The begin position of the reference to consider.
     * \param[in] end_pos The end position of the reference to consider.
     */
    journal_sequence_tree_context_enumerator(jst_t const * jst,
                                             size_t const context_size,
                                             std::ptrdiff_t begin_pos = 0,
                                             std::ptrdiff_t end_pos = std::numeric_limits<std::ptrdiff_t>::max())
        noexcept :
            base_t{jst, context_size, begin_pos, end_pos}
    {}

    /*!\brief Constructs the context enumerator from a given traverser model and a context size.
     *
     * \param[in] model The model to construct the traverser for.
     * \param[in] context_size The context size to use for the cursor.
     */
    journal_sequence_tree_context_enumerator(model_t model, size_t const context_size) noexcept :
        base_t{std::move(model), context_size}
    {}
    //!\}

    /*!\name Iterator
     * \{
     */
    //!\brief Returns an input iterator to the begin of the context enumerator.
    iterator begin()
    {
        return iterator{this};
    }

    //!\brief Const qualified iteration is not allowed.
    iterator begin() const = delete;

    //!\brief Returns the sentinel denoting the end of the context enumerator.
    std::default_sentinel_t end() noexcept
    {
        return std::default_sentinel;
    }

    //!\brief Const qualified iteration is not allowed.
    std::default_sentinel_t end() const = delete;
    //!\}

private:
    //!\brief NOOP function which does nothing.
    void notify_push() const noexcept
    {}

    //!\brief NOOP function which does nothing.
    void notify_pop() const noexcept
    {}
};

/*!\brief The iterator of the context enumerator.
 * \implements std::input_iterator
 *
 * \details
 *
 * This iterator models an input iterator interface to enumerate all contexts.
 * This iterator is move only, i.e. to generate multiple iterator from the same journal sequence tree, one needs to
 * get different instances of the context enumerator. This is not an output iterator, i.e. the referenced context
 * cannot be modified by the caller.
 */
template <typename jst_t>
class journal_sequence_tree_context_enumerator<jst_t>::iterator
{
public:
    /*!\name Associated types
     * \{
     */
    using context_positions_type = std::vector<libjst::context_position>; //!< The context positions type.
    using value_type = sequence_context_type; //!< The value type is a range over the context.
    using reference = value_type; //!< The reference type is also a range over the context and is not modifiable.
    using pointer = void; //!< No pointer available.
    using difference_type = std::ptrdiff_t; //!< The difference type.
    using iterator_category = std::input_iterator_tag; //!< The iterator category.
    //!\}

private:
    //!\brief The pointer to the underlying host.
    journal_sequence_tree_context_enumerator * _host{};
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
    explicit constexpr iterator(journal_sequence_tree_context_enumerator * host) : _host{host}
    {
        if (!_host->at_end() && !_host->has_full_context_in_branch())
            ++(*this);
    }
    //!\}

    /*!\name Element access
     * \{
     */
    //!\brief Returns the current context.
    reference operator*() const noexcept(noexcept(_host->current_context()))
    {
        return _host->current_context();
    }

    //!\brief Returns a vector with the positions valid for the current context. Can be empty.
    context_positions_type const & positions() noexcept
    {
        assert(_host != nullptr);

        _context_positions.clear();

        coverage_type branch_coverage = _host->determine_supported_context_coverage();
        size_type context_position = _host->context_begin_position();
        size_type sequence_id{};
        for (auto && [offset, is_covered] : seqan3::views::zip(_host->_sequence_offsets, branch_coverage))
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
        while (!_host->next_context())
        {}

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
