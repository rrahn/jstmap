// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a journal implementation using contiguous memory to store the elements.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <iterator>
#include <vector>
#include <ranges>
#include <type_traits>

#include <libjst/sequence/journal_position.hpp>
#include <libjst/reference_sequence/reference_sequence_concept.hpp>

namespace libjst
{

    template <libjst::preserving_reference_sequence source_t>
    class inline_sequence_journal
    {
    private:
        class record_impl;

        class breakend_impl;

        ///@name Member types
        ///@{
    private:

        using journal_type = std::vector<record_impl>;

    public:
        using source_type = source_t;
        using sequence_type = typename record_impl::sequence_type;
        using breakend_type = breakend_impl;
        using breakpoint_type = std::pair<breakend_type, breakend_type>;
        using size_type = std::ranges::range_size_t<source_type>;
        using key_type = size_type;
        using iterator = std::ranges::iterator_t<journal_type>;
        using const_iterator = std::ranges::iterator_t<journal_type const>;
        ///@}

        /// @name Member Variables
        /// @{
    private:
        source_type _source;
        journal_type _journal{};
        /// @}

        /// @name Member functions
        /// @{
    public:

        constexpr inline_sequence_journal()
            noexcept(noexcept(initialize_journal()))
            requires std::default_initializable<source_type>
            : _journal{}
        {
            initialize_journal();
        }

        constexpr explicit inline_sequence_journal(source_type source)
            noexcept(noexcept(initialize_journal()))
            : _source{std::move(source)},
              _journal{}
        {
            initialize_journal();
        }

        constexpr source_type const & source() const & noexcept
        {
            return _source;
        }

        constexpr source_type source() && noexcept
        {
            return std::move(_source);
        }
    protected:

        constexpr journal_type & journal() noexcept
        {
            return _journal;
        }

        constexpr journal_type const & journal() const noexcept
        {
            return _journal;
        }
        /// @}

        /// @name Iterators
        /// @{
    public:
        iterator begin() noexcept
        {
            return _journal.begin();
        }

        const_iterator begin() const noexcept
        {
            return _journal.begin();
        }

        iterator end() noexcept
        {
            return _journal.end() - 1;
        }

        const_iterator end() const noexcept
        {
            return _journal.end() - 1;
        }
        /// @}

        /// @name Capacity
        /// @{
    public:
        size_type size() const noexcept
        {
            return _journal.size() - 1;
        }

        size_type max_size() const noexcept
        {
            return _journal.max_size() - 1;
        }

        constexpr bool empty() const noexcept
        {
            return size() == 0;
        }
        /// @}

        /// @name Modifiers
        /// @{
    public:
        void clear() noexcept(noexcept(initialize_journal()))
        {
            _journal.clear();
            initialize_journal();
        }

        constexpr iterator record(breakpoint_type breakpoint, sequence_type sequence)
        {
            return record_inline(std::move(breakpoint), std::move(sequence));
        }
        /// @}

        /// @name Lookup
        /// @{
    public:
        iterator lower_bound(key_type const key) noexcept
        {
            return std::ranges::lower_bound(begin(), end(), key, std::ranges::less{}, &record_impl::position);
        }

        const_iterator lower_bound(key_type const key) const noexcept
        {
            return std::ranges::lower_bound(begin(), end(), key, std::ranges::less{}, &record_impl::position);
        }

        iterator upper_bound(key_type const key) noexcept
        {
            return std::ranges::upper_bound(begin(), end(), key, std::ranges::less{}, &record_impl::position);
        }

        const_iterator upper_bound(key_type const key) const noexcept
        {
            return std::ranges::upper_bound(begin(), end(), key, std::ranges::less{}, &record_impl::position);
        }

        iterator find(key_type const key) noexcept
        {
            return std::ranges::find(begin(), end(), key, std::ranges::less{}, &record_impl::position);
        }

        const_iterator find(key_type const key) const noexcept
        {
            return std::ranges::find(begin(), end(), key, std::ranges::less{}, &record_impl::position);
        }
        /// @}

        /// @name Utilities
        /// @{
    private:

        template <typename sequence_t>
        constexpr static auto get_breakpoint_slice(sequence_t && segment,
                                                   std::ranges::iterator_t<sequence_t> from,
                                                   std::ranges::iterator_t<sequence_t> to)
        {
            auto breakpoint = libjst::to_breakpoint(segment, std::move(from), std::move(to));
            return libjst::breakpoint_slice(segment, std::move(breakpoint));
        }

        /**
         * @brief Splits the journal entry at the given position.
         *
         * @param split_it The iterator pointing to the split position.
         * @retval std::pair<record_impl, record_impl> A pair of journal entries, where the first entry covers the
         *         interval `[segment().begin(), split_it)` and the second entry covers the interval
         *         `[split_it, segment().end())`.
         * @rst
         * .. note:: If ``split_it`` is out of bounds, either the first entry or the second entry will have an
         *           empty segment.
         * @endrst
         *
         * @rst
         * .. danger:: If ``split_it`` is not a valid iterator of the sequence segment, the behaviour is undefined.
         * @endrst
         */
        constexpr std::pair<record_impl, record_impl> split_at(breakend_type breakend) const noexcept
        {
            auto [journal_it, split_it] = std::move(breakend).base();
            auto && segment = journal_it->sequence();

            auto record_prefix = record_impl{journal_it->position(),
                                             get_breakpoint_slice(segment, std::ranges::begin(segment), split_it)};

            auto split_offset = std::ranges::distance(std::ranges::begin(segment), std::move(split_it));
            auto record_suffix = record_impl{journal_it->position() + split_offset,
                                             get_breakpoint_slice(segment, split_it, std::ranges::end(segment))};

            return {std::move(record_prefix), std::move(record_suffix)};
        }

        constexpr void force_overwrite_through(std::ranges::iterator_t<journal_type const> it, record_impl record)
            noexcept(std::is_nothrow_move_constructible_v<record_impl>)
        {
            auto out_it = std::ranges::begin(_journal) + std::ranges::distance(_journal.cbegin(), it);
            *out_it = std::move(record);
        }

        iterator record_inline(breakpoint_type breakpoint, sequence_type new_sequence)
        {
            std::ptrdiff_t deletion_size = libjst::breakend_span(breakpoint);
            std::ptrdiff_t insertion_size = std::ranges::ssize(new_sequence);

            std::array<record_impl, 2> entries_marked_for_insertion{};
            size_t marked_insert_entries{};

            auto [low_prefix, low_suffix] = split_at(libjst::low_breakend(breakpoint));
            auto [high_prefix, high_suffix] = split_at(libjst::high_breakend(breakpoint));

            if (!std::ranges::empty(low_prefix.sequence()))
                entries_marked_for_insertion[marked_insert_entries++] = std::move(low_prefix);

            if (insertion_size > 0)
                entries_marked_for_insertion[marked_insert_entries++] =
                        record_impl{low_suffix.position(), std::move(new_sequence)};

            auto from_journal_it = std::get<0>(libjst::low_breakend(breakpoint).base());
            auto to_journal_it = std::get<0>(libjst::high_breakend(breakpoint).base());
            force_overwrite_through(to_journal_it, std::move(high_suffix));
            auto insert_position = _journal.erase(from_journal_it, to_journal_it);
            auto insert_entries_it = std::make_move_iterator(entries_marked_for_insertion.begin());
            auto first_inserted = _journal.insert(insert_position,
                                                  insert_entries_it,
                                                  insert_entries_it + marked_insert_entries);

            update_positions_of_remaining_entries(first_inserted + marked_insert_entries, insertion_size - deletion_size);

            assert(check_journal_invariants());

            return std::ranges::next(first_inserted, marked_insert_entries - (insertion_size > 0));
        }

        void update_positions_of_remaining_entries(journal_type::iterator journal_it, std::ptrdiff_t offset) noexcept
        {
            if (offset == 0)
                return;

            for (; journal_it != std::ranges::end(_journal); ++journal_it)
                *journal_it = record_impl{static_cast<size_type>(journal_it->position() + offset),
                                          std::move(journal_it->sequence())};
        }

        bool check_journal_invariants() const noexcept {
            // 1. check if first record starts at position 0
            // 2. check if all adjacent entries are non-overlapping such that end position of smaller record is equal to begin position of successor
            if (_journal.empty() || _journal.begin()->position() != 0)
                return false;

            for (auto it = begin(); it != end(); ++it)
            {
                auto const &record = *it;
                auto const &next_record = *std::ranges::next(it);
                if (record.position() + std::ranges::size(record.sequence()) != next_record.position())
                    return false;
            }
            return true;
        }

        constexpr void initialize_journal()
            noexcept(std::is_nothrow_constructible_v<record_impl, key_type, sequence_type>)
        {
            if (!std::ranges::empty(source())) {
                auto src_bpt = libjst::to_breakpoint(source(), std::ranges::begin(source()), std::ranges::end(source()));
                _journal.emplace_back(0, libjst::breakpoint_slice(source(), std::move(src_bpt)));
            }

            _journal.emplace_back(std::ranges::size(source()), sequence_type{});
        }
        /// @}
    };

    template <libjst::preserving_reference_sequence source_t>
    class inline_sequence_journal<source_t>::record_impl
    {
        ///@name Public member types
        ///@{
    public:
        ///@brief The breakpoint type to specify the inline positions of the segments.
        using size_type = std::ranges::range_size_t<source_t>;
        ///@brief The sequence type is implemented as std::ranges::subrange.*/
        using sequence_type = libjst::breakpoint_slice_t<source_t const>;

        ///@}

        /** @name Member variables */
        ///@{
    private:
        /** @brief The begin position of the referenced segment in the journaled sequence. */
        size_type _position{};
        /** @brief The referenced segment in the journaled sequence. */
        sequence_type _sequence{};
        ///@}

        ///@name Constructors, assignment and destructor */
        ///@{
    public:
        /**
         * @brief Construct a new journal entry object.
         *
         * The position will be zero initialized and the subrange will be default constructed.
         * This constructor is only available if the sequence type is std::default_initializable.
         */
        constexpr record_impl()
            requires std::default_initializable<sequence_type>
        = default;

        /**
         * @brief Constructs a journal entry from a position and a sequence.
         *
         * @param position The start position of the referenced sequence segment in the journaled sequence, which is also used as the key to sort the records.
         * @param sequence The sequence referenced by this journal entry.
         */
        constexpr record_impl(size_type position, sequence_type sequence) :
            _position{position},
            _sequence{std::move(sequence)}
        {
        }
        ///@}

        /** @name Accessors */
        ///@{
    public:
        /**
         * @brief Returns the begin position of the referenced segment in the journaled sequence.
         *
         * @returns size_type
         */
        constexpr size_type position() const noexcept
        {
            return _position;
        }

        /**
         * @brief Returns the represented sequence slice.
         *
         * @returns sequence_type
         * @exception implementation_defined If the sequence type is nothrow copy constructible this function never throws.
         */
        constexpr sequence_type sequence() const noexcept(std::is_nothrow_copy_constructible_v<sequence_type>)
        {
            return _sequence;
        }
        ///@}

        /// @name Non-member functions
        /// @{
    public:
        /**
         * @brief Equality comparison.
         *
         * @param lhs The left hand side of the comparison.
         * @param rhs The right hand side of the comparison.
         *
         * @result `true` if both entries refer to the same segment, `false` otherwise.
         * @rst
         * .. note::
         *    This operator does not compare the segments lexicographically, but only if they point to the same
         *    memory range in the original sequence.
         * @endrst
         */
        constexpr friend bool operator==(record_impl const &lhs, record_impl const &rhs) noexcept
        {
            return lhs.position() == rhs.position() &&
                   std::ranges::begin(lhs.sequence()) == std::ranges::begin(rhs.sequence()) &&
                   std::ranges::size(lhs.sequence()) == std::ranges::size(rhs.sequence());
        }

        /**
         * @brief Three-way comparison.
         *
         * @param lhs The left hand side of the comparison.
         * @param rhs The right hand side of the comparison.
         *
         * This operator only compares the begin positions of the journal entries to determine their order.
         *
         * @retval std::weak_ordering The result of the comparison:
         * @rst
         * |   ``std::weak_ordering::less`` -- if the begin position of ``lhs`` is less than the begin position of ``rhs``.
         * |   ``std::weak_ordering::greater`` -- if the begin position of ``lhs`` is greater than the begin position of ``rhs``.
         * |   ``std::weak_ordering::equivalent`` -- if the begin position of ``lhs`` is equal to the begin position of ``rhs``.
         * @endrst
         */
        constexpr friend std::weak_ordering operator<=>(record_impl const &lhs, record_impl const &rhs) noexcept
        {
            return lhs.position() <=> rhs.position();
        }
        ///@}
    };

    template <libjst::preserving_reference_sequence source_t>
    class inline_sequence_journal<source_t>::breakend_impl
    {
        /// @name Member types
        /// @{
    private:
        using journal_iterator = std::ranges::iterator_t<journal_type const>;
        using sequence_iterator = std::ranges::iterator_t<sequence_type const>;
        /// @}

        /// @name Member variables
        /// @{
    private:
        journal_iterator _journal_it{};
        sequence_iterator _sequence_it{};
        /// @}

        /// @name Member functions
        /// @{
    public:
        constexpr breakend_impl()
            noexcept(std::is_nothrow_default_constructible_v<journal_iterator> &&
                     std::is_nothrow_default_constructible_v<sequence_iterator>) = default;

        constexpr explicit breakend_impl(journal_iterator journal_it, sequence_iterator sequence_it)
            noexcept(std::is_nothrow_move_constructible_v<journal_iterator> &&
                     std::is_nothrow_move_constructible_v<sequence_iterator>) :
            _journal_it{std::move(journal_it)},
            _sequence_it{std::move(sequence_it)}
        {
        }

        constexpr std::pair<journal_iterator, sequence_iterator> base() const &
            noexcept(std::is_nothrow_copy_constructible_v<journal_iterator> &&
                     std::is_nothrow_copy_constructible_v<sequence_iterator>)
        {
            return {_journal_it, _sequence_it};
        }

        constexpr std::pair<journal_iterator, sequence_iterator> base() && noexcept
        {
            return {std::move(_journal_it), std::move(_sequence_it)};
        }
        /// @}

        /// @name Conversion
        /// @{
    public:
        template <std::integral integral_t>
        constexpr operator integral_t() const noexcept
        {
            // how can we guarantee that this always works, although journal_it may point to end of journal?

            return _journal_it->position() +
                   std::ranges::distance(std::ranges::begin(_journal_it->sequence()), _sequence_it);
        }
        /// @}

        /// @name Non-member functions
        /// @{
    public:
        constexpr friend std::ptrdiff_t operator-(breakend_impl const & lhs, breakend_impl const & rhs) noexcept
        {
            return static_cast<std::ptrdiff_t>(lhs) - static_cast<std::ptrdiff_t>(rhs);
        }
        constexpr friend bool operator==(breakend_impl const &, breakend_impl const &) = default;
        constexpr friend std::strong_ordering operator<=>(breakend_impl const &, breakend_impl const &) = default;
        /// @}
    };

}  // namespace libjst
