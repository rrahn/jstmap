// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides implementation of libjst::journaled_sequence.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <iterator>
#include <limits>
#include <ranges>

#include <libjst/journal/inline_sequence_journal.hpp>

namespace libjst
{

    template <libjst::preserving_reference_sequence source_t>
    class journaled_sequence
    {
        // ----------------------------------------------------------------------------
        // Member classes
        // ----------------------------------------------------------------------------
    public:
        template <bool>
        class iterator_impl;

        // ----------------------------------------------------------------------------
        // Member Types
        // ----------------------------------------------------------------------------
    private:
        using journal_type = inline_sequence_journal<source_t>;
        using source_type = typename journal_type::source_type;
        using journal_breakend_type = journal_type::breakend_type;
        using journal_breakpoint_type = journal_type::breakpoint_type;
        using entry_type = std::ranges::range_value_t<journal_type>;

    public:
        using size_type = std::ranges::range_size_t<journal_type>;
        using sequence_type = typename journal_type::sequence_type;
        using iterator = iterator_impl<false>;
        using const_iterator = iterator_impl<true>;

        // ----------------------------------------------------------------------------
        // Member Variables
        // ----------------------------------------------------------------------------
    private:
        journal_type _journal{};

        // ----------------------------------------------------------------------------
        // Member functions
        // ----------------------------------------------------------------------------
    public:
        // ----------------------------------------------------------------------------
        // Constructors, assignment and destructor

        constexpr journaled_sequence()
            noexcept(std::is_nothrow_default_constructible_v<journal_type>)
            requires std::default_initializable<journal_type>
            : _journal{}
        {
        }

        constexpr journaled_sequence(source_type const & source)
            noexcept(std::is_nothrow_constructible_v<journal_type, source_type const &>)
            : _journal{source}
        {
        }

        constexpr journaled_sequence(source_type && source)
            noexcept(std::is_nothrow_constructible_v<journal_type, source_type &&>)
            : _journal{std::move(source)}
        {
        }

        constexpr source_type const & source() const & noexcept
        {
            return _journal.source();
        }

        constexpr source_type source() && noexcept(std::is_nothrow_move_constructible_v<source_type>)
        {
            return std::move(_journal).source();
        }

        // ----------------------------------------------------------------------------
        // Iterators

        iterator begin() noexcept
        {
            return iterator{std::addressof(_journal), _journal.begin()};
        }

        const_iterator begin() const noexcept
        {
            return const_iterator{std::addressof(_journal), _journal.begin()};
        }

        iterator end() noexcept
        {
            return iterator{std::addressof(_journal), _journal.end()};
        }

        const_iterator end() const noexcept
        {
            return const_iterator{std::addressof(_journal), _journal.end()};
        }

        // ----------------------------------------------------------------------------
        // Capacity

        constexpr bool empty() const noexcept
        {
            return _journal.empty();
        }

        constexpr size_type size() const noexcept
        {
            return _journal.end()->position();
        }

        constexpr size_type max_size() const noexcept
        {
            return std::numeric_limits<size_type>::max();
        }

        // ----------------------------------------------------------------------------
        // Modifiers

        void clear() noexcept
        {
            _journal.clear();
        }

        iterator insert(const_iterator position, sequence_type segment)
        {
            return replace(position, position, std::move(segment));
        }

        iterator erase(const_iterator from)
        {
            return erase(from, from + 1);
        }

        iterator erase(const_iterator from, const_iterator to)
        {
            assert(from <= to);
            return replace(std::move(from), std::move(to), sequence_type{});
        }

        iterator replace(const_iterator from, const_iterator to, sequence_type segment)
        {
            return iterator{std::addressof(_journal),
                            _journal.record(journal_breakpoint_type{std::move(from), std::move(to)}, std::move(segment))};
        }
    };

    // deduction guide
    template <libjst::preserving_reference_sequence source_t>
    journaled_sequence(source_t &&) -> journaled_sequence<std::remove_reference_t<source_t>>;

    template <libjst::preserving_reference_sequence source_t>
    template <bool is_const>
    class journaled_sequence<source_t>::iterator_impl
    {
        friend journaled_sequence;

        template <bool>
        friend class iterator_impl;

        // ----------------------------------------------------------------------------
        // Member Types
        // ----------------------------------------------------------------------------
    private:

        using maybe_const_journal_type = std::conditional_t<is_const, journal_type const, journal_type>;
        using journal_iterator = std::ranges::iterator_t<maybe_const_journal_type>;

        using journal_record_t = std::iter_reference_t<journal_iterator>;
        using sequence_type = decltype(std::declval<journal_record_t>().sequence());
        using sequence_iterator = std::ranges::iterator_t<sequence_type>;

    public:
        using value_type = std::iter_value_t<sequence_iterator>;
        using reference = std::iter_reference_t<sequence_iterator>;
        using difference_type = std::iter_difference_t<journal_iterator>;
        using iterator_category = std::random_access_iterator_tag;
        using pointer = typename std::iterator_traits<sequence_iterator>::pointer;

        // ----------------------------------------------------------------------------
        // Member Variables
        // ----------------------------------------------------------------------------
    private:
        maybe_const_journal_type *_journal{};
        journal_iterator _journal_it{};
        sequence_iterator _sequence_it{};

        // ----------------------------------------------------------------------------
        // Member Functions
        // ----------------------------------------------------------------------------

        // ----------------------------------------------------------------------------
        // Constructors, assignment and destructor
    private:

        explicit iterator_impl(maybe_const_journal_type *journal,
                               journal_iterator journal_it) :
            _journal{journal},
            _journal_it{std::move(journal_it)},
            _sequence_it{std::ranges::begin(_journal_it->sequence())}
        {
        }

    public:

        iterator_impl() = default;

        iterator_impl(iterator_impl<!is_const> other) requires is_const :
            _journal{std::move(other._journal)},
            _journal_it{std::move(other._journal_it)},
            _sequence_it{std::move(other._sequence_it)}
        {
        }

        // ----------------------------------------------------------------------------
        // Element access

        reference operator*() const noexcept
        {
            return *_sequence_it;
        }

        pointer operator->() const noexcept
        {
            return std::addressof(this->operator*());
        }

        reference operator[](difference_type const offset) const noexcept
        {
            return *(*this + offset);
        }
    private:

        constexpr operator journal_breakend_type() const & noexcept
        {
            return journal_breakend_type{_journal_it, _sequence_it};
        }

        constexpr operator journal_breakend_type() && noexcept
        {
            return journal_breakend_type{std::move(_journal_it), std::move(_sequence_it)};
        }

        // ----------------------------------------------------------------------------
        // Arithmetic operators
    public:

        iterator_impl &operator++() noexcept
        {
            assert(_sequence_it < std::ranges::end(_journal_it->sequence()));
            if (++_sequence_it == std::ranges::end(_journal_it->sequence()))
            {
                ++_journal_it;
                _sequence_it = std::ranges::begin(_journal_it->sequence());
            }
            return *this;
        }

        iterator_impl operator++(int) noexcept
        {
            iterator_impl tmp{*this};
            ++(*this);
            return tmp;
        }

        iterator_impl &operator+=(difference_type count) noexcept
        {
            if (auto next_position = current_position() + count; current_record_covers(next_position)) {
                _sequence_it += count;
            } else {
                _journal_it = _journal->lower_bound(next_position);
                if (static_cast<difference_type>(_journal_it->position()) > next_position)
                    --_journal_it;

                assert(static_cast<difference_type>(_journal_it->position()) <= next_position);
                assert(next_position - _journal_it->position() <= std::ranges::size(_journal_it->sequence()));
                _sequence_it = std::ranges::next(std::ranges::begin(_journal_it->sequence()),
                                                 next_position - _journal_it->position());
            }

            return *this;
        }

        friend iterator_impl operator+(iterator_impl const &lhs, difference_type const rhs) noexcept
        {
            iterator_impl tmp{lhs};
            tmp += rhs;
            return tmp;
        }

        friend iterator_impl operator+(difference_type const lhs, iterator_impl const &rhs) noexcept
        {
            return rhs + lhs;
        }

        iterator_impl &operator--() noexcept
        {
            if (_sequence_it == std::ranges::begin(_journal_it->sequence()))
            {
                --_journal_it;
                _sequence_it = std::ranges::end(_journal_it->sequence());
            }
            --_sequence_it;
            return *this;
        }

        iterator_impl operator--(int) noexcept
        {
            iterator_impl tmp{*this};
            --(*this);
            return tmp;
        }

        iterator_impl &operator-=(difference_type const count) noexcept
        {
            return (*this) += -count;
        }

        friend iterator_impl operator-(iterator_impl const &lhs, difference_type const rhs) noexcept
        {
            iterator_impl tmp{lhs};
            tmp -= rhs;
            return tmp;
        }

        friend difference_type operator-(iterator_impl const &lhs, iterator_impl const &rhs) noexcept
        {
            return lhs.current_position() - rhs.current_position();
        }

        // ----------------------------------------------------------------------------
        // Comparison operators

        friend bool operator==(iterator_impl const &lhs, iterator_impl const &rhs) noexcept
        {
            return lhs._journal_it == rhs._journal_it && lhs._sequence_it == rhs._sequence_it;
        }

        friend std::strong_ordering operator<=>(iterator_impl const &lhs, iterator_impl const &rhs) noexcept
        {
            return lhs.current_position() <=> rhs.current_position();
        }

        // ----------------------------------------------------------------------------
        // Utility functions
    private:
        difference_type current_position() const noexcept
        {
            assert(static_cast<difference_type>(_journal_it->position()) >= 0);

            return _journal_it->position() +
                   std::ranges::distance(std::ranges::begin(_journal_it->sequence()), _sequence_it);
        }

        constexpr bool current_record_covers(difference_type next_position) const noexcept
        {
            using position_t = typename std::remove_cvref_t<journal_record_t>::size_type;
            return _journal_it->position() <= static_cast<position_t>(next_position) &&
                   static_cast<position_t>(next_position) < (_journal_it + 1)->position();
        }
    };

} // namespace libjst
